import argparse
import yaml
import numpy as np
import ROOT
import ctypes
import uproot
import pandas as pd
import array
import os
import re
from ROOT import TFile, TKDE, TCanvas, TH1D, TF1

def templ_producer_kde(tree, pt_min, pt_max, mass_min, mass_max, name, outfile='', bkg_min=0, bkg_max=1, fd_min=0, fd_max=1, var='fM'):
    
    if fd_max != 0:
        query = f"{pt_min} <= fPt < {pt_max} and fMlScore0 < {bkg_max} and fMlScore1 >= {fd_min} and fMlScore1 < {fd_max}"
        # query = f"{mass_min} <= fM and fM < {mass_max} and {pt_min} <= fPt < {pt_max} and fMlScore0 < {bkg_max} and fMlScore1 >= {fd_min} and fMlScore1 < {fd_max}"
    else:
        query = f"{pt_min} <= fPt and fPt < {pt_max}"
        # query = f"{mass_min} <= fM and fM < {mass_max} and {pt_min} <= fPt and fPt < {pt_max}"

    print(f"Producing KDE from {name} for var {var}, query: {query}")
    var_values = tree.query(query)[var].tolist()
    # print(f"var_values: {var_values}")
    # quit()
    if len(var_values) > 10:
        kde = TKDE(len(var_values), np.asarray(var_values, 'd'), mass_min, mass_max)
        
        binned_var_values = TH1D(f'hBinned', f'hBinned', 2000, mass_min, mass_max)
        for var_value in var_values:
            binned_var_values.Fill(var_value)

        if outfile != '':
            max_content = 0
            kde_func = kde.GetFunction(500)
            print(f"kde_func.Integral(1,3): {kde_func.Integral(1,3)}")
            for bin_idx in range(1, binned_var_values.GetNbinsX() + 1):
                bin_content = binned_var_values.GetBinContent(bin_idx)
                if bin_content > max_content:
                    max_content = bin_content
                    max_bin = bin_idx
            binned_var_values.Scale(1 / binned_var_values.Integral(), 'width')
            # binned_var_values.Scale(kde_func.GetMaximum() / binned_var_values.GetBinContent(max_bin))

            cOverlap = TCanvas('cOverlap', 'cOverlap', 600, 600)
            cOverlap.cd()
            binned_var_values.Draw()
            kde_func.Draw('same')
            outfile.mkdir(f'KDE_pT_{pt_min}_{pt_max}_{name}')
            outfile.cd(f'KDE_pT_{pt_min}_{pt_max}_{name}')
            kde.Write('kde')
            binned_var_values.Write()
            kde_func.Write()
            cOverlap.Write()
        # quit()
        
        return kde
        # return kde, kde_func, binned_var_values
    else:
        print('CIAO ELSE')
        # quit()
        return None
    
def templ_producer_histo(tree_file, var, pt_min, pt_max, queries, names, relweights=[], outfile='', tree_name='O2hfcanddplite'):

    print(f"Producing KDE from {tree_file} for var {var}, {pt_min} <= pt < {pt_max}, names {names}")
    # convert the tree_file to a pandas dataframe
    dfsData = []
    print(f"tree_file: {tree_file}")
    with uproot.open(f'{tree_file}') as f:
        for key in f.keys():
            if tree_name in key:
                dfData = f[key].arrays(library='pd')
                dfsData.append(dfData)      
    df = pd.concat([df for df in dfsData], ignore_index=True)
    histos_templ = []
    for query, name in zip(queries, names):
        print(f"query: {query}")
        print(f"{pt_min} < fPt < {pt_max} and {query}")
        templ_df = df.query(f"{pt_min} < fPt < {pt_max} and {query}")[var].to_numpy()
        histos_templ.append(ROOT.TH1D(
            f"hist_templ_{name}_pt{pt_min:.1f}_{pt_max:.1f}",
            "#it{M}(K#pi#pi) (GeV/#it{c})", 1500, 1.00, 2.50))
        for var_value in templ_df:
            histos_templ[-1].Fill(var_value)

    histo_comb = ROOT.TH1D(
        f"hist_templ_combined_pt{pt_min:.1f}_{pt_max:.1f}",
        "#it{M}(K#pi#pi) (GeV/#it{c})", 1500, 1.00, 2.50)

    if relweights != []:
        for irelweight, histo_templ in zip(relweights, histos_templ):
            histo_comb.Add(histo_templ, irelweight)
    else:
        for irelweight, histo_templ in zip(relweights, histos_templ):
            histo_comb.Add(histo_templ, 1)

    histo_comb_smoothened = histo_comb.Clone(f"{histo_comb.GetName()}_smooth")
    histo_comb_smoothened.Smooth(100)

    if outfile != '':
        outfile.mkdir(f'hTempl_pT_{pt_min}_{pt_max}')
        outfile.cd(f'hTempl_pT_{pt_min}_{pt_max}')
        for hist in histos_templ:
            hist.Write()
        histo_comb.Write()
        histo_comb_smoothened.Write()

    return histo_comb

def extract_template_weights(config, correlated):

    with open(config, 'r') as cfg:
        config = yaml.safe_load(cfg)

    os.makedirs(f"{config['out_dir']}/cutvar_{config['suffix']}/ry/", exist_ok=True)
    weights_file = TFile(f"{config['out_dir']}/cutvar_{config['suffix']}/ry/weights.root", 'recreate')

    ###### MC decay tables
    # D+ decay table from https://github.com/AliceO2Group/O2DPG/blob/master/MC/config/PWGHF/pythia8/generator/pythia8_charmhadronic_with_decays_Mode2.cfg
    # 411:oneChannel = 1 0.0752 0 -321 211 211
    # 411:addChannel = 1 0.0104 0 -313 211
    # 411:addChannel = 1 0.0156 0 311 211
    # 411:addChannel = 1 0.0752 0 333 211, same amount of D+->KKpi and D+->Kpipi
    
    # Ds decay table in MC --> all in Ds --> KKpi
        
    ##### PDG branching ratios
    # D+ -> Kpipi: 9.38e-2
    # D+ -> KKpi: 9.68e-3
    # Ds+ -> KKpi: 5.37e-2
        
    # Reweight contributions with (BR_PDG / BR_MC)
    BRDplusTotMC = 0.0752 + 0.0104 + 0.0156 + 0.0752
    BRDplusKPiPiMC = 0.0752 + 0.0156 + 0.0104
    BRDplusKKPiMC = 0.0752
    BRDplusKPiPiPDG = 9.38e-2
    BRDplusKKPiPDG = 9.68e-3
    BRDsKKPiPDG = 5.37e-2
    BRDsKKPiMC = 1.
    BRDsPhiPiKKPiPDG = 2.21e-2 + 0.29e-2  # 0.29 is added so that the BRs of Ds --> Phi Pi and Ds --> K*K add to 5.37e-2
    BRDsPhiPiKKPiMC = 4.4e-2 / (4.e-2 + 4.4e-2)   # Same BR of native pythia is kept, but only Ds -> Phi Pi 
                                                  # and Ds -> K*K are implemented, thus they are rescaled
    BRDsKStarKKKPiPDG = 2.58e-2 + 0.29e-2
    BRDsKStarKKKPiMC = 1 - BRDsPhiPiKKPiMC
    BRD0KPiPDG = 3.89e-2
    BRD0KPiMC = 3.89e-2
    BRD0KKPiMC = 3.89e-3
    
    signalBRNorm = BRDplusKPiPiPDG / (BRDplusKPiPiMC / BRDplusTotMC)
    templatesBRNorms = []
    for _, templ in enumerate(config['TemplsNames']):
        if templ == "DsKKPi":
            # Ds/D+ is underestimated in pythia CRMode2 --> multiply by 1.25
            templatesBRNorms.append( (BRDsKKPiPDG / BRDsKKPiMC) * 1.25)
        elif templ == "DsPhiPi":
            # Ds/D+ is underestimated in pythia CRMode2 --> multiply by 1.25
            templatesBRNorms.append( (BRDsPhiPiKKPiPDG / BRDsPhiPiKKPiMC) * 1.25)
        elif templ == "DsKStarK":
            # Ds/D+ is underestimated in pythia CRMode2 --> multiply by 1.25
            templatesBRNorms.append( (BRDsKStarKKKPiPDG / BRDsKStarKKKPiMC) * 1.25)
        elif templ == "DplusKKPi":
            templatesBRNorms.append(BRDplusKKPiPDG / (BRDplusKKPiMC / BRDplusTotMC))
        elif templ == "DstarKPiPi":
            # the decay table of D* is not modified in the MC, thus BR_PDG / BR_MC = 1
            # but the modification of the decay table of D0 needs to be taken into account
            templatesBRNorms.append(1. * (BRD0KPiPDG / (BRD0KPiMC / (BRD0KPiPDG + BRD0KKPiMC))))
        else:
            templatesBRNorms.append(1.)

    ### load the trees for bkg and signal
    templatesYieldsDfs = [pd.read_parquet(config['SignalPath'])] + [pd.read_parquet(templPath) for templPath in config['TemplsPaths']]
    templatesYieldsNames = ["Signal"] + config['TemplsNames']
    templatesBRNorms = [signalBRNorm] + templatesBRNorms
    
    ### Loop over the cutsets
    config_files = [f for f in os.listdir(f"{config['out_dir']}/cutvar_{config['suffix']}/config/") if os.path.isfile(os.path.join(f"{config['out_dir']}/cutvar_{config['suffix']}/config/", f))]
    for config_file in config_files:
        match = re.search(r"(\d+)\.yml$", config_file)
        if not match:
            continue
        cutset = match.group(1)

        with open(f"{config['out_dir']}/cutvar_{config['suffix']}/config/{config_file}", 'r') as cfg:
            config_cut = yaml.safe_load(cfg)

        pt_mins = config_cut['cutvars']['Pt']['min']
        pt_maxs = config_cut['cutvars']['Pt']['max']
        pt_bins = array.array('d', pt_mins + [pt_maxs[-1]])

        # Initialize histograms
        histos = {
            name: {
                "mass": None,
                "massbr": None,
                "massbryield": None,
                "yield": TH1D(f"reco_yields_{name}", f";p_T;{name}/Counts", len(pt_bins)-1, pt_bins),
                "br": TH1D(f"br_{name}", f";p_T;BR", len(pt_bins)-1, pt_bins),
                "rew": TH1D(f"rew_factor_{name}", f";p_T;BR x N_{{reco}}", len(pt_bins)-1, pt_bins),
                "rel": TH1D(f"rel_weight_to_sgn_{name}", f";p_T;Rel. Sgn. Weight", len(pt_bins)-1, pt_bins)
            } for name in templatesYieldsNames
        }

        for iPt, (ptmin, ptmax) in enumerate(zip(pt_mins, pt_maxs)):
            pt_dir = f"cutset_{cutset}/pt_{ptmin}_{ptmax}/"
            weights_file.mkdir(pt_dir)
            weights_file.cd(pt_dir)

            fit_min = config_cut['fitrangemin'][iPt]
            fit_max = config_cut['fitrangemax'][iPt]
            nbins = int((fit_max - fit_min) * 1000)

            histo_total_br = TH1D(f"histo_total_br", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, fit_min, fit_max)
            histo_total_br_yield = TH1D(f"histo_total_br_yield", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, fit_min, fit_max)
            # Evaluate templates
            for name, df, br in zip(templatesYieldsNames, templatesYieldsDfs, templatesBRNorms):
                df_sel = df.query(f"{ptmin} <= fPt < {ptmax} and {fit_min} <= fM < {fit_max}")

                if config.get('MlDiffWeights'):
                    score_bkg = config_cut['cutvars'].get('score_bkg')
                    score_fd = config_cut['cutvars'].get('score_FD')
                    if score_bkg:
                        df_sel = df_sel.query(f"fMlScore0 < {score_bkg['max'][iPt]}")
                    if score_fd:
                        if correlated:
                            df_sel = df_sel.query(f"fMlScore1 >= {score_fd['min'][iPt]}")
                        else:
                            df_sel = df_sel.query(f"{score_fd['min'][iPt]} <= fMlScore1 < {score_fd['max'][iPt]}")

                n_reco = len(df_sel)
                histos[name]["yield"].SetBinContent(iPt+1, n_reco)
                histos[name]["br"].SetBinContent(iPt+1, br)
                histos[name]["rew"].SetBinContent(iPt+1, n_reco * br)
                histos[name]["rel"].SetBinContent(iPt+1, histos[name]["rew"].GetBinContent(iPt+1) / histos["Signal"]["rew"].GetBinContent(iPt+1))
                print(f"n_reco: {n_reco}, br: {br}")

                # Fill and write mass histogram
                histos[name]["mass"] = TH1D(f"htempl_{name}", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, fit_min, fit_max)
                histos[name]["massbr"] = TH1D(f"htempl_{name}", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, fit_min, fit_max)
                histos[name]["massbryield"] = TH1D(f"htempl_{name}", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, fit_min, fit_max)
                for mass in df_sel["fM"].to_numpy():
                    histos[name]["mass"].Fill(mass)
                    histos[name]["massbr"].Fill(mass)
                    histos[name]["massbryield"].Fill(mass)
                histos[name]["mass"].Smooth(10)
                histos[name]["mass"].Write(f"h{name}Raw")
                histos[name]["massbr"].Smooth(10)
                histos[name]["massbr"].Scale(histos[name]["br"].GetBinContent(iPt+1) / histos["Signal"]["br"].GetBinContent(iPt+1))
                histos[name]["massbryield"].Smooth(10)
                histos[name]["massbryield"].Scale(histos[name]["rel"].GetBinContent(iPt+1))

                if name != "Signal":
                    histo_total_br.Add(histos[name]["massbr"])
                    histo_total_br_yield.Add(histos[name]["massbryield"])

            histo_total_br.Write(f"histo_total_br_rew")
            histo_total_br_yield.Write(f"histo_total_br_yield_rew")
                
        # Save summary histograms
        for name, histset in histos.items():
            summary_dir = f"cutset_{cutset}/{name}/"
            weights_file.mkdir(summary_dir)
            weights_file.cd(summary_dir)
            histset["mass"].Write("hMassRaw")
            histset["massbr"].Write("hMassBrWrtSgn")
            histset["massbryield"].Write("hMassBrYieldWrtSgn")
            histset["yield"].Write("hYield")
            histset["br"].Write("hBR")
            histset["rew"].Write("hRewFactor")
            histset["rel"].Write("hRelWeightToSgn")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(descriPtion="Arguments")
    parser.add_argument("--config", "-cfg", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("--var", "-v", metavar="text",
                        default="fM", help="variable of interest")
    parser.add_argument("--ptmin", "-pmin", metavar="text",
                        default="2.", help="min pt")
    parser.add_argument("--ptmax", "-pmax", metavar="text",
                        default="4.", help="max pt")
    parser.add_argument("--flag", "-f", metavar="chn flag",
                        default="2", help="channel flag")
    parser.add_argument("--input", "-in", metavar="path/input.root",
                        default="AnalysisResults.root", help="path to file containing tree")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--correlated", "-c", action="store_true", 
                        help="Produce yml files for correlated cuts")
    
    args = parser.parse_args()
    
    KDEs = []
    binned_histos = []
    if args.config != parser.get_default("config"):
        with open(args.config, 'r') as ymlCfgFile:
            config = yaml.load(ymlCfgFile, yaml.FullLoader)
        
    output_dir = config["outputdir"] if args.outputdir == parser.get_default("outputdir") else args.outputdir
    suffix = config["suffix"] if args.suffix == parser.get_default("suffix") else args.suffix
    outfile = ROOT.TFile(f'{output_dir}/kde_{suffix}.root', 'RECREATE')
        
    if args.config != parser.get_default("config"):
        for pt_low, pt_max in zip(config['pt_mins'], config['pt_maxs']):
            KDE, _, histo = templ_producer_kde(config['input'], config['variable'], pt_low, 
                                           pt_max, config['chn_flag'], outfile)
    else:
        KDE, _, histo = templ_producer_kde(args.input, args.var, args.ptmin,
                                       args.ptmax, args.flag, outfile)
    outfile.Close()    
