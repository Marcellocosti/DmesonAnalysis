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
        kde = TKDE(len(var_values), np.asarray(var_values, 'd'), 1, 3)
        
        binned_var_values = TH1D(f'hBinned', f'hBinned', 2000, 1, 3)
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

def extract_template_weights(config):

    with open(config, 'r') as cfg:
        config = yaml.safe_load(cfg)

    weights_file = TFile(config['WeightsFile'], 'recreate')

    templatesBRNorms = []
    if config['Dmeson'] == 'Dplus':
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
        BRD0KPiPDG = 3.89e-2
        BRD0KPiMC = 3.89e-2
        BRD0KKPiMC = 3.89e-3
        
        print(config['TemplsNames'])
        for iTemplate, templ in enumerate(config['TemplsNames']):
            if templ == "DsKKPi":
                # Ds/D+ is underestimated in pythia CRMode2 --> multiply by 1.25
                templatesBRNorms.append( (BRDsKKPiPDG / BRDsKKPiMC) * 1.25)
                if iTemplate == 0:
                    firstTemplBRNorm = (BRDsKKPiPDG / BRDsKKPiMC) * 1.25
            elif templ == "DplusKKPi":
                templatesBRNorms.append(BRDplusKKPiPDG / (BRDplusKKPiMC / BRDplusTotMC))
                if iTemplate == 0:
                    firstTemplBRNorm = BRDplusKKPiPDG / (BRDplusKKPiMC / BRDplusTotMC)
            elif templ == "DstarKPiPi":
                # the decay table of D* is not modified in the MC, thus BR_PDG / BR_MC = 1
                # but the modification of the decay table of D0 needs to be taken into account
                templatesBRNorms.append(1. * (BRD0KPiPDG / (BRD0KPiMC / (BRD0KPiPDG + BRD0KKPiMC))))
            else:
                templatesBRNorms.append(1.)
        signalBRNorm = BRDplusKPiPiPDG / (BRDplusKPiPiMC / BRDplusTotMC)

    ### load the trees for bkg and signal
    templatesDfs = [pd.read_parquet(templPath) for templPath in config['TemplsPaths']]
    signalDf = pd.read_parquet(config['SignalPath'])

    templatesYieldsDfs = [signalDf] + templatesDfs 
    templatesYieldsNames = ["Signal"] + config['TemplsNames']
    templatesBRNorms = [signalBRNorm] + templatesBRNorms
    config_files = [f for f in os.listdir(f"{config['out_dir']}/cutvar_{config['suffix']}/config/") if os.path.isfile(os.path.join(f"{config['out_dir']}/cutvar_{config['suffix']}/config/", f))]
    print(f"config_files: {config_files}")
    
    ### Loop over the cutsets
    for config_file in config_files:
        match = re.search(r"(\d+)\.yml$", os.path.basename(config_file))
        if match:
            cutset = match.group(1)
        with open(f"{config['out_dir']}/cutvar_{config['suffix']}/config/{config_file}", 'r') as cfg:
            config_cut = yaml.safe_load(cfg)
        pt_bins = array.array('d', config_cut['cutvars']['Pt']['min'] + [config_cut['cutvars']['Pt']['max'][-1]])
        hist_frac_templ_to_signal = {}
        hist_frac_templ_to_firsttempl = {}
        for iTempl, templName in enumerate(templatesYieldsNames):
            hist_frac_templ_to_signal[templName] = ROOT.TH1D(f"hist_{templName}_over_signal", f";p_T;{templName}/Signal", len(pt_bins)-1, pt_bins)
            if iTempl >= 2:
                hist_frac_templ_to_firsttempl[templName] = ROOT.TH1D(f"hist_{templName}_over_{templatesYieldsNames[1]}", f";p_T;{templName}/{templatesYieldsNames[1]}", len(pt_bins)-1, pt_bins)
        hist_frac_templ_to_signal["TotalBkg"] = ROOT.TH1D(f"hist_totalbkg_over_signal", f";p_T;Bkg/Signal", len(pt_bins)-1, pt_bins)
        print(f"hist_frac_templ_to_signal: {hist_frac_templ_to_signal}")

        ### Obtain the raw histo template and the one reweighted with the respective BR
        for iPt, (ptmin, ptmax) in enumerate(zip(config_cut['cutvars']['Pt']['min'], config_cut['cutvars']['Pt']['max'])):
            nbins = int( (config_cut['fitrangemax'][iPt] - config_cut['fitrangemin'][iPt]) * 1000)
            hist_templ_total = ROOT.TH1D(f"hist_templ_total", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, config_cut['fitrangemin'][iPt], config_cut['fitrangemax'][iPt])
            hist_signal = ROOT.TH1D(f"hist_signal", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, config_cut['fitrangemin'][iPt], config_cut['fitrangemax'][iPt])
            hist_first_templ = ROOT.TH1D(f"hist_first_templ", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, config_cut['fitrangemin'][iPt], config_cut['fitrangemax'][iPt])
            for iTemplate, (templName, templDf, BRnorm) in enumerate(zip(templatesYieldsNames, templatesYieldsDfs, templatesBRNorms)):
                weights_file.mkdir(f"cutset_{cutset}/{templName}/pt_{ptmin}_{ptmax}/")
                weights_file.cd(f"cutset_{cutset}/{templName}/pt_{ptmin}_{ptmax}/")
                templDfPt = templDf.query(f"fPt >= {ptmin} and fPt < {ptmax}")
                if config.get('MlDiffWeights'):
                    if config_cut['cutvars'].get('score_bkg'):
                        templDfPt = templDfPt.query(f"fMlScore0 >= {config_cut['cutvars']['score_bkg']['min'][iPt]} and fMlScore0 < {config_cut['cutvars']['score_bkg']['max'][iPt]}")
                    if config_cut['cutvars'].get('score_FD'):
                        templDfPt = templDfPt.query(f"fMlScore1 >= {config_cut['cutvars']['score_FD']['min'][iPt]} and fMlScore1 < {config_cut['cutvars']['score_FD']['max'][iPt]}")

                hist_templ = ROOT.TH1D(f"h{templName}Raw", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, config_cut['fitrangemin'][iPt], config_cut['fitrangemax'][iPt])
                hist_templ_BR_rew = ROOT.TH1D(f"h{templName}BRRew", ";#it{M}(K#pi#pi) (GeV/#it{c})", nbins, config_cut['fitrangemin'][iPt], config_cut['fitrangemax'][iPt])
                for mass in templDfPt["fM"].to_numpy():
                    hist_templ.Fill(mass)
                hist_templ.Write(f"h{templName}Raw")
                if templName == 'Signal':
                    hist_signal.Add(hist_templ, signalBRNorm)
                    hist_templ_BR_rew.Add(hist_templ, signalBRNorm)
                elif iTemplate == 1:
                    hist_templ_BR_rew.Add(hist_templ, BRnorm)
                    hist_first_templ.Add(hist_templ, BRnorm)                    
                else:
                    hist_templ_BR_rew.Add(hist_templ, BRnorm)
                    hist_templ_total.Add(hist_templ, BRnorm)

                hist_templ_BR_rew.Write(f"h{templName}RescaledBR")
                hist_frac_templ_to_signal[templName].SetBinContent(iPt+1, hist_templ_BR_rew.Integral() / hist_signal.Integral())
                if iTemplate >= 2:
                    hist_frac_templ_to_firsttempl[templName].SetBinContent(iPt+1, hist_templ_BR_rew.Integral() / hist_first_templ.Integral())

            hist_frac_templ_to_signal["TotalBkg"].SetBinContent(iPt+1, hist_templ_total.Integral() / hist_signal.Integral())
            weights_file.mkdir(f"cutset_{cutset}/CombinedSpectra/pt_{ptmin}_{ptmax}/")
            weights_file.cd(f"cutset_{cutset}/CombinedSpectra/pt_{ptmin}_{ptmax}/")
            hist_signal.Smooth(10)
            hist_signal.Write()
            hist_templ_total.Smooth(10)
            hist_templ_total.Write()

        weights_file.cd(f"cutset_{cutset}")
        for iTemplFrac in hist_frac_templ_to_signal.values():
            iTemplFrac.Write()
        for iTemplFrac in hist_frac_templ_to_firsttempl.values():
            iTemplFrac.Write()

    weights_file.Close()

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
