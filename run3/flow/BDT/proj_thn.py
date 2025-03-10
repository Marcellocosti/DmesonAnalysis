'''
Script to project the MC distributions and apply the pt weights from the AnRes.root of Dtask
python3 proj_thn_mc.py config_flow.yml config_cutset.yml -o path/to/output -s text
                                                        --ptWeights path/to/file histName 
                                                        --ptWeightsB path/to/file histName
'''
import ROOT
import uproot
import yaml
import argparse
import sys
import os
from ROOT import TFile, TObject
from alive_progress import alive_bar
from scipy.interpolate import InterpolatedUnivariateSpline
from sparse_dicts import get_sparses 
sys.path.append('..')
from flow_analysis_utils import get_vn_versus_mass, get_centrality_bins, reweight_histo

### please fill your path of DmesonAnalysis
sys.path.append('../../..')

def proj_data(sparse_flow, ptMin, ptMax, centMin, centMax, axes, inv_mass_bins, reso, writeopt, syst=False):

    print(f"sparse_flow: {sparse_flow}")
    if isinstance(sparse_flow, dict):
        for isparse, (_, sparse) in enumerate(sparse_flow.items()):
            hist_mass_temp = sparse.Projection(axes['Flow']['Mass'])
            # REVIEW: in case the Potential memory leak
            hist_mass_temp.SetName(f'hist_mass_{isparse}')
            hist_mass_temp.SetDirectory(0)
            # REVIEW: I would suggest to keep th fd score distribution of a dedicated pt bin,
            # from my experience, it could help us to choose a proper cutset
            if not syst:
                hist_fd_temp = sparse.Projection(axes['Flow']['score_FD'])
                hist_fd_temp.SetName(f'hist_fd_cent{cent_min}_{cent_max}_pt{ptMin}_{ptMax}_{isparse}')
                hist_bkg_temp = sparse.Projection(axes['Flow']['score_bkg'])
                hist_bkg_temp.SetName(f'hist_bkg_cent{cent_min}_{cent_max}_pt{ptMin}_{ptMax}_{isparse}')

            if isparse == 0:
                hist_mass = hist_mass_temp.Clone('hist_mass')
                hist_mass.SetDirectory(0)
                hist_mass.Reset()
                if not syst:
                    hist_fd = hist_fd_temp.Clone('hist_fd')
                    hist_fd.SetDirectory(0)
                    hist_fd.Reset()
                    hist_bkg = hist_fd_temp.Clone('hist_bkg')
                    hist_bkg.SetDirectory(0)
                    hist_bkg.Reset()

            hist_mass.Add(hist_mass_temp)
            if not syst:
                hist_fd.Add(hist_fd_temp)
                hist_bkg.Add(hist_bkg_temp)
            
        print(f"list(sparse_flow.values()): {list(sparse_flow.values())}")
        print(f"axes['Flow']['Mass']: {axes['Flow']['Mass']}")
        print(f"axes['Flow']['sp']: {axes['Flow']['sp']}")
        hist_vn_sp = get_vn_versus_mass(list(sparse_flow.values()), inv_mass_bins, axes['Flow']['Mass'], axes['Flow']['sp'])
        hist_vn_sp.SetDirectory(0)
        if reso > 0:
            hist_vn_sp.Scale(1./reso)
    else:
        hist_mass = sparse_flow.Projection(axes['Flow']['Mass'])
        hist_mass.SetDirectory(0)
        if not syst:
            hist_fd = sparse_flow.Projection(axes['Flow']['score_FD'])
            hist_fd.SetDirectory(0)
            hist_bkg = sparse_flow.Projection(axes['Flow']['score_bkg'])
            hist_bkg.SetDirectory(0)
        hist_vn_sp = get_vn_versus_mass(sparse_flow, inv_mass_bins, axes['Flow']['Mass'], axes['Flow']['sp'])
        hist_vn_sp.SetDirectory(0)
        if reso > 0:
            hist_vn_sp.Scale(1./reso)

    hist_mass.Write(f'hist_mass_cent{centMin}_{centMax}_pt{ptMin}_{ptMax}', writeopt)
    hist_vn_sp.Write(f'hist_vn_sp_pt{ptMin}_{ptMax}', writeopt)
    if not syst:
        hist_fd.Write(f'hist_fd_cent{cent_min}_{cent_max}_pt{ptMin}_{ptMax}', writeopt)
        hist_bkg.Write(f'hist_bkg_cent{cent_min}_{cent_max}_pt{ptMin}_{ptMax}', writeopt)

def proj_mc_reco(sparsesReco, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB, writeopt):
    
    for key, sparse in sparsesReco.items():
        if key != 'RecoPrompt' and key != 'RecoFD':
            for iProjVar in ('Mass', 'Pt'):
                sparse.Projection(axes[key][iProjVar]).Write(f'h{key}{iProjVar}')

    hMassPrompt = sparsesReco['RecoPrompt'].Projection(axes['RecoPrompt']['Mass'])
    hMassPrompt.SetName(f'hPromptMass_{ptMin}_{ptMax}')
    hMassFD = sparsesReco['RecoFD'].Projection(axes['RecoFD']['Mass'])    
    hMassFD.SetName(f'hFDMass_{ptMin}_{ptMax}')

    ### project pt prompt
    hPtPrompt = sparsesReco['RecoPrompt'].Projection(axes['RecoPrompt']['Pt'])
    if ptWeights:
        hPtPrompt = reweight_histo(hPtPrompt, sPtWeights, 'hPromptPt') 
    ### project pt FD
    if ptWeightsB:
        if Bspeciesweights:
            hPtFD = reweight_histo(sparsesReco['RecoFD'].Projection(axes['RecoFD']['Pt'], axes['RecoFD']['pt_bmoth'], axes['RecoFD']['flag_bhad']), 
                                   sPtWeightsB, 'hFDPt', Bspeciesweights)
        else:
            hPtFD = reweight_histo(sparsesReco['RecoFD'].Projection(axes['RecoFD']['pt_bmoth'], axes['RecoFD']['Pt']), sPtWeightsB, 'hFDPt') # 2D projection: Projection(ydim, xdim)
    elif ptWeights:
        hPtFD = reweight_histo(sparsesReco['RecoFD'].Projection(axes['RecoFD']['Pt']), sPtWeights, 'hFDPt')
    elif Bspeciesweights:
        hPtFD = reweight_histo(sparsesReco['RecoFD'].Projection(axes['RecoFD']['flag_bhad'], axes['RecoFD']['Pt']), [], 'hFDPt', Bspeciesweights) # 2D projection: Projection(ydim, xdim)
    else:
        hPtFD = sparsesReco['RecoFD'].Projection(axes['RecoFD']['Pt'])

    ## write the output      
    hMassPrompt.Write('hPromptMass', writeopt)
    hMassFD.Write('hFDMass', writeopt)
    hPtPrompt.Write('hPromptPt', writeopt)
    hPtFD.Write('hFDPt', writeopt)

def proj_mc_gen(sparsesGen, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB, writeopt):

    for key, sparse in sparsesGen.items():
        if key != 'GenPrompt' and key != 'GenFD':
            sparse.Projection(axes[key]['Pt']).Write(f'h{key}Pt')

    ### prompt
    hGenPtPrompt = sparsesGen['GenPrompt'].Projection(axes['GenPrompt']['Pt'])
    if ptWeights:
        hGenPtPrompt = reweight_histo(hGenPtPrompt, sPtWeights, 'hPromptGenPt')
    ### FD
    if ptWeightsB:
        if Bspeciesweights:
            hGenPtFD = reweight_histo(sparsesGen['GenFD'].Projection(axes['GenFD']['Pt'], axes['GenFD']['pt_bmoth'], axes['GenFD']['flag_bhad']), 
                                      sPtWeightsB, 'hFDGenPt', Bspeciesweights)
        else:
            hGenPtFD = reweight_histo(sparsesGen['GenFD'].Projection(axes['GenFD']['pt_bmoth'], axes['GenFD']['Pt']), sPtWeightsB, 'hFDGenPt') # 2D projection: Projection(ydim, xdim)
    elif ptWeights:
        hGenPtFD = reweight_histo(sparsesGen['GenFD'].Projection(axes['GenFD']['Pt']), sPtWeights, 'hFDGenPt')
    elif Bspeciesweights:
        hGenPtFD = reweight_histo(sparsesGen['GenFD'].Projection(axes['GenFD']['flag_bhad'], axes['GenFD']['Pt']), [], 'hFDPt', Bspeciesweights) # 2D projection: Projection(ydim, xdim)
    else:
        hGenPtFD = sparsesGen['GenFD'].Projection(axes['GenFD']['Pt'])

    ## write the output
    hGenPtPrompt.Write('hPromptGenPt', writeopt)
    hGenPtFD.Write('hFDGenPt', writeopt)

def pt_weights_info(ptweights, ptweightsB):
    """Get pt weights and return weights flags with spline

    Args:
        ptweights (list): [file path, histogram name] for pt weights
        ptweightsB (list): [file path, histogram name] for B pt weights

    Outputs:
        ptWeights (bool): ptWeights flag
        ptWeightsB (bool): ptWeightsB flag
        Bspeciesweights (str): B species weights #TODO
        sPtWeights (spline): Spline for ptWeights interpolation
        sPtWeightsB (spline): Spline for ptWeightsB weights interpolation
    """

# REVIEW: the ptWeights inputed is a list, but the ptWeights outputed is a TH1D object
# and actually ptweights is used as a flag
    # compute info for pt weights
    if ptweights != []:
        with uproot.open(ptweights[0]) as f:
            hPtWeights = f[ptweights[1]]
            bins = hPtWeights.axis(0).edges()
            ptCentW = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
            sPtWeights = InterpolatedUnivariateSpline(ptCentW, hPtWeights.values())
        ptWeights = True
    else:
        print('\033[91m WARNING: pt weights will not be provided! \033[0m')
        ptWeights = False
        sPtWeights = None

    if ptweightsB != []:
        with uproot.open(ptweightsB[0]) as f:
            hPtWeightsB = f[ptweightsB[1]]
            bins = hPtWeightsB.axis(0).edges()
            ptCentWB = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
            sPtWeightsB = InterpolatedUnivariateSpline(ptCentWB, hPtWeightsB.values())
        ptWeightsB = True
    else:
        print('\033[91m WARNING: B weights will not not be provided! \033[0m')
        ptWeightsB = False
        sPtWeightsB = None

    if config.get('Bspeciesweights'):
        Bspeciesweights = config['Bspeciesweights']
    else:
        print('\033[91m WARNING: B species weights will not be provided! \033[0m')
        Bspeciesweights = None
    
    return ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="flow configuration file")
    parser.add_argument('anres_dir', metavar='text', 
                        nargs='*', help='input ROOT files with anres')
    parser.add_argument('--cutsetConfig', "-cc", metavar='text', type=str, nargs='?',
                        const=None, default='cutsetConfig.yaml',
                        help='Optional cutset configuration file (default: cutsetConfig.yaml)')    
    parser.add_argument("--proj_data", action="store_true", 
                        help="Flag to project data")
    parser.add_argument("--proj_mc", action="store_true", 
                        help="Flag to project MC")
    parser.add_argument('--preprocessed', action='store_true', 
                        help='Determines whether the sparses are pre-processed')
    parser.add_argument("--systematics", action="store_true",
                        help="cutset based AnRes files")
    parser.add_argument("--ptweights", "-w", metavar="text", nargs=2, required=False,
                        default=[], help="path to pt weights file and histogram name")
    parser.add_argument("--ptweightsB", "-wb", metavar="text", nargs=2, required=False,
                        default=[], help="path to pt weightsB file and histogram name")
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k3050", help="centrality class")
    parser.add_argument("--resolution", "-r", metavar="text",
                        default="reso.root", help="resolution file")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()
    
    print(f"args.pre_processed: {args.preprocessed}")
    print(f"args.proj_data: {args.proj_data}")
    print(f"args.proj_mc: {args.proj_mc}")
    with open(args.config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)


    cent, (cent_min, cent_max) = get_centrality_bins(args.centrality)
    print(f"cent_min: {cent_min}")
    print(f"args.suffix: {args.suffix}")
    os.makedirs(f'{args.outputdir}/proj', exist_ok=True)
    outfilename = f'{args.outputdir}/proj/proj_{args.suffix}'
    create_new_file = True
    write_opt_data = 0
    write_opt_mc = 0
    write_opt_cent_reso = 0
    if args.proj_data and args.proj_mc:
        print(f"No existing previous projections, creating new file and project data and mc!")
        outfile = TFile(outfilename + '.root', 'RECREATE')
    else:
        suffixCode = (args.suffix).split("_")[-1]
        projFiles = [f'{outfilename}*' for file in os.listdir(f'{args.outputdir}/proj/') if file.endswith(f'_{suffixCode}.root')]
        if len(projFiles) == 0:
            print(f"No existing previous projections, creating new file and project data ({args.proj_data}) or mc ({args.proj_mc})!")
            outfile = TFile(outfilename + '.root', 'RECREATE')
        else:
            create_new_file = False
            print(f"Found previous projections, updating existing file!")
            outfile = TFile.Open(outfilename + '.root', 'UPDATE')
            if args.proj_data:
                write_opt_data = TObject.kOverwrite 
                write_opt_cent_reso = TObject.kOverwrite
            if args.proj_mc:
                write_opt_mc = TObject.kOverwrite
                write_opt_cent_reso = TObject.kOverwrite

    print(f"create_new_file: {create_new_file}")
    print(f"outfilename + '.root': {outfilename+ '.root'}")
    outfile_dir = 'hf-candidate-creator-2prong' if config['Dmeson'] == 'Dzero' else 'hf-candidate-creator-3prong'
    infilemc = TFile.Open(config['eff_filename'], 'r')
    histo_cent = infilemc.Get(f'{outfile_dir}/hSelCollisionsCent')
    histo_cent.GetXaxis().SetRangeUser(cent_min, cent_max)
    resofile = TFile.Open(args.resolution, 'r')
    try:
        det_A = config['detA']
        det_B = config['detB']
        det_C = config['detC']
        histo_reso = resofile.Get(f'{det_A}_{det_B}_{det_C}/histo_reso_delta_cent')
        histo_reso.SetDirectory(0)
        reso = histo_reso.GetBinContent(1)
    except:
        histo_reso = resofile.Get(f'hf-task-flow-charm-hadrons/spReso/hSpReso{det_B}{det_C}')
        histo_reso.SetDirectory(0)
        reso = histo_reso.GetBinContent(1)
    
    outfile.cd()
    histo_reso.Write('hist_reso', write_opt_cent_reso)
    if create_new_file:
        outfile.mkdir(outfile_dir)
    outfile.cd(outfile_dir)
    histo_cent.Write('hSelCollisionsCent', write_opt_cent_reso)
    resofile.Close()
    infilemc.Close()

    # with open(args.cutsetConfig, 'r') as ymlCutSetFile:
    #     cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
    # cutVars = cutSetCfg['cutvars']

    # iCut = '00'
    # print(f"args.cutsetConfig: {args.cutsetConfig}")
    # if args.cutsetConfig != 'cutsetConfig.yaml':
    #     print(f"Entered loop")
    with open(args.cutsetConfig, 'r') as ymlCutSetFile:
        cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
        iCut = f"{int(cutSetCfg['icutset']):02d}"
    cutVars = cutSetCfg['cutvars']

    # load thnsparse
    # # REVIEW chuntai: 
    # # for the main workflow, only the config_flow
    print(f"args.proj_data: {args.proj_data}")
    print(f"args.proj_mc: {args.proj_mc}")
    sparsesFlow, sparsesReco, sparsesGen, axes = get_sparses(config, args.proj_data, args.proj_mc, args.proj_mc, config.get('anresdir', []), 
                                                             args.preprocessed, f'{config.get("skim_out_dir", "")}', args.systematics, iCut)
    print(f"\n")
    print(f"sparsesFlow: {sparsesFlow}")
    print(f"\n")
    if not args.preprocessed:
        for key, iSparse in sparsesFlow.items():
            iSparse.GetAxis(axes['Flow']['cent']).SetRangeUser(cent_min, cent_max)
    for key, iSparse in sparsesGen.items():
        iSparse.GetAxis(axes[key]['cent']).SetRangeUser(cent_min, cent_max)
    for key, iSparse in sparsesReco.items():
        iSparse.GetAxis(axes[key]['cent']).SetRangeUser(cent_min, cent_max)

    # compute info for pt weights
    if args.proj_mc:
        ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB = pt_weights_info(args.ptweights, args.ptweightsB)

    with alive_bar(len(cutVars['Pt']['min']), title='Processing pT bins') as bar:
        for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
            print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
            ptLowLabel = ptMin * 10
            ptHighLabel = ptMax * 10
            ptcentdir = f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}'  
            if create_new_file:    
                print(f"creating new directory: {ptcentdir}")      
                outfile.mkdir(ptcentdir)
            outfile.cd(ptcentdir)

            if args.proj_data:
                if args.preprocessed:
                    print('PREPROCESSED')
                    print('Taking pre-processed AnRes for data!')
                    if not args.systematics:
                        sparsesFlow[f"Flow_{ptLowLabel}_{ptHighLabel}"].GetAxis(axes['Flow']['score_FD']).SetRangeUser(cutVars['score_FD']['min'][iPt], cutVars['score_FD']['max'][iPt])
                        if 'score_bkg' in config['axestokeep']:
                            print(f"Cutting on bkg on pre-processed AnRes!")
                            sparsesFlow[f"Flow_{ptLowLabel}_{ptHighLabel}"].GetAxis(axes['Flow']['score_bkg']).SetRangeUser(cutVars['score_bkg']['min'][iPt], cutVars['score_bkg']['max'][iPt])
                    proj_data(sparsesFlow[f"Flow_{ptLowLabel}_{ptHighLabel}"], ptMin, ptMax, cent_min, cent_max, axes, config['inv_mass_bins'][iPt], reso, write_opt_data, args.systematics)
                    outfile.cd(ptcentdir)
                    print(f"Projected data!")

                else:
                    for iSparse, (key, sparse) in enumerate(sparsesFlow.items()):
                        for iVar in cutVars:
                            sparse.GetAxis(axes['Flow'][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                    proj_data(sparsesFlow, ptMin, ptMax, cent_min, cent_max, axes, config['inv_mass_bins'][iPt], reso, write_opt_data)
                    print("Projected data!")
            else:
                print("Kept data from previous projections!")

            if args.proj_mc:
                for iVar in cutVars:
                    for key, iSparse in sparsesReco.items():
                        iSparse.GetAxis(axes[key][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                    if iVar == 'Pt':
                        for key, iSparse in sparsesGen.items():
                            iSparse.GetAxis(axes[key][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                    if iVar == 'score_FD' or iVar == 'score_bkg':
                        print(f'{iVar}: {cutVars[iVar]["min"][iPt]} < {iVar} < {cutVars[iVar]["max"][iPt]}')

                proj_mc_reco(sparsesReco, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB, write_opt_mc)
                print("Projected mc reco!")
                proj_mc_gen(sparsesGen, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB, write_opt_mc)
                print("Projected mc gen!")
            else:
                print("Kept mc from previous projections!")

            # if args.systematics and not args.proj_mc:
            # if args.systematics:
            #     mc_histos, mc_histos_names = [], []
            #     cutsetConfig = args.cutsetConfig
            #     icutset = f"{cutSetCfg['icutset']:02d}"

            #     proj_file = cutsetConfig.replace('config', 'proj').replace(f'cutset_uncorr_{icutset}.yml', f'proj_uncorr_{icutset}.root')
                
            #     proj = TFile.Open(proj_file, 'read')
            #     print(f"Opening {proj_file}")
            #     proj.cd(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}')
            #     histos = proj.GetDirectory(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}').GetListOfKeys()
            #     print(f"histos: {histos}")
            #     for histo in histos:
            #         print(f"histo: {histo.GetName()}")
            #         if 'hist' not in histo.GetName():
            #             mc_histos_names.append(histo.GetName())
            #             histo = histo.ReadObj() 
            #             mc_histos.append(histo)
            #             mc_histos[-1].SetDirectory(0)
            #     proj.Close()
            #     outfile.cd(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}')
            #     for histo, name in zip(mc_histos, mc_histos_names):
            #         histo.Write(name)
            #     print(f"Projected systematics!")

            print('\n')
            bar()
    
    outfile.Close()
