'''
Script for the variation with the BDT cuts.
python cut_variation.py config.yaml an_res.root -c k3050 -r resolution.root -o path/to/output -s text
'''
import ROOT
import yaml
import uproot
import argparse
import numpy as np
import os
import sys
from alive_progress import alive_bar
sys.path.append('../')
sys.path.append('../../../')
from flow_analysis_utils import get_vn_versus_mass, get_centrality_bins, get_cut_sets, apply_pt_weights
from utils.TaskFileLoader import LoadSparseFromTaskHard
from scipy.interpolate import InterpolatedUnivariateSpline
from ROOT import TH1F

def sel_projection(sparseData, sparseReco, varnames, sel, config):
    for var, (low_range, upp_range) in zip(varnames, sel):
        sparseData.GetAxis(config[f'axes'][var]).SetRangeUser(low_range, upp_range)
        for key, sparse in sparseReco.items():
            sparse.GetAxis(config[f'axes_{key}'][var]).SetRangeUser(low_range, upp_range)
    
    return sparseData, sparseReco

def proj_mc_gen(config, sparseGen, ptweights, ptweightsB, Bspeciesweights, sPtWeights, sPtWeightsB):
    histosGen = []
    hGenPtPrompt, hGenPtFD = apply_pt_weights(ptweights, ptweightsB, Bspeciesweights, sparseGen, 
                                              'Pt', sPtWeights, sPtWeightsB)
    histosGen.append(hGenPtPrompt)
    histosGen[-1].SetName('hGenPtPrompt')
    histosGen.append(hGenPtFD)
    histosGen[-1].SetName('hGenPtFD')
    if config.get('enableSecPeak'):
        hGenPtPromptSecPeak = sparseGen['GenSecPeakPrompt'].Projection(0)
        hGenPtFDSecPeak = sparseGen['GenSecPeakFD'].Projection(0)
        histosGen.append(hGenPtPromptSecPeak)
        histosGen[-1].SetName('hGenPtPromptSecPeak')
        histosGen.append(hGenPtFDSecPeak)
        histosGen[-1].SetName('hGenPtFDSecPeak')

    return histosGen

def proj_mc_reco(config, sparseReco, ptweights, ptweightsB, Bspeciesweights, sPtWeights, sPtWeightsB):
    histosReco = []
    for iProjVar in ('InvMass', 'Pt'):
        for key, sparse in sparseReco.items():
            if key != 'RecoPrompt' and key != 'RecoFD':
                histosReco.append(sparse.Projection(config[f'axes_{key}'][iProjVar]))
                histosReco[-1].SetName(f'{key}_{iProjVar}')
            elif key == 'RecoPrompt':
                if iProjVar == 'Pt':
                    hPrompt, hFD = apply_pt_weights(ptweights, ptweightsB, Bspeciesweights, sparseReco, iProjVar, 
                                                    sPtWeights, sPtWeightsB, config['axes_RecoPrompt'][iProjVar], 'Reco')
                else:
                    hPrompt = sparseReco['RecoPrompt'].Projection(config['axes_RecoPrompt'][iProjVar])
                    hFD = sparseReco['RecoFD'].Projection(config['axes_RecoFD'][iProjVar])

                histosReco.append(hPrompt)
                histosReco[-1].SetName(f'RecoPrompt_{iProjVar}')
                histosReco.append(hFD)
                histosReco[-1].SetName(f'RecoFD_{iProjVar}')

    return histosReco

def proj_data(config, sparse, reso, inv_mass_bins):
    histos_data = []
    hist_mass = sparse.Projection(config['axes']['InvMass'])
    hist_mass.SetDirectory(0)
    histos_data.append(hist_mass)
    histos_data[-1].SetName('Data_InvMass')
    hist_vn_sp = get_vn_versus_mass(sparse, inv_mass_bins, config['axes']['InvMass'], config['axes']['sp'])
    hist_vn_sp.SetDirectory(0)
    if reso > 0:
        hist_vn_sp.Scale(1./reso)
    histos_data.append(hist_vn_sp)
    histos_data[-1].SetName('Data_VnVsMass')

    return histos_data

def cut_var(config, an_res_file, centrality, resolution, outputdir, suffix):
    with open(config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)

    # get resolution
    resoFile = ROOT.TFile(resolution, 'READ')
    try:
        det_A = config['detA']
        det_B = config['detB']
        det_C = config['detC']
        reso_hist = resoFile.Get(f'{det_A}_{det_B}_{det_C}/histo_reso')
        reso_hist.SetName('hist_reso')
        reso_hist.SetDirectory(0)
        reso = reso_hist.GetBinContent(1)
    except:
        reso_hist = resoFile.Get(f'hf-task-flow-charm-hadrons/spReso/hSpReso{det_B}{det_C}')
        reso_hist.SetName('hist_reso')
        reso_hist.SetDirectory(0)
        reso = reso_hist.GetBinContent(1)
    reso = reso_hist.GetBinContent(1)

    # load data and mc files and their cent histos
    infilemc = ROOT.TFile(config['MC_filename'], 'read')
    infiledata = ROOT.TFile(an_res_file, 'read')
    nprongs = 2 if config['Dmeson'] == 'Dzero' else 3
    cent_hist_mc = infilemc.Get(f'hf-candidate-creator-{nprongs}prong/hSelCollisionsCent')
    cent_hist_data = infiledata.Get(f'hf-candidate-creator-{nprongs}prong/hSelCollisionsCent')

    os.makedirs(f'{outputdir}/proj', exist_ok=True)
    outfile_resocent = ROOT.TFile(f'{outputdir}/proj/proj_{suffix}.root', 'recreate')
    cent_hist_mc.Write()
    cent_hist_data.Write()
    reso_hist.Write()
    outfile_resocent.Close()

    # compute info for pt weights
    if config.get('ptWeights'):
        ptWeights = uproot.open(config['ptweightPath'])[config['ptweightName']]
        bins = ptWeights.axis(0).edges()
        ptCentW = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
        sPtWeights = InterpolatedUnivariateSpline(ptCentW, ptWeights.values())
    else:
        print('\033[91m WARNING: pt weights will not be provided! \033[0m')
        ptWeights = None
        sPtWeights = None

    if config.get('ptWeightsB'):
        ptWeightsB = uproot.open(config['ptweightBPath'])[config['ptweightBName']]
        bins = ptWeightsB.axis(0).edges()
        ptCentWB = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
        sPtWeightsB = InterpolatedUnivariateSpline(ptCentWB, ptWeightsB.values())
    else:
        print('\033[91m WARNING: B weights will not not be provided! \033[0m')
        ptWeightsB = None
        sPtWeightsB = None

    if config.get('Bspeciesweights'):
        Bspeciesweights = config['Bspeciesweights']
    else:
        print('\033[91m WARNING: B species weights will not be provided! \033[0m')
        Bspeciesweights = None

    # load sparses to be projected
    sparseData = infiledata.Get('hf-task-flow-charm-hadrons/hSparseFlowCharm')
    sparseReco, sparseGen = LoadSparseFromTaskHard(infilemc, config)

    # crosscheck Dzero projections
    if config['Dmeson'] == 'Dzero':
        sparseReco['RecoPrompt'].GetAxis(6).SetRange(2, 2) # make sure it is prompt
        sparseReco['RecoFD'].GetAxis(6).SetRange(3, 3)  # make sure it is non-prompt
        sparseReco['RecoPrompt'].GetAxis(8).SetRange(1, 2)  # make sure it is signal
        sparseReco['RecoFD'].GetAxis(8).SetRange(1, 2)  # make sure it is signal
        sparseGen['GenPrompt'].GetAxis(3).SetRange(2, 2)  # make sure it is prompt
        sparseGen['GenFD'].GetAxis(3).SetRange(3, 3)  # make sure it is non-prompt
    if config.get('enableRef'):
        sparseReco['RecoRefl'].GetAxis(8).SetRange(3, 4)  # make sure it is reflection
        sparseReco['RecoReflPrompt'].GetAxis(8).SetRange(3, 4)  # make sure it is reflection
        sparseReco['RecoReflPrompt'].GetAxis(6).SetRange(2, 2)  # make sure it is prompt reflection
        sparseReco['RecoReflFD'].GetAxis(8).SetRange(3, 4)  # make sure it is reflection
        sparseReco['RecoReflFD'].GetAxis(6).SetRange(3, 3)  # make sure it is non-prompt reflection
    #TODO: safety checks for Dmeson reflecton and secondary peak

    # centrality projection
    _, cent_bins = get_centrality_bins(centrality)
    cent_min = cent_bins[0]
    cent_max = cent_bins[1]
    sparseData.GetAxis(config['axes']['cent']).SetRangeUser(cent_min, cent_max)
    for key, sparse in sparseReco.items():
        sparse.GetAxis(config[f'axes_{key}']['cent']).SetRangeUser(cent_min, cent_max)
    for key, sparse in sparseGen.items():
        sparse.GetAxis(config[f'axes_{key}']['cent']).SetRangeUser(cent_min, cent_max)
    
    # obtain cut set for cut variation
    pt_mins = config['ptmins']
    pt_maxs = config['ptmaxs']
    nCutSets, sels, cutVars = get_cut_sets(pt_mins, config['cut_variation']['bdt_cut'], 
                                           config['minimisation']['correlated'])

    # run over pt bins
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        print(f"Processing pt bin: {pt_min} -- {pt_max}")
        ptLowLabel = pt_min * 10
        ptHighLabel = pt_max * 10
        
        # pt projection
        sparseData.GetAxis(config['axes']['Pt']).SetRangeUser(pt_min, pt_max)
        for key, sparse in sparseReco.items():
            sparse.GetAxis(config[f'axes_{key}']['Pt']).SetRangeUser(pt_min, pt_max)
        for key, sparse in sparseGen.items():
            sparse.GetAxis(config[f'axes_{key}']['Pt']).SetRangeUser(pt_min, pt_max)
        
        # gen histos just projected over pt
        histos_mc_gen = proj_mc_gen(config, sparseGen, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB)

        # cut variation for data and reco mc
        with alive_bar(nCutSets, title='Processing BDT cuts') as bar:
            for sel in sels[ipt]:
                cutset_dir = f'pt_{ptLowLabel}_{ptHighLabel}_cent_{cent_min}_{cent_max}'
                sparseData, sparseReco = sel_projection(sparseData, sparseReco, cutVars, sel, config)
                for iVarName, iVarRange in zip(cutVars, sel):
                    cutset_dir = cutset_dir + f'_{iVarName}_{iVarRange[0]:.2f}_{iVarRange[1]:.2f}'
                print(f"Processing BDT cut: {cutset_dir}")
                # histos_data = proj_data(config, sparseData, reso, config['inv_mass_bins'][ipt])
                # histos_mc_reco = proj_mc_reco(config, sparseReco, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB)

                outfile = ROOT.TFile(f'{outputdir}/proj/proj_{cutset_dir}.root', 'RECREATE')
                # for histo in histos_data:
                #     histo.Write()
                # for histo in histos_mc_reco:
                #     histo.Write()

                # outfile.mkdir('gen/')
                # outfile.cd('gen/')
                # for histo in histos_mc_gen:
                #     histo.Write()
            
                outfile.Close()
                bar()

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("an_res_file", metavar="text",
                        default="an_res.root", help="input ROOT file with anres")
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k3050", help="centrality class")
    parser.add_argument("--resolution", "-r", metavar="text",
                        default="reso.root", help="resolution file")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    cut_var(
        config=args.config,
        an_res_file=args.an_res_file,
        centrality=args.centrality,
        resolution=args.resolution,
        outputdir=args.outputdir,
        suffix=args.suffix
    )
