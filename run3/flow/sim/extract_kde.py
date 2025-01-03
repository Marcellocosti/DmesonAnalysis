import uproot
import pandas as pd
import numpy as np
import ROOT
from ROOT import TFile, TH1D

import sys
sys.path.append('../')
from kde_producer import kde_producer, kde_producer_sim

ptMins = [0, 1, 2, 4, 6, 8, 12, 16]
ptMaxs = [1, 2, 4, 6, 8, 12, 16, 50]

DsFilePath = "DsToKKPiMisidKPi.root"
DsTreePath = 'DsKKPiMisidKPi'
TDsData = []
with uproot.open(DsFilePath) as f:
    for iKey, key in enumerate(f.keys()):
        if DsTreePath in key:
            TDsData = f[key].arrays(library='pd')
            

outFileDs = TFile('Templ_DsKKPiMisidKPi.root', 'recreate')

for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    outFileDs.mkdir(f'Ds_pt_{ptMin}_{ptMax}')
    outFileDs.cd(f'Ds_pt_{ptMin}_{ptMax}')
    kde_phi = kde_producer_sim("DsToKKPiMisidKPi.root", 'DsKKPiMisidKPi', ptMin, ptMax, 1)
    kde_K892 = kde_producer_sim("DsToKKPiMisidKPi.root", 'DsKKPiMisidKPi', ptMin, ptMax, 2)
    fkde_phi = kde_phi.GetFunction(1000)
    fkde_phi.SetName('fDsToPhiPi')
    fkde_phi.Write()
    fkde_K892 = kde_K892.GetFunction(1000)
    fkde_K892.SetName('fDsToK892')
    fkde_K892.Write()
    
    hDsToPhiPi = TH1D('hDsToPhiPi', 'hDsToPhiPi', 1500, 1.3, 2.5)
    hDsToK892 = TH1D('hDsToK892', 'hDsToK892', 1500, 1.3, 2.5)
    pTTDsData = TDsData.query(f"{ptMin} <= pt < {ptMax}")
    TDsDataToPhi = pTTDsData.query(f"decayflag == 1")
    TDsDataToK892 = pTTDsData.query(f"decayflag == 2")
    hDsToPhiPi.FillN(len(TDsDataToPhi['mass']), np.asarray(TDsDataToPhi['mass'], dtype=np.float64), np.asarray(TDsDataToPhi['decayflag'], dtype=np.float64))
    hDsToK892.FillN(len(TDsDataToK892['mass']), np.asarray(TDsDataToK892['mass'], dtype=np.float64), np.asarray(TDsDataToK892['decayflag'], dtype=np.float64))
    hDsToPhiPi.Write()
    hDsToK892.Write()
    
outFileDs.Close()

DplusFilePath = "DplusToKKPiMisidKPi.root"
DplusTreePath = 'DplusKKPiMisidKPi'
TDplusData = []
with uproot.open(DplusFilePath) as f:
    for iKey, key in enumerate(f.keys()):
        if DplusTreePath in key:
            TDplusData = f[key].arrays(library='pd')
            

outFileDplus = TFile('Templ_DplusKKPiMisidKPi.root', 'recreate')

for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    outFileDplus.mkdir(f'Dplus_pt_{ptMin}_{ptMax}')
    outFileDplus.cd(f'Dplus_pt_{ptMin}_{ptMax}')
    kde_phi = kde_producer_sim("DplusToKKPiMisidKPi.root", 'DplusKKPiMisidKPi', ptMin, ptMax, 1)
    kde_K892 = kde_producer_sim("DplusToKKPiMisidKPi.root", 'DplusKKPiMisidKPi', ptMin, ptMax, 2)
    fkde_phi = kde_phi.GetFunction(1000)
    fkde_phi.SetName('fDplusToPhiPi')
    fkde_phi.Write()
    fkde_K892 = kde_K892.GetFunction(1000)
    fkde_K892.SetName('fDplusToK892')
    fkde_K892.Write()
    
    hDplusToPhiPi = TH1D('hDplusToPhiPi', 'hDplusToPhiPi', 1500, 1.3, 2.5)
    hDplusToK892 = TH1D('hDplusToK892', 'hDplusToK892', 1500, 1.3, 2.5)
    pTTDplusData = TDplusData.query(f"{ptMin} <= pt < {ptMax}")
    TDplusDataToPhi = pTTDplusData.query(f"decayflag == 1")
    TDplusDataToK892 = pTTDplusData.query(f"decayflag == 2")
    hDplusToPhiPi.FillN(len(TDplusDataToPhi['mass']), np.asarray(TDplusDataToPhi['mass'], dtype=np.float64), np.asarray(TDplusDataToPhi['decayflag'], dtype=np.float64))
    hDplusToK892.FillN(len(TDplusDataToK892['mass']), np.asarray(TDplusDataToK892['mass'], dtype=np.float64), np.asarray(TDplusDataToK892['decayflag'], dtype=np.float64))
    hDplusToPhiPi.Write()
    hDplusToK892.Write()
    
outFileDplus.Close()
