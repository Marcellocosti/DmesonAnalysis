import sys
import os
import ctypes
import itertools
from ROOT import TFile, TCanvas, TLegend, TH1D, TH1F # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, gPad, gInterpreter, kBlack, kRed, kBlue, kFullCircle # pylint: disable=import-error,no-name-in-module
import numpy as np
gInterpreter.ProcessLine(f'#include "/home/mdicosta/alice/DmesonAnalysis/run3/flow/invmassfitter/InvMassFitter.cxx"')
gInterpreter.ProcessLine(f'#include "/home/mdicosta/alice/DmesonAnalysis/run3/flow/invmassfitter/VnVsMassFitter.cxx"')
from ROOT import InvMassFitter, VnVsMassFitter
sys.path.append("/home/mdicosta/alice/DmesonAnalysis/")
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, DivideCanvas

gROOT.SetBatch(1)
SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, padtopmargin=0.12, opttitle=1)
cent = '3040'
centMinMax = [30, 40]

# pt bins
ptMins = [
           1.0,
           1.5,
           2.0,
           2.5,
           3.0,
           3.5,
           4.0,
           5.0,
           6.0,
           7.0,
           8.0,
           10.0,
           12.0,
           16.0
       ]
ptMaxs = [
           1.5,
           2.0,
           2.5,
           3.0,
           3.5,
           4.0,
           5.0,
           6.0,
           7.0,
           8.0,
           10.0,
           12.0,
           16.0,
           24.0
       ]
massMins = [
          1.75, 
          1.75,
          1.75, 
          1.75,
          1.75, 
          1.75, 
          1.70, 
          1.75, 
          1.70, 
          1.70,
          1.65, 
          1.65, 
          1.65, 
          1.65, 
       ]
massMaxs = [
          2.00,
          2.00,
          2.00, 
          1.95, 
          2.00, 
          2.00, 
          2.00, 
          1.95,
          1.95,
          1.95,
          1.97, 
          2.04, 
          2.17, 
          2.17, 
       ]
rebins = [
          2,
          2,
          1, 
          1, 
          1, 
          1, 
          1, 
          1,
          1,
          1,
          1, 
          1, 
          1, 
          1, 
       ]

# sanity check of fit configuration
SgnFunc = []
SgnFuncStr = [
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            "kDoubleCBAsymm",
            ]

for iPt, sgnStr in enumerate(SgnFuncStr):
    if sgnStr == 'kGaus':
        SgnFunc.append(InvMassFitter.kGaus)
    elif sgnStr == 'kDoubleCBAsymm':
        SgnFunc.append(InvMassFitter.kDoubleCBAsymm)
    elif sgnStr == 'kDoubleCBSymm':
        SgnFunc.append(InvMassFitter.kDoubleCBSymm)
    else:
        print('ERROR: only kGaus, k2Gaus and k2GausSigmaRatioPar signal functions supported! Exit!')
        sys.exit()

print(f"SgnFunc: {SgnFunc}")

# load histos
infile = TFile.Open("/home/mdicosta/FlowDplus/FinalResults/crystalball/cutvar_projmcinclusive/proj/proj_projmcinclusive_00.root")
if not infile or not infile.IsOpen():
    print(f'ERROR: file "/home/mdicosta/FlowDplus/FinalResults/crystalball/cutvar_projmcinclusive/proj/proj_projmcinclusive_00.root" cannot be opened! Exit!')
    sys.exit()
hMasses = []
fTotFuncMass, fSgnFuncMass = [], []

for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    print(f'loading: cent_bins{centMinMax[0]}_{centMinMax[1]}/pt_bins{ptMin}_{ptMax}/hPromptMass')
    hMasses.append(infile.Get(f'cent_bins{centMinMax[0]}_{centMinMax[1]}/pt_bins{ptMin}_{ptMax}/hPromptMass'))
    hMasses[iPt].SetDirectory(0)
    SetObjectStyle(hMasses[iPt], color=kBlack, markerstyle=kFullCircle)
infile.Close()

print(f"hMasses: {hMasses}")

ptLims = list(ptMins)
nPtBins = len(ptMins)
ptLims.append(ptMaxs[-1])
ptBinsArr = np.asarray(ptLims, 'd')
ptTit = '#it{p}_{T} (GeV/#it{c})'
hRawYields = TH1D('hRawYields', f';{ptTit};raw yield', nPtBins, ptBinsArr)
hSigma = TH1D('hSigma', f';{ptTit};#sigma', nPtBins, ptBinsArr)
hMean = TH1D('hMean', f';{ptTit};mean', nPtBins, ptBinsArr)
hRedChi2 = TH1D('hRedChi2', f';{ptTit};#chi^{{2}}/#it{{ndf}}', nPtBins, ptBinsArr)

SetObjectStyle(hRawYields, color=kRed, markerstyle=kFullCircle)
SetObjectStyle(hSigma, color=kRed, markerstyle=kFullCircle)
SetObjectStyle(hMean, color=kRed, markerstyle=kFullCircle)
SetObjectStyle(hRedChi2, color=kRed, markerstyle=kFullCircle)

### Asymmetric case
CBInitParsAsymm = [
    # SgnInt, mu, width, alpha1, n1, alpha2, n2 --> Asymm
    [["sgnInt", 50, 10, 200], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 200, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 300, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 450, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 700, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 700, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 700, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 700, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 400, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 350, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 750, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 800, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
    [["sgnInt", 1300, 10, 2400], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha1", 2, 1, 3], ["N1", 4, 1.1, 10], ["Alpha2", 5, 3, 10], ["N2", 4, 1.1, 10]],
]

### Symmetric case
CBInitParsSymm = [
    # SgnInt, mu, width, alpha, n --> Symm
    [["sgnInt", 100, 10, 200], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 100, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
    [["sgnInt", 1000, 10, 2000], ["Mean", 1.87, 1.86, 1.88], ["Sigma", 0.0125, 0.005, 0.025], ["Alpha", 2, 1, 3], ["N", 4, 1.1, 10]],
]

#_____________________________________________________
# Mass fit
outFile = TFile(f"/home/mdicosta/FlowDplus/FinalResults/crystalball/FitMC_{SgnFuncStr[0]}.root", "recreate")
massFitters = []
for iPt, (hMass, ptMin, ptMax, reb, sgnEnum, massMin, massMax) in enumerate(zip(hMasses, ptMins, ptMaxs, rebins, SgnFunc, massMins, massMaxs)):
    iCanv = iPt
    hMass.Rebin(reb)
    binWidth = hMass.GetBinWidth(1)
    hMass.SetTitle((f'{ptMin:0.1f} < #it{{p}}_{{T}} < {ptMax:0.1f} GeV/#it{{c}};M(KKPi);'
                                    f'Counts per {binWidth*1000:.0f} MeV/#it{{c}}^{{2}}'))
    hMass.SetName(f'Mass{iPt}')
    SetObjectStyle(hMass, color=kRed-3, markerstyle=kFullCircle, markersize=0.8)

    print(f'Fitting {ptMin} - {ptMax} GeV/c')
    hMassForFit = TH1F()
    hMass.Copy(hMassForFit)
    massFitters.append(InvMassFitter(hMassForFit,  massMin, massMax, InvMassFitter.kNoBk, sgnEnum))
    massFitters[iPt].SetUseLikelihoodFit()
    massFitters[iPt].SetBoundGaussianMean(1.87, 1.86, 1.88)
    if sgnEnum == 3:
        massFitters[iPt].SetInitPars(CBInitParsAsymm[iPt])
    if sgnEnum == 4:
        massFitters[iPt].SetInitPars(CBInitParsSymm[iPt])
    massFitters[iPt].MassFitter(False)

    # collect fit results
    rawyield = massFitters[iPt].GetRawYield()
    rawyielderr = massFitters[iPt].GetRawYieldError()
    sigma = massFitters[iPt].GetSigma()
    sigmaerr = massFitters[iPt].GetSigmaUncertainty()
    mean = massFitters[iPt].GetMean()
    meanerr = massFitters[iPt].GetMeanUncertainty()
    redchi2 = massFitters[iPt].GetReducedChiSquare()
    signif, signiferr = ctypes.c_double(), ctypes.c_double()
    sgn, sgnerr = ctypes.c_double(), ctypes.c_double()
    bkg, bkgerr = ctypes.c_double(), ctypes.c_double()
    massFitters[iPt].Significance(3, signif, signiferr)
    massFitters[iPt].Signal(3, sgn, sgnerr)
    massFitters[iPt].Background(3, bkg, bkgerr)

    hRawYields.SetBinContent(iPt+1, rawyield)
    hRawYields.SetBinError(iPt+1, rawyielderr)
    hSigma.SetBinContent(iPt+1, sigma)
    hSigma.SetBinError(iPt+1, sigmaerr)
    hMean.SetBinContent(iPt+1, mean)
    hMean.SetBinError(iPt+1, meanerr)
    hRedChi2.SetBinContent(iPt+1, redchi2)
    hRedChi2.SetBinError(iPt+1, 1.e-20)

    # plot the results
    outFile.mkdir(f"cent_bins{centMinMax[0]}_{centMinMax[1]}/pt_bins{ptMin}_{ptMax}")
    outFile.cd(f"cent_bins{centMinMax[0]}_{centMinMax[1]}/pt_bins{ptMin}_{ptMax}")
    hMassForFit.Write(f'hMass_pt{ptMin*10:.0f}_{ptMax*10:.0f}')
    fTotFuncMass = massFitters[iPt].GetMassFunc()
    print(f"fTotFuncMass.GetParameter(1): {fTotFuncMass.GetParameter(1)}")
    print(f"Integral of CB: {fTotFuncMass.Integral(0,4)}")
    print(f"Integral of CB normalized: {fTotFuncMass.Integral(0,4)/fTotFuncMass.GetParameter(1)}")
    SetObjectStyle(fTotFuncMass, color=kBlue, linestyle=2, linewidth=3)
    fTotFuncMass.Write(f'fTotFuncMass_pt{ptMin*10:.0f}_{ptMax*10:.0f}')

    # first parameter is global normalization
    hSgnMCFuncPars = TH1D('hSgnMCFuncPars', f';{ptTit};#chi^{{2}}/#it{{ndf}}', fTotFuncMass.GetNpar()-1, 0, fTotFuncMass.GetNpar()-1)
    for iBin in range(hSgnMCFuncPars.GetNbinsX()):
        hSgnMCFuncPars.SetBinContent(iBin+1, fTotFuncMass.GetParameter(iBin+1))
        hSgnMCFuncPars.SetBinError(iBin+1, fTotFuncMass.GetParError(iBin+1))
    hSgnMCFuncPars.Write()

    c1 = TCanvas(f"cMass_{ptMin*10:.0f}_{ptMax*10:.0f}", f"Mass Fit {ptMin}-{ptMax} GeV/c", 800, 600)
    hMassForFit.SetStats(0)
    hMassForFit.Draw("E")
    fTotFuncMass.SetLineColor(kBlue)
    fTotFuncMass.SetLineWidth(2)
    fTotFuncMass.Draw("same")
    c1.Write()
    c1.SaveAs(f"/home/mdicosta/FlowDplus/FinalResults/crystalball/MassFit_{ptMin*10:.0f}_{ptMax*10:.0f}.png")


#save output histos
outFile.cd()
hRawYields.Write()
hSigma.Write()
hMean.Write()
hRedChi2.Write()

outFile.Close()
