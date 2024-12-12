'''
python script with helper functions to load objects from task
'''

import sys
from ROOT import TFile, TAxis, THnSparseF, TH2F  # pylint: disable=import-error,no-name-in-module
import numpy as np

# pylint: disable=too-many-branches,too-many-statements, too-many-return-statements
def LoadSparseFromTask(infilename, inputCfg, particleName, no_List_isMC=False):
    '''
    Method to retrieve sparses from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file
    - flag to whether load 'isMC' and list from config
        - if True, load sparses for MC and no list
        - if False, load 'isMC' from config and load sparse from list
        - default is False

    Returns
    ----------
    - list of sparses with reconstructed quantities for all candidates,
      and prompt, FD D mesons and bkg candidates (only if MC)
    - list of sparses with generated quantities for prompt and FD D mesons (only if MC)
    '''
    print('Loading THnSparses from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None, None
    
    if not no_List_isMC:
        inlistData = indirData.Get(inputCfg['listname'])
        if not inlistData:
            print(f'List {inputCfg["listname"]} not found!')
            return None, None
        isMC = inputCfg['isMC'] if 'isMC' in inputCfg else False
    else:
        isMC = True

    sparses, sparsesGen = {}, {}
    if inputCfg['sparsenameAll']:
        if not no_List_isMC:
            sparses['RecoAll'] = inlistData.FindObject(inputCfg['sparsenameAll']) # not mandatory for MC
        else:
            sparses['RecoAll'] = indirData.Get(inputCfg['sparsenameAll'])
        if not sparses['RecoAll']:
            print(f'ERROR: sparse {inputCfg["sparsenameAll"]} not found!')
            return None, None
    if isMC:
            
        if ((inputCfg['sparsenamePrompt'] is not None and inputCfg['sparsenameAll'] == inputCfg['sparsenamePrompt']) or
                (inputCfg['sparsenameFD'] is not None and inputCfg['sparsenameAll'] == inputCfg['sparsenameFD'])):
            print('ERROR: do not use the same object for different sparses, this gives an error when merged! Exit')
            sys.exit()
        if inputCfg['sparsenamePrompt'] is not None:
            if not no_List_isMC:
                sparses['RecoPrompt'] = inlistData.FindObject(inputCfg['sparsenamePrompt'])
            else:
                sparses['RecoPrompt'] = indirData.Get(inputCfg['sparsenamePrompt'])
                sparses['RecoPrompt'].SetName("SparseRecoPrompt")
            if not sparses['RecoPrompt']:
                print(f'ERROR: sparse {inputCfg["sparsenamePrompt"]} not found!')
                return None, None
        if inputCfg['sparsenameFD'] is not None:
            if not no_List_isMC:
                sparses['RecoFD'] = inlistData.FindObject(inputCfg['sparsenameFD'])
            else:
                sparses['RecoFD'] = sparses['RecoPrompt'].Clone('SparseRecoFD')
            if not sparses['RecoFD']:
                print(f'ERROR: sparse {inputCfg["sparsenameFD"]} not found!')
                return None, None
        if particleName == "Dplus":
            print("GETTING TH2s FOR DPLUS")
            if not no_List_isMC:
                print("LIST")
                print(f"type(inlistData): {type(inlistData)}")
                GenPromptTh2 = inlistData.FindObject(inputCfg['TH2nameGenPrompt'])
                GenFDTh2 = inlistData.FindObject(inputCfg['TH2nameGenFD'])
                sparsesGen['GenPrompt'] = GetSparseFromTH2(GenPromptTh2)
                sparsesGen['GenFD'] = GetSparseFromTH2(GenFDTh2)
            else:
                print("NO LIST")
                # print(f"type(indirData): {type(indirData)}")
                # GenPromptTh2 = TH2F()
                # obj = indirData.Get('hPt')
                # print(f"type(obj): {type(obj)}")
                # obj.Copy(GenPromptTh2)
                # indirData.FindObject('hPt').Copy(GenPromptTh2)
                # indirData.FindObject('hPt').Copy(GenPromptTh2)
                GenPromptTh2 = indirData.Get(inputCfg['TH2nameGenPrompt'])
                GenFDTh2 = indirData.Get(inputCfg['TH2nameGenFD'])
                # print(f"type(GenPromptTh2): {type(GenPromptTh2)}")
                # print(f"inputCfg['TH2nameGenPrompt']: {inputCfg['TH2nameGenPrompt']}")
                # print(f"inputCfg['TH2nameGenFD']: {inputCfg['TH2nameGenFD']}")
                # GenFDTh2 = TH2F() 
                # indirData.FindObject(inputCfg['TH2nameGenFD']).Copy(GenFDTh2)
                # print(f"type(GenFDTh2): {type(GenFDTh2)}")
                # print(f"indirData.GetListOfKeys(): {indirData.GetListOfKeys()}")
                # print(f"indirData.ls(): {indirData.ls()}")
                sparsesGen['GenPrompt'] = GetSparseFromTH2(GenPromptTh2)
                sparsesGen['GenFD'] = GetSparseFromTH2(GenFDTh2)
                # sparsesGen['GenPrompt'] = indirData.Get(inputCfg['TH2nameGenPrompt'])
                # sparsesGen['GenPrompt'].SetName("TH2GenPrompt")
                # sparsesGen['GenFD'] = indirData.Get(inputCfg['TH2nameGenFD'])
                # sparsesGen['GenFD'].SetName("TH2GenFD")
            
            print('DONEEEE')
            print(f"type(sparsesGen['GenPrompt']): {type(sparsesGen['GenPrompt'])}")
            print(f"type(sparsesGen['GenFD']): {type(sparsesGen['GenFD'])}")
        else: 
            if not no_List_isMC:
                sparsesGen['GenPrompt'] = inlistData.FindObject(inputCfg['sparsenameGenPrompt'])
            else:
                sparsesGen['GenPrompt'] = indirData.Get(inputCfg['sparsenameGenPrompt'])
                sparsesGen['GenPrompt'].SetName("SparseGenPrompt")
            if not sparsesGen['GenPrompt']:
                print(f'ERROR: sparse {inputCfg["sparsenameGenPrompt"]} not found!')
                return None, None
            if not no_List_isMC:
                sparsesGen['GenFD'] = inlistData.FindObject(inputCfg['sparsenameGenFD'])
            else:
                sparsesGen['GenFD'] = sparsesGen['GenPrompt'].Clone('SparseGenFD')
            if not sparsesGen['GenFD']:
                print(f'ERROR: sparse {inputCfg["sparsenameGenFD"]} not found!')
                return None, None
        
        # For D0, the reflection and reconstruction are in the same sparse
        if inputCfg.get('enableRef'):
            if inputCfg['sparsenameRefl'] is not None:
                if not no_List_isMC:
                    sparses['RecoRefl'] = inlistData.FindObject(inputCfg['sparsenameRefl'])
                else:
                    sparses['RecoRefl'] = indirData.Get(inputCfg['sparsenameRefl'])
                    sparses['RecoRefl'].SetName("SparseRefl")
                    sparses['RecoReflPrompt'] = sparses['RecoRefl'].Clone('SparseReflPrompt')
                    sparses['RecoReflFD'] = sparses['RecoRefl'].Clone('SparseReflFD')
            if not sparses['RecoRefl']:
                print(f'ERROR: sparse {inputCfg["sparsenameRefl"]} not found!')
                return None, None

        if inputCfg.get('enableSecPeak'):
            if inputCfg['sparsenamePromptSecPeak'] is not None:
                if not no_List_isMC:
                    sparses['RecoSecPeakPrompt'] = inlistData.FindObject(inputCfg['sparsenamePromptSecPeak'])
                else:
                    sparses['RecoSecPeakPrompt'] = indirData.Get(inputCfg['sparsenamePromptSecPeak'])
                if not sparses['RecoSecPeakPrompt']:
                    print(f'ERROR: sparse {inputCfg["sparsenamePromptSecPeak"]} not found!')
                    return None, None
            if inputCfg['sparsenameFDSecPeak'] is not None:
                if not no_List_isMC:
                    sparses['RecoSecPeakFD'] = inlistData.FindObject(inputCfg['sparsenameFDSecPeak'])
                else:
                    sparses['RecoSecPeakFD'] = indirData.Get(inputCfg['sparsenameFDSecPeak'])
                if not sparses['RecoSecPeakFD']:
                    print(f'ERROR: sparse {inputCfg["sparsenameFDSecPeak"]} not found!')
                    return None, None
            if not no_List_isMC:
                sparsesGen['GenSecPeakPrompt'] = inlistData.FindObject(inputCfg['sparsenameGenPromptSecPeak'])
            else:
                sparsesGen['GenSecPeakPrompt'] = indirData.Get(inputCfg['sparsenameGenPromptSecPeak'])
            if not sparsesGen['GenSecPeakPrompt']:
                print(f'ERROR: sparse {inputCfg["sparsenameGenPromptSecPeak"]} not found!')
                return None, None
            if not no_List_isMC:
                sparsesGen['GenSecPeakFD'] = inlistData.FindObject(inputCfg['sparsenameGenFDSecPeak'])
            else:
                sparsesGen['GenSecPeakFD'] = indirData.Get(inputCfg['sparsenameGenFDSecPeak'])
            if not sparsesGen['GenSecPeakFD']:
                print(f'ERROR: sparse {inputCfg["sparsenameGenFDSecPeak"]} not found!')
                return None, None
    infileData.Close()

    return sparses, sparsesGen


def LoadSingleSparseFromTask(infilename, inputCfg, sparsetype='sparsenameBkg'):
    '''
    Method to retrieve single sparse from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file
    - sparse that should be returned

    Returns
    ----------
    - selected sparse from file
    '''
    print('Loading THnSparse from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None
    inlistData = indirData.Get(inputCfg['listname'])
    if not inlistData:
        print(f'List {inputCfg["listname"]} not found!')
        return None
    sparse = inlistData.FindObject(inputCfg[sparsetype])
    if not sparse:
        print(f'ERROR: sparse {inputCfg[sparsetype]} not found!')
        return None

    return sparse


def LoadNormObjFromTask(infilename, inputCfg):
    '''
    Method to retrieve normalisation objects from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - histo with event info and normalisation counter
    '''
    print('Loading norm objects from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None, None
    inlistData = indirData.Get(inputCfg['listname'])
    if not inlistData:
        print(f'List {inputCfg["listname"]} not found!')
        return None, None
    normCounter = indirData.Get(inputCfg['normname'])
    if not normCounter:
        print(f'Norm counter {inputCfg["normname"]} not found!')
        return None, None
    hEv = inlistData.FindObject(inputCfg['histoevname'])
    if not hEv:
        print(f'Histogram {inputCfg["histoevname"]} not found!')
        return None, None

    return hEv, normCounter


def LoadListFromTask(infilename, inputCfg):
    '''
    Method to retrieve list of objects from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - TList of objects from input file
    '''
    print('Loading TList from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None
    inlistData = indirData.Get(inputCfg['listname'])
    if not inlistData:
        print(f'List {inputCfg["listname"]} not found!')
        return None

    return inlistData


def LoadCutObjFromTask(infilename, inputCfg):
    '''
    Method to retrieve D-meson cut object from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - D-meson cut object
    '''
    print('Loading cut object from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None, None
    cutobjname = inputCfg['listname'].replace('coutputDs', 'coutputDsCuts')
    cutobjname = cutobjname.replace('coutputDplus', 'coutputDplusCuts')
    cutobj = indirData.Get(cutobjname)
    if not cutobj:
        print(f'Cut object {cutobjname} not found!')
        return None, None

    return cutobj, cutobjname


def LoadPIDSparses(infilename, inputCfg):
    '''
    Method to retrieve PID sparses from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - sparse of NsigmaTPC and NsigmaTOF variables
    - sparse of NsigmaComb variables
    - dictionary with variable : sparse-axis number
    '''
    print('Loading PID THnSparses from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    if not indirData:
        print(f'Directory {inputCfg["dirname"]} not found!')
        return None, None, None
    inlistData = indirData.Get(inputCfg['listname'])
    if not inlistData:
        print(f'List {inputCfg["listname"]} not found!')
        return None, None, None
    sparsePIDNsigma = inlistData.FindObject('fnSparsePID')
    if not sparsePIDNsigma:
        print('ERROR: sparse fnSparsePID not found!')
        return None, None, None
    sparsePIDNsigmaComb = inlistData.FindObject('fnSparsePIDcomb')
    if not sparsePIDNsigma:
        print('ERROR: sparse fnSparsePIDcomb not found!')
        return None, None, None

    # dictionary of sparse axes with detectors, mass hypothesis, daughter number
    axes = {'TPC': {'Pi': {'0': 2, '1': 6, '2': 10}, 'K': {'0': 3, '1': 7, '2': 11}},
            'TOF': {'Pi': {'0': 4, '1': 8, '2': 12}, 'K': {'0': 5, '1': 9, '2': 13}},
            'Comb': {'Pi': {'0': 2, '1': 4, '2': 6}, 'K': {'0': 3, '1': 5, '2': 7}}}

    return sparsePIDNsigma, sparsePIDNsigmaComb, axes


def LoadSparseFromTaskV2(inputCfg):
    '''
    Method to retrieve sparses from D-meson vn output task file

    Inputs
    ----------
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - list of sparses with reconstructed quantities for all candidates
    '''
    infilename = inputCfg['filename']
    print(f'Loading THnSparses from file {infilename}')
    sparses = []
    infileData = TFile(infilename)

    for dirname, listname in zip(inputCfg['dirname'], inputCfg['listname']):
        indirData = infileData.Get(dirname)
        if not indirData:
            print(f'Directory {dirname} not found!')
            return []
        inlistData = indirData.Get(listname)
        if not inlistData:
            print(f'List {listname} not found!')
            return []
        sparse = inlistData.FindObject(inputCfg['sparsename'])
        if not sparse:
            print('Sparse not found!')
            return []
        sparses.append(sparse)

    return sparses


def LoadListFromTaskV2(infilename, dirname, listname):
    '''
    Method to retrieve list of objects from D-meson vn output task file

    Inputs
    ----------
    - input root file name
    - name of directory in input root file
    - name of list in directory in input root file

    Returns
    ----------
    - TList of objects from input file
    '''
    print('Loading TList from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(dirname)
    if not indirData:
        print(f'Directory {dirname} not found!')
        return None
    inlistData = indirData.Get(listname)
    if not inlistData:
        print(f'List {listname} not found!')
        return None

    return inlistData

def GetSparseFromTH2(th2histo):
    '''
    Method to convert the pt vs Y TH2 histo of the D+ task
    to a ThNSparse having 3 axes (pt, pt_B, Y), for consistency
    with the D0 task

    Inputs
    ----------
    - input histogram of type TH2

    Returns
    ----------
    - ThNSparse
    '''
    
    # print('CONVERTING SPARSE FROM TH2')
    # print(f"TYPE TH2HISTO: {type(th2histo)}")
    # axisPt = TAxis() 
    # th2histo.GetXaxis().Copy(axisPt)
    # axisY = TAxis() 
    # th2histo.GetYaxis().Copy(axisY)
    # axisPtB = TAxis(1, 0, 1)
    # axesVector = [axisPt, axisPtB, axisY]
    # print(f"axesVector: {axesVector}")
    # sparseFromTh2 = THnSparseF(f"{th2histo.GetTitle()}", f"{th2histo.GetName()}", axesVector)
    # for iptbin in axisPt.GetNbins():
    #     ptBinCenter = axisPt.GetBinCenter(iptbin)
    #     for iYbin in axisY.GetNbins():
    #         YBinCenter = axisY.GetBinCenter(iYbin)
    #         bin = sparseFromTh2.GetBin(np.as_array([ptBinCenter, 0.5, YBinCenter]), True)
    #         sparseFromTh2.SetBinContent(th2histo.GetBinContent(iptbin, iYbin))
    
    # print('ENDED CONVERSION')
    # return sparseFromTh2 
    
    # Get bin edges for the X and Y axes
    n_bins_x = th2histo.GetNbinsX()
    n_bins_y = th2histo.GetNbinsY()
    x_edges = [th2histo.GetXaxis().GetBinLowEdge(i) for i in range(1, n_bins_x + 2)]
    y_edges = [th2histo.GetYaxis().GetBinLowEdge(i) for i in range(1, n_bins_y + 2)]

    # Create nbins, xmin, and xmax arrays
    nbins = np.array([n_bins_x, 1, n_bins_y], dtype=np.int32)
    xmin = np.array([min(x_edges), 0.0, min(y_edges)], dtype=np.float64)
    xmax = np.array([max(x_edges), 1.0, max(y_edges)], dtype=np.float64)

    # Define the THnSparse
    sparse = THnSparseF(
        th2histo.GetName(),
        th2histo.GetTitle(),
        3,  # Number of dimensions
        nbins,  # Number of bins per dimension
        xmin,  # Lower bounds
        xmax   # Upper bounds
    )

    # Fill the sparse histogram with content from the TH2
    for ix in range(1, n_bins_x + 1):
        for iy in range(1, n_bins_y + 1):
            content = th2histo.GetBinContent(ix, iy)
            bin_center_x = th2histo.GetXaxis().GetBinCenter(ix)
            bin_center_y = th2histo.GetYaxis().GetBinCenter(iy)
            sparse.Fill(np.array([bin_center_x, 0.5, bin_center_y], dtype=np.float64), content)

    print('ENDED CONVERSION')
    return sparse