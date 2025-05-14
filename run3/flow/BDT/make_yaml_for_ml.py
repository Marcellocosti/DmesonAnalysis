'''
Script to create a yaml file with a set of cuts for ML
python3 make_yaml_for_ml.py config_flow.yml -o path/to/output -s text

'''
import yaml
import argparse
import os
import numpy as np
import sys
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get script's directory
sys.path.append(os.path.abspath(os.path.join(script_dir, '..')))  # Append parent directory

def pad_to_length(list, target_len):
    '''
        Function to pad a list to a target length
        Args:
            lst (list): list to be padded
            target_len (int): target length of the list
        Returns:
            list: padded list
    '''
    return list + [list[-1]] * (target_len - len(list)) if len(list) < target_len else list

def make_yaml(flow_config, outputdir, suffix):
    '''
        Function to create a yaml file with a set of cuts for ML
        Args:
            flow_config (str): path to the flow config file
            outputdir (str): path to the output directory
            suffix (str): suffix for the output files
    '''
    with open(flow_config, 'r') as f:
        cfg = yaml.safe_load(f)

    # os.makedirs(outputdir, exist_ok=True)
    ptmins = cfg['ptmins']
    ptmaxs = cfg['ptmaxs']
    massmins = cfg['MassMin']
    massmaxs = cfg['MassMax']
    if cfg.get('select_bin'):
        print(f"Using pt binning from the cfg file")
        npt_bins = 1
        ptmins = [ptmins[cfg['select_bin']-1]]
        ptmaxs = [ptmaxs[cfg['select_bin']-1]]
        massmins = [massmins[cfg['select_bin']-1]]
        massmaxs = [massmaxs[cfg['select_bin']-1]]
    else:
        if len(ptmins) != len(ptmaxs):
            raise ValueError(f'''The number of pt bins({len(ptmins)}, {len(ptmaxs)} are not the same''')
        npt_bins = len(ptmins)
    ptBinIdxs = [cfg['ptmins'].index(pt) for pt in ptmins]

    if cfg['minimisation'].get('correlated'):
        sig_cut = cfg['cut_variation']['corr_bdt_cut']['sig']
        sig_cuts_lower = [list(np.arange(sig_cut['min'][iPt], sig_cut['max'][iPt], sig_cut['step'][iPt])) for iPt in ptBinIdxs]
        sig_cuts_upper = [[1.0] * len(sig_low_edge) for sig_low_edge in sig_cuts_lower]
        nCutSets = [len(sig_low_edge) for sig_low_edge in sig_cuts_lower]
        bkg_cuts_upper = [[cfg['cut_variation']['corr_bdt_cut']['bkg_max'][idx]] * nCutSets[iPt] for iPt, idx in enumerate(ptBinIdxs)]
    else:
        sig_cut = cfg['cut_variation']['uncorr_bdt_cut']['sig']
        sig_cuts_lower = [sig_cut[iPt]['min'] for iPt in ptBinIdxs]
        sig_cuts_upper = [sig_cut[iPt]['max'] for iPt in ptBinIdxs]
        nCutSets = [len(sig_cuts_lower[iPt]) for iPt in ptBinIdxs]
        bkg_cuts_upper = [cfg['cut_variation']['uncorr_bdt_cut']['bkg_max'][iPt] for iPt in ptBinIdxs]

    for iPt in range(npt_bins):
        assert len(sig_cuts_lower[iPt]) == len(sig_cuts_upper[iPt]) == len(bkg_cuts_upper[iPt]) == nCutSets[iPt], (
            f"Mismatch in lengths for pt bin {iPt}: \n"
            f"sig_low:{len(sig_cuts_lower[iPt])}, \n"
            f"sig_up: {len(sig_cuts_upper[iPt])}, \n"
            f"bkg_up: {len(bkg_cuts_upper[iPt])}, \n"
            f"nCutSets: {nCutSets[iPt]}"
        )

    maxCutSets = max(nCutSets)

    sig_cuts_lower = [pad_to_length(cuts, maxCutSets) for cuts in sig_cuts_lower]
    sig_cuts_upper = [pad_to_length(cuts, maxCutSets) for cuts in sig_cuts_upper]
    bkg_cuts_upper = [pad_to_length(cuts, maxCutSets) for cuts in bkg_cuts_upper]

    os.makedirs(f'{outputdir}/config', exist_ok=True)
    for iCut in range(maxCutSets):
        score_bkg_max = [float(bkg_cuts_upper[i][iCut]) for i in range(len(ptBinIdxs))]
        score_fd_min  = [float(sig_cuts_lower[i][iCut]) for i in range(len(ptBinIdxs))]
        score_fd_max  = [float(sig_cuts_upper[i][iCut]) for i in range(len(ptBinIdxs))]

        combinations = {
            'icutset': iCut,
            'cutvars': {
                'Pt': {'min': ptmins, 'max': ptmaxs, 'name': 'pt_cand'},
                'score_bkg': {'min': [0.0] * len(ptmins), 'max': score_bkg_max, 'name': 'score_bkg'},
                'score_FD': {'min': score_fd_min, 'max': score_fd_max, 'name': 'score_FD'},
            },
            'fitrangemin': massmins,
            'fitrangemax': massmaxs,
        }

        with open(f'{outputdir}/config/cutset_{suffix}_{iCut:02}.yml', 'w') as file:
            yaml.dump(combinations, file, default_flow_style=False)

    print(f'Yaml files are saved in {outputdir}/config')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('flow_config', metavar='text', default='config_flow.yml')
    parser.add_argument('--preprocessed', action='store_true', help='Flag to indicate preprocessing of the sparses')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text", default="", help="suffix for output files")
    args = parser.parse_args()

    make_yaml(args.flow_config, args.outputdir, args.suffix)
