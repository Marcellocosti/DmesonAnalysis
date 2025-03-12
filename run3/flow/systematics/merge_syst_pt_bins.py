import os
import argparse
from ROOT import TFile, TH1D, TCanvas, kAzure, kOrange
import array
import sys
sys.path.append('../../../')
import fitz
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle

def merge_syst_info(syst_outdir):
    pt_bins_directories = [d for d in os.listdir(syst_outdir) if os.path.isdir(os.path.join(syst_outdir, d))]
    print(f"pt_bins_directories: {pt_bins_directories}")
    pt_bins_directories = sorted(pt_bins_directories, key=lambda x: (int(x.split('_')[1]), int(x.split('_')[2])))
    print(f"pt_bins_directories: {pt_bins_directories}")
    pt_bins = array.array('d', sorted(set(float(x)/10 for d in pt_bins_directories for x in d.split('_')[1:])))
    print(f"pt_bins: {pt_bins}")
    
    hvn_syst_histos = {}
    hvn_syst_histos_summary = {}
    histos_names = ['hsyst_final_prompt', 'hsyst_final_prompt_rel', 'hsyst_final_fd', 'hsyst_final_fd_rel']
    histos_axes_labels = ['Syst. Unc. Prompt D^{+}', 'Rel. Syst. Unc. Prompt D^{+} (%)', 'Syst. Unc. NonPrompt D^{+}', 'Rel. Syst. Unc. NonPrompt D^{+} (%)']
    hvn_syst_histos_ntrials = TH1D(f"hist_ntrials", f";p_{{T}};# trials", len(pt_bins)-1, pt_bins)
    for histo_name, histos_axes_label in zip(histos_names, histos_axes_labels):
        hvn_syst_histos[histo_name] = []
        hvn_syst_histos_summary[histo_name] = TH1D(f"hist_{histo_name}_over_signal", f";p_{{T}};{histos_axes_label}", len(pt_bins)-1, pt_bins)

    for ipt, pt_bin_dir in enumerate(pt_bins_directories):
        syst_file = TFile.Open(f"{syst_outdir}/{pt_bin_dir}/syst_summary/syst_v2.root", 'r')
        hvn_syst_histos_ntrials.SetBinContent(ipt+1, syst_file.Get('hvn_prompt').GetNbinsX())
        for histo_name in histos_names:
            hvn_syst_histos_summary[histo_name].SetBinContent(ipt+1, syst_file.Get(histo_name).GetBinContent(1))
            
            hvn_syst_histos[histo_name].append(syst_file.Get(histo_name))
            syst_file.Get(histo_name).SetDirectory(0)
        syst_file.Close()

    canvsyst_prompt = TCanvas('canvsyst_prompt', 'canvsyst_prompt', 800, 800)
    canvsyst_prompt.cd()
    canvsyst_prompt.SetGrid()
    SetObjectStyle(hvn_syst_histos_summary['hsyst_final_prompt'], markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hvn_syst_histos_summary['hsyst_final_prompt'].SetStats(0)
    hvn_syst_histos_summary['hsyst_final_prompt'].Draw('same hist')
    
    canvsyst_fd = TCanvas('canvsyst_fd', 'canvsyst_fd', 800, 800)
    canvsyst_fd.cd()
    canvsyst_fd.SetGrid()
    SetObjectStyle(hvn_syst_histos_summary['hsyst_final_fd'], markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hvn_syst_histos_summary['hsyst_final_fd'].SetStats(0)
    hvn_syst_histos_summary['hsyst_final_fd'].Draw('same hist')

    canvsyst_prompt_relunc = TCanvas('canvsyst_prompt_relunc', 'canvsyst_prompt_relunc', 800, 800)
    canvsyst_prompt_relunc.cd()
    canvsyst_prompt_relunc.SetGrid()
    SetObjectStyle(hvn_syst_histos_summary['hsyst_final_prompt_rel'], markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hvn_syst_histos_summary['hsyst_final_prompt_rel'].SetStats(0)
    hvn_syst_histos_summary['hsyst_final_prompt_rel'].Draw('same hist')
    
    canvsyst_fd_relunc = TCanvas('canvsyst_fd_relunc', 'canvsyst_fd_relunc', 800, 800)
    canvsyst_fd_relunc.cd()
    canvsyst_fd_relunc.SetGrid()
    SetObjectStyle(hvn_syst_histos_summary['hsyst_final_fd_rel'], markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hvn_syst_histos_summary['hsyst_final_fd_rel'].SetStats(0)
    hvn_syst_histos_summary['hsyst_final_fd_rel'].Draw('same hist')

    canvsyst_ntrials = TCanvas('canvsyst_ntrials', 'canvsyst_ntrials', 800, 800)
    canvsyst_ntrials.cd()
    canvsyst_ntrials.SetGrid()
    SetObjectStyle(hvn_syst_histos_ntrials, markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hvn_syst_histos_ntrials.SetStats(0)
    hvn_syst_histos_ntrials.Draw('same hist')

    canvsyst_prompt.SetLeftMargin(0.15)
    canvsyst_prompt.SaveAs(f'{syst_outdir}/SystPromptv2_vs_pt.pdf')
    canvsyst_prompt_relunc.SetLeftMargin(0.15)
    canvsyst_prompt_relunc.SaveAs(f'{syst_outdir}/SystPromptv2_vs_pt_relunc.pdf')
    canvsyst_fd.SetLeftMargin(0.15)
    canvsyst_fd.SaveAs(f'{syst_outdir}/SystFDv2_vs_pt.pdf')
    canvsyst_fd_relunc.SetLeftMargin(0.15)
    canvsyst_fd_relunc.SaveAs(f'{syst_outdir}/SystFDv2_vs_pt_relunc.pdf')
    canvsyst_ntrials.SetLeftMargin(0.15)
    canvsyst_ntrials.SaveAs(f'{syst_outdir}/Systv2_ntrials.pdf')

    print(f"Writing to: {syst_outdir}/SystAllPt.root")
    outfile = TFile(f"{syst_outdir}/SystAllPt.root", "recreate")
    hvn_syst_histos_ntrials.Write()
    for histo_summary in hvn_syst_histos_summary.values():
        histo_summary.Write()
    outfile.Close()

    trials_summary_all_ptbins = fitz.open()
    for ipt, pt_bin_dir in enumerate(pt_bins_directories):
        doc = fitz.open(f"{syst_outdir}/{pt_bin_dir}/syst_summary/Systv2.pdf")
        trials_summary_all_ptbins.insert_pdf(doc)
    print(f"Writing to: {syst_outdir}/SystTrialsAllPt.pdf")
    trials_summary_all_ptbins.save(f"{syst_outdir}/SystTrialsAllPt.pdf")
    trials_summary_all_ptbins.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('syst_outdir', metavar='text', default='path/to/syst/output')
    args = parser.parse_args()

    merge_syst_info(args.syst_outdir)
