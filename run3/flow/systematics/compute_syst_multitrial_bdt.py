import sys
import os
import numpy as np
import argparse
import re
import glob
from ROOT import TFile, TCanvas, TH1F, TGraphAsymmErrors, TLegend, kOrange, kAzure, kBlack
sys.path.append('../../../')
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle

class SystMultitrial:
    prompt = False
    ptmin = 0
    ptmax = 0

class Trails:
    def __init__(self, trail_path, default=False):
        self.trail_path = trail_path
        self.name = trail_path.split('_')[-1]
        self.v2vsfracfile = self._get_v2_vs_frac_file(trail_path)
        self.ryfiles = self._get_ry_files(trail_path)
        self.hv2vsptprompt = self._get_hv2vsptprompt(self.v2vsfracfile)
        self.hv2vsptfd = self._get_hv2vsptfd(self.v2vsfracfile)
        self.gvnsimfits, self.hchi2s, self.hsignificances = self._get_fit_results(self.ryfiles, default)
        # self.gvnsimfits = self._get_gvnSimFit(self.ryfiles)
        if default:
            self.ncut = len(self.ryfiles)

    def _get_pt_bins(self):
        pt_axis = self.hv2vsptprompt.GetXaxis()
        npt_bins = self.hv2vsptprompt.GetNbinsX()
        return [pt_axis.GetBinLowEdge(i) for i in range(1, npt_bins + 2)]
    
    def _get_v2_vs_frac_file(self, trail_path):
        # print(f"\n")
        print(f"Looking for V2VsFrac file in {trail_path}")
        try:
            # v2vsfrac_files = glob.glob(f'{trail_path}/**/V2VsFrac*.root', recursive=True)
            # print(f"v2vsfrac_files: {v2vsfrac_files}")
            for file in os.listdir(path=f'{trail_path}/V2VsFrac'):
                if file.endswith('.root'):
                    print(f"Found V2VsFrac file: {file}")
                    return os.path.join(trail_path, 'V2VsFrac', file)
            return v2vsfrac_files
        except FileNotFoundError:
            print(f'No V2VsFrac directory found in {trail_path}')
            self.v2vsfracfile = None
            return None
        
    def _get_ry_files(self, trail_path):
        # print(f"\n")
        print(f"Looking for ry file in {trail_path}")
        files = []
        try:
            # ry_files = glob.glob(f'{trail_path}/**/ry*.root', recursive=True)
            # ry_files.sort()
            # ry_files = glob.glob(f'{trail_path}/**/ry*.root', recursive=True)
            
            for raw_file in os.listdir(path=f'{trail_path}/ry'):
                if raw_file.endswith('.root'):
                    print(f"Found ry file: {raw_file}")
                    files.append(os.path.join(trail_path, 'ry', raw_file))
                # make sure the files are sorted
                files.sort()
        except FileNotFoundError:
            print(f'No ry directory found in {trail_path}')
            return []
        return files
    
    def _get_hv2vsptprompt(self, v2vsfracfile):# -> None | Any:
        if v2vsfracfile is None:
            return None

        print(f"\n\n")
        print(f"v2vsfracfile: {v2vsfracfile}")
        tfile = TFile().Open(v2vsfracfile)
        if not tfile.GetListOfKeys().Contains('hV2VsPtPrompt'):
            print(f'Warning: {v2vsfracfile} does not contain hV2VsPtPrompt')
            return None
        hv2vsptprompt = tfile.Get('hV2VsPtPrompt')
        hv2vsptprompt.SetDirectory(0)
        tfile.Close()
        return hv2vsptprompt
    
    def _get_hv2vsptfd(self, v2vsfracfile):# -> None | Any:
        if v2vsfracfile is None:
            return None
        tfile = TFile().Open(v2vsfracfile)
        if not tfile.GetListOfKeys().Contains('hV2VsPtFD'):
            print(f'Warning: {v2vsfracfile} does not contain hV2VsPtFD')
            return None
        hv2vsptfd = tfile.Get('hV2VsPtFD')
        hv2vsptfd.SetDirectory(0)
        tfile.Close()
        return hv2vsptfd
    
    def _get_fit_results(self, ryfiles, default):
        gvnsimfits = []
        hchi2s = []
        hsignificances = []
        for file in ryfiles:
            tfile = TFile().Open(file)
            if not tfile.GetListOfKeys().Contains('gvnSimFit'):
                print(f'Warning: {file} does not contain gvnSimFit')
                gvnsimfits.append(None)
                hchi2s.append(None)
                hsignificances.append(None)
                tfile.Close()
                continue
            gvnsimfits.append(tfile.Get('gvnSimFit'))
            hchi2s.append(tfile.Get('hRedChi2SimFit'))
            hchi2s[-1].SetDirectory(0)
            hsignificances.append(tfile.Get('hRawYieldsSignificanceSimFit'))
            hsignificances[-1].SetDirectory(0)
            tfile.Close()
        return gvnsimfits, hchi2s, hsignificances

SetGlobalStyle(titleoffsety=1.1, maxdigits=3, topmargin=0.1, bottommargin=0.4, leftmargin=0.3, rightmargin=0.15,
               labelsizey=0.04, setoptstat=0, setopttitle=0, setdecimals=True,titleoffsetx=0.74)

def compute_syst_multitrial(rypathsyst, ry_default, outputdir):
    pt_dirs = [d for d in os.listdir(rypathsyst) if os.path.isdir(os.path.join(rypathsyst, d)) and re.match(r'^pt_\d{2}_\d{2}$', d)]
    # default trail i.e. central value
    default_trail = Trails(ry_default, default=True)
    pt_bins = default_trail._get_pt_bins()
    for pt_dir in pt_dirs:
        trails_dirs = os.listdir(f"{rypathsyst}/{pt_dir}/trails/")
        trails = []
        for trail_dir in trails_dirs:
            trails.append(Trails(os.path.join(f"{rypathsyst}/{pt_dir}/trails/", trail_dir)))    

        print(f"pt_bins: {pt_bins}")
        
        folder_pt_value = float(re.search(r'pt_(\d+)_\d+', pt_dir).group(1)) / 10  # Extract first number
        default_pt_index = pt_bins.index(folder_pt_value) 
        print(f"folder_pt_value: {folder_pt_value}")
        print(f"default_pt_index: {default_pt_index}")        
        
        outputdir = os.path.join(outputdir, 'syst_multitrial')
        os.makedirs(outputdir, exist_ok=True)
        
        ntrials_syst = len(trails[0].ryfiles)
        print(f"ntrials_syst: {ntrials_syst}")
        if SystMultitrial.prompt:
            for iFile in range(ntrials_syst):
                print(f"iFile: {iFile}")
                compute_systematics(f"{rypathsyst}/{pt_dir}/syst_summary/", default_trail.gvnsimfits[iFile], \
                                    [trail.gvnsimfits[iFile] for trail in trails if len(trail.gvnsimfits) > iFile], \
                                    [trail.hchi2s[iFile] for trail in trails if len(trail.hchi2s) > iFile], \
                                    [trail.hsignificances[iFile] for trail in trails if len(trail.hsignificances) > iFile], \
                                    default_pt_index, trails, iFile)
            compute_systematics_prompt(f"{rypathsyst}/{pt_dir}/syst_summary/", default_trail, trails, default_pt_index)
        else:
            compute_systematics(f"{rypathsyst}/{pt_dir}/syst_summary/", default_trail.gvnsimfits, \
                    [trail.gvnsimfits for trail in trails], \
                            [trail.hchi2s for trail in trails], \
                    [trail.hsignificances for trail in trails], default_pt_index)

def compute_systematics(outputdir, gvn_vs_mass_default, gvn_vs_mass_trials, hchi2_trails, hsignificance_trails, iPtBin, trails, iFile=0):
    print("ENTERED IN compute_systematics")
    hsyst_final = hchi2_trails[0].Clone('hsyst_final')
    hsyst_final.SetTitle('hsyst_final;#it{p}_{T} (GeV/#it{c});syst (%)')
    hsyst_final.Reset()

    hvn = TH1F('hvn', 'hvn;trial;#it{v}_{n}', len(trails), 0, len(trails)+1)
    hsignificance_vs_trial = TH1F('hsignificance_vs_trial', 'hsignificance_vs_trial;trial;significance', len(trails), 0, len(trails)+1)
    hchi2_vs_trial = TH1F('hchi2_vs_trial', 'hchi2_vs_trial;trial;#chi^{2}', len(trails), 0, len(trails)+1)
    hsyst = TH1F('hsyst', 'hsyst;#it{v}_{n}^{trial} - #it{v}_{n}^{default};"', 200, -0.2, 0.2)
    canvas = TCanvas(f'syst_multitrail_pt{iPtBin}', f'syst_multitrail_pt{iPtBin}_{iPtBin}', 800, 800)
    canvas.cd().Divide(2, 2)
        
    SetObjectStyle(hvn, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(hsignificance_vs_trial, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(hchi2_vs_trial, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(hsyst, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    
    for itrial, (trial, gvn, hchi2, hsign) in enumerate(zip(trails, gvn_vs_mass_trials, hchi2_trails, hsignificance_trails)):    # loop over trials
        if trial.gvnsimfits == None:
            print(f'No prompt v2 vs pt found for trial {j}: {trial.trail_path}')
            continue
        significance = hsign.GetBinContent(iPtBin+1)
        chi2 = hchi2.GetBinContent(iPtBin+1)

        if len(trial.ryfiles) <= iFile:
            print(f'No ry file found for trial {j}: {trial.trail_path}')
            continue
        # Skip chi2 higher than 5 and significance lower than 3
        if (chi2 > 4 and chi2 != 0) or (significance < 6 or significance > 600):
            print(f'Skipping trial {itrial}: {trial.ryfiles[iFile]} for pt bin {iPtBin} with chi2 = {chi2} and significance = {significance}')
            continue
        hchi2_vs_trial.SetBinContent(itrial+1, chi2)
        hchi2_vs_trial.SetBinError(itrial+1, hchi2.GetBinError(iPtBin+1))
        hsignificance_vs_trial.SetBinContent(itrial+1, significance)
        hsignificance_vs_trial.SetBinError(itrial+1, hsign.GetBinError(iPtBin+1))

        hvn.SetBinContent(itrial+1, gvn.GetY()[0])
        hvn.SetBinError(itrial+1, gvn.GetEYlow()[0])
        
        # the gvn of trials only have one point
        hsyst.Fill(gvn.GetY()[0] - gvn_vs_mass_default.GetY()[iPtBin])
        
        # Compute systematic uncertainty
        rms = hsyst.GetRMS()
        mean = hsyst.GetMean()
        syst = np.sqrt(rms**2 + mean**2)
        print(f'pt bin {iPtBin}, mean = {mean}, rms = {rms}, syst = {syst}')
        gsyst = TGraphAsymmErrors()
        gsyst.SetPoint(0, 0, hsyst.GetMaximum()*0.5)
        gsyst.SetPointError(0, syst, syst,
                                hsyst.GetMaximum()*0.5, hsyst.GetMaximum()*0.5)
        SetObjectStyle(gsyst, markerstyle=20, markercolor=kOrange+2,
                       markersize=1, linecolor=kOrange+2,
                          linewidth=2, fillcolor=kOrange+2, fillstyle=3153,
                          fillalpha=0.5, linestyle=9)

        # Pad 1: vn vs trial
        canvas.cd(1).SetGrid()
        # Define reference line
        gref = TGraphAsymmErrors()
        gref.SetPoint(0, 0, gvn_vs_mass_default.GetY()[iPtBin])
        gref.SetPointError(0, 0, 0, gvn_vs_mass_default.GetEYlow()[iPtBin],
                           gvn_vs_mass_default.GetEYhigh()[iPtBin])
        gref.SetPoint(1, len(trails), gvn_vs_mass_default.GetY()[iPtBin])
        gref.SetPointError(1, 0, 0, gvn_vs_mass_default.GetEYlow()[iPtBin],
                           gvn_vs_mass_default.GetEYhigh()[iPtBin])
        SetObjectStyle(gref, markerstyle=20, markercolor=kAzure+2,
                       markersize=0, linecolor=kAzure+2,
                       linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        hvn.Draw('same')
        gref.Draw('c3 same')
        canvas.cd(2)

        # Legend
        leg = TLegend(0.6, 0.6, 0.9, 0.9)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        leg.SetHeader(f'ptmin < #it{{p}}_{{T}} < ptmax GeV/#it{{c}}')
        leg.AddEntry(hvn, f'vn vs trial ({hvn.GetEntries()} trials)', 'p')
        leg.AddEntry(gsyst, f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst:.3f}', 'f')
        leg.AddEntry(gref, 'vn stat. unc.', 'f')

        # Define reference vertical line at 1
        gref_two = gref.Clone()
        gref_two.SetPoint(0, 0, hsyst.GetMaximum()*0.5)
        gref_two.SetPointError(0, gvn_vs_mass_default.GetEYlow()[iPtBin],
                                   gvn_vs_mass_default.GetEYhigh()[iPtBin],
                                   hsyst.GetMaximum()*0.5, hsyst.GetMaximum()*0.5)
        SetObjectStyle(gref_two, markerstyle=20, markercolor=kAzure+2,
                       markersize=1, linecolor=kAzure+2,
                       linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        hsyst.GetYaxis().SetRangeUser(0, hsyst.GetMaximum()*1.8)
        hsyst.GetXaxis().SetRangeUser(hsyst.GetXaxis().GetXmin(), hsyst.GetXaxis().GetXmax())
        hsyst.Draw('same')
        gsyst.Draw('2')
        gref_two.Draw('2')
        gsyst.Draw('2')
        leg.Draw()
        # Pad 3: chi2 vs trial
        canvas.cd(3)
        hchi2_vs_trial.Draw('same')
        # Pad 4: significance vs trial
        canvas.cd(4)
        hsignificance_vs_trial.Draw('same')
        hsyst_final.SetBinContent(1, syst)
        canvas.Update()

    canvsyst = TCanvas('canvsyst', 'canvsyst', 800, 800)
    canvsyst.cd()
    canvsyst.SetGrid()
    SetObjectStyle(hsyst_final, markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hsyst_final.Draw('same e')

    canvas.SaveAs(f'{outputdir}/SystRy_{iFile}.pdf')
    canvsyst.SaveAs(f'{outputdir}/SystRy_vs_pt_{iFile}.pdf)')

    outdir = os.path.join(outputdir, 'syst_multitrial')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfilename = os.path.join(outdir, f'syst_ry_{iFile}.root')
    outfile = TFile(outfilename, 'recreate')
    hvn.Write()
    hsyst.Write()
    hsyst_final.Write()
    outfile.Close()

def compute_systematics_prompt(outputdir, default_trail, trails, iPtBin=0):

    hsyst_final_prompt = trails[0].hchi2s[0].Clone('hsyst_final_prompt')
    hsyst_final_prompt_rel = trails[0].hchi2s[0].Clone('hsyst_final_prompt_rel')
    hsyst_final_fd = trails[0].hchi2s[0].Clone('hsyst_final_fd')
    hsyst_final_fd_rel = trails[0].hchi2s[0].Clone('hsyst_final_fd_rel')
    hsyst_final_prompt.SetTitle('hsyst_final_prompt;#it{p}_{T} (GeV/#it{c});syst (abs)')
    hsyst_final_prompt_rel.SetTitle('hsyst_final_prompt_rel;#it{p}_{T} (GeV/#it{c});syst (%)')
    hsyst_final_fd.SetTitle('hsyst_final_fd;#it{p}_{T} (GeV/#it{c});syst (abs)')
    hsyst_final_fd_rel.SetTitle('hsyst_final_fd_rel;#it{p}_{T} (GeV/#it{c});syst (%)')
    hsyst_final_prompt.Reset()
    hsyst_final_prompt_rel.Reset()
    hsyst_final_fd.Reset()
    hsyst_final_fd_rel.Reset()
    
    hvn_prompt = TH1F('hvn_prompt', 'hvn_prompt;trial;prompt #it{v}_{2}', len(trails), 0, len(trails)+1)
    hvn_fd = TH1F('hvn_fd', 'hvn_fd;trial;non-prompt #it{v}_{2}', len(trails), 0, len(trails)+1)
    hsyst_prompt = TH1F('hsyst_prompt', 'hsyst_prompt;prompt #it{v}_{2}^{trial} - prompt #it{v}_{2}^{default};"', 200, -0.2, 0.2)
    hsyst_fd = TH1F('hsyst_fd', 'hsyst_fd;non-prompt #it{v}_{2}^{trial} - non-prompt #it{v}_{2}^{default};"', 200, -0.2, 0.2)
    canvas = TCanvas(f'syst_multitrail_pt{iPtBin}', f'syst_multitrail_pt{iPtBin}', 800, 800)
    canvas.cd().Divide(2, 2)
        
    SetObjectStyle(hvn_prompt, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(hvn_fd, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(hsyst_prompt, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(hsyst_fd, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
        
    hvn_prompt, hsyst_prompt = process_trial(default_trail, trails, hvn_prompt, hsyst_prompt, iPtBin, prompt=True)
    hvn_fd, hsyst_fd = process_trial(default_trail, trails, hvn_fd, hsyst_fd, iPtBin, prompt=False)

    # Compute systematic uncertainty
    rms_prompt = hsyst_prompt.GetRMS()
    mean_prompt = hsyst_prompt.GetMean()
    syst_prompt = np.sqrt(rms_prompt**2 + mean_prompt**2)
    print(f'pt bin {iPtBin}, mean = {mean_prompt}, rms = {rms_prompt}, syst = {syst_prompt}')
    gsyst_prompt = TGraphAsymmErrors()
    gsyst_prompt.SetPoint(0, 0, hsyst_prompt.GetMaximum()*0.5)
    gsyst_prompt.SetPointError(0, syst_prompt, syst_prompt,
                            hsyst_prompt.GetMaximum()*0.5, hsyst_prompt.GetMaximum()*0.5)
    SetObjectStyle(gsyst_prompt, markerstyle=20, markercolor=kOrange+2, markersize=1, linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3153, fillalpha=0.5, linestyle=9)

    rms_fd = hsyst_fd.GetRMS()
    mean_fd = hsyst_fd.GetMean()
    syst_fd = np.sqrt(rms_fd**2 + mean_fd**2)
    print(f'pt bin {iPtBin}, mean = {mean_fd}, rms = {rms_fd}, syst = {syst_fd}')
    gsyst_fd = TGraphAsymmErrors()
    gsyst_fd.SetPoint(0, 0, hsyst_fd.GetMaximum()*0.5)
    gsyst_fd.SetPointError(0, syst_fd, syst_fd,
                            hsyst_fd.GetMaximum()*0.5, hsyst_fd.GetMaximum()*0.5)
    SetObjectStyle(gsyst_fd, markerstyle=20, markercolor=kOrange+2, markersize=1, linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3153, fillalpha=0.5, linestyle=9)
    
    # Pad 1: prompt v2 vs trial
    canvas.cd(1).SetGrid()
    # Define reference line
    gref_prompt = TGraphAsymmErrors()
    gref_prompt.SetPoint(0, 0, default_trail.hv2vsptprompt.GetBinContent(iPtBin+1))
    gref_prompt.SetPointError(0, 0, 0, default_trail.hv2vsptprompt.GetBinError(iPtBin+1),
                       default_trail.hv2vsptprompt.GetBinError(iPtBin+1))
    gref_prompt.SetPoint(1, len(trails), default_trail.hv2vsptprompt.GetBinContent(iPtBin+1))
    gref_prompt.SetPointError(1, 0, 0, default_trail.hv2vsptprompt.GetBinError(iPtBin+1),
                       default_trail.hv2vsptprompt.GetBinError(iPtBin+1))
    SetObjectStyle(gref_prompt, markerstyle=20, markercolor=kAzure+2,
                   markersize=0, linecolor=kAzure+2,
                   linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    hvn_prompt.Draw('same')
    gref_prompt.Draw('c3 same')
        
    # Pad 2: non-prompt v2 vs trial
    canvas.cd(2)
    # Define reference line
    gref_fd = TGraphAsymmErrors()
    gref_fd.SetPoint(0, 0, default_trail.hv2vsptfd.GetBinContent(iPtBin+1))
    gref_fd.SetPointError(0, 0, 0, default_trail.hv2vsptfd.GetBinError(iPtBin+1),
                       default_trail.hv2vsptfd.GetBinError(iPtBin+1))
    gref_fd.SetPoint(1, len(trails), default_trail.hv2vsptfd.GetBinContent(iPtBin+1))
    gref_fd.SetPointError(1, 0, 0, default_trail.hv2vsptfd.GetBinError(iPtBin+1),
                       default_trail.hv2vsptfd.GetBinError(iPtBin+1))
    SetObjectStyle(gref_fd, markerstyle=20, markercolor=kAzure+2,
                   markersize=0, linecolor=kAzure+2,
                   linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    hvn_fd.Draw('same')
    gref_fd.Draw('c3 same')
        
    # Pad 3: prompt v2 systematic uncertainty
    canvas.cd(3)
    leg_prompt = TLegend(0.6, 0.6, 0.9, 0.9)
    leg_prompt.SetBorderSize(0)
    leg_prompt.SetFillStyle(0)
    leg_prompt.SetTextSize(0.03)
    leg_prompt.SetHeader(f'ptmin < #it{{p}}_{{T}} < ptmax')
    leg_prompt.AddEntry(hvn_prompt, f'prompt v2 vs trial ({hvn_prompt.GetEntries()} trials)', 'p')
    leg_prompt.AddEntry(gsyst_prompt, f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst_prompt:.3f}', 'f')
    leg_prompt.AddEntry(gref_prompt, 'prompt v2 stat. unc.', 'f')
        
    gref_two_prompt = gref_prompt.Clone()
    gref_two_prompt.SetPoint(0, 0, hsyst_prompt.GetMaximum()*0.5)
    gref_two_prompt.SetPointError(0, default_trail.hv2vsptprompt.GetBinError(iPtBin+1),
                               default_trail.hv2vsptprompt.GetBinError(iPtBin+1),
                               hsyst_prompt.GetMaximum()*0.5, hsyst_prompt.GetMaximum()*0.5)
    SetObjectStyle(gref_two_prompt, markerstyle=20, markercolor=kAzure+2,
                   markersize=1, linecolor=kAzure+2,
                   linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    
    hsyst_prompt.GetYaxis().SetRangeUser(0, hsyst_prompt.GetMaximum()*1.8)
    hsyst_prompt.GetXaxis().SetRangeUser(hsyst_prompt.GetXaxis().GetXmin(), hsyst_prompt.GetXaxis().GetXmax())
    hsyst_final_prompt.SetBinContent(iPtBin, syst_prompt)
    hsyst_final_prompt_rel.SetBinContent(iPtBin, syst_prompt/default_trail.hv2vsptprompt.GetBinContent(iPtBin+1))
    hsyst_prompt.Draw('same')
    gsyst_prompt.Draw('2')
    gref_prompt.Draw('2')
    gref_two_prompt.Draw('2')
    leg_prompt.Draw()

    # Pad 4: non-prompt v2 systematic uncertainty
    canvas.cd(4)
    leg_fd = TLegend(0.6, 0.6, 0.9, 0.9)
    leg_fd.SetBorderSize(0)
    leg_fd.SetFillStyle(0)
    leg_fd.SetTextSize(0.03)
    leg_fd.SetHeader(f'ptmin < #it{{p}}_{{T}} < ptmax')
    leg_fd.AddEntry(hvn_fd, f'non-prompt v2 vs trial ({hvn_fd.GetEntries()} trials)', 'p')
    leg_fd.AddEntry(gsyst_fd, f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst_fd:.3f}', 'f')
    leg_fd.AddEntry(gref_fd, 'non-prompt v2 stat. unc.', 'f')
        
    gref_two_fd = gref_fd.Clone()
    gref_two_fd.SetPoint(0, 0, hsyst_fd.GetMaximum()*0.5)
    gref_two_fd.SetPointError(0, default_trail.hv2vsptfd.GetBinError(iPtBin+1),
                               default_trail.hv2vsptfd.GetBinError(iPtBin+1),
                               hsyst_fd.GetMaximum()*0.5, hsyst_fd.GetMaximum()*0.5)
    SetObjectStyle(gref_two_fd, markerstyle=20, markercolor=kAzure+2,
                   markersize=1, linecolor=kAzure+2,
                   linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        
    hsyst_fd.GetYaxis().SetRangeUser(0, hsyst_fd.GetMaximum()*1.8)
    hsyst_fd.GetXaxis().SetRangeUser(hsyst_fd.GetXaxis().GetXmin(), hsyst_fd.GetXaxis().GetXmax())
    hsyst_final_fd.SetBinContent(iPtBin, syst_fd)
    hsyst_final_fd_rel.SetBinContent(iPtBin, syst_fd/default_trail.hv2vsptfd.GetBinContent(iPtBin+1))
    hsyst_fd.Draw('same')
    gsyst_fd.Draw('2')
    gref_fd.Draw('2')
    gref_two_fd.Draw('2')
    leg_fd.Draw()
    canvas.Update()

    canvsyst_prompt = TCanvas('canvsyst_prompt', 'canvsyst_prompt', 800, 800)
    canvsyst_prompt.cd()
    canvsyst_prompt.SetGrid()
    SetObjectStyle(hsyst_final_prompt, markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hsyst_final_prompt.Draw('same e')
    SetObjectStyle(hsyst_final_prompt_rel, markerstyle=20, markercolor=kAzure+3,
                   markersize=1.,linecolor=kOrange+3,
                   linewidth=2, fillcolor=kOrange+3, fillstyle=3135, fillalpha=0.5)
    
    canvsyst_fd = TCanvas('canvsyst_fd', 'canvsyst_fd', 800, 800)
    canvsyst_fd.cd()
    canvsyst_fd.SetGrid()
    SetObjectStyle(hsyst_final_fd, markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hsyst_final_fd.Draw('same e')
    SetObjectStyle(hsyst_final_fd_rel, markerstyle=20, markercolor=kAzure+3,
                   markersize=1.,linecolor=kOrange+3,
                   linewidth=2, fillcolor=kOrange+3, fillstyle=3135, fillalpha=0.5)

    canvas.SaveAs(f'{outputdir}/Systv2.pdf')
    canvsyst_prompt.SaveAs(f'{outputdir}/SystPromptv2_vs_pt.pdf)')
    canvsyst_fd.SaveAs(f'{outputdir}/SystFDv2_vs_pt.pdf)')

    outdir = os.path.join(outputdir, 'syst_multitrial')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfilename = os.path.join(outdir, 'syst_v2.root')
    outfile = TFile(outfilename, 'recreate')
    hvn_prompt.Write(f'hvn_prompt')
    hvn_fd.Write(f'hvn_fd')
    hsyst_prompt.Write(f'hsyst_prompt')
    hsyst_fd.Write(f'hsyst_fd')
    hsyst_final_prompt.Write()
    hsyst_final_prompt_rel.Write()
    hsyst_final_fd.Write()
    hsyst_final_fd_rel.Write()
    outfile.Close()

def process_trial(default_trail, trails, hvn, hsyst, ipt, prompt=True):        
    for j, trail in enumerate(trails):    # loop over trials
        if prompt:
            hv2vspt = trail.hv2vsptprompt
            default_hv2vspt = default_trail.hv2vsptprompt
            
        else:
            hv2vspt = trail.hv2vsptfd
            default_hv2vspt = default_trail.hv2vsptfd

        if hv2vspt is None:
            print(f'No v2 vs pt found for trial {j}: {trail.trail_path}')
            continue

        skip_trial = False
            
        for iFile in range(len(trail.ryfiles)):
            significance = trail.hsignificances[iFile].GetBinContent(ipt)
            chi2 = trail.hchi2s[iFile].GetBinContent(ipt)

            if (chi2 > 5 and chi2 != 0) or (significance < 6 or significance > 600):
                print(f'Skipping trial {j}: {trail.ryfiles[iFile]} for pt bin {ipt} with chi2 = {chi2} and significance = {significance}')
                skip_trial = True
                break
            
        if skip_trial:
            continue
        
        if hv2vspt.GetBinContent(ipt) == 0:
            print(f'v2 in pt bin {ipt} for trial {j} is 0: {trail.trail_path}')
            continue
        if prompt:
            if hv2vspt.GetBinContent(ipt) < 0:
                print(f'v2 in pt bin {ipt} for trial {j} is negative: {trail.trail_path}')
                print(f'v2 = {hv2vspt.GetBinContent(ipt)}')
                continue
            
        hvn[-1].SetBinContent(j, hv2vspt.GetBinContent(ipt))
        hvn[-1].SetBinError(j, hv2vspt.GetBinError(ipt))
        hsyst[-1].Fill(hv2vspt.GetBinContent(ipt) - default_hv2vspt.GetBinContent(ipt))

    return hvn, hsyst

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('rypathsyst', metavar='text', default='path to the directory containing the .root files from multitrial')
    parser.add_argument('ry_default', metavar='text', default='default .root file')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    parser.add_argument("--prompt", "-p", action="store_true", help="compute systematics for prompt D mesons")
    args = parser.parse_args()

    if args.prompt:
        SystMultitrial.prompt = True
        
    print(f"CIAO MULTITRIAL BDT")

    print(f"args.rypathsyst: {args.rypathsyst}")
    print(f"args.ry_default: {args.ry_default}")
    print(f"args.outputdir: {args.outputdir}")
    # quit()

    compute_syst_multitrial(args.rypathsyst,
                            args.ry_default,
                            args.outputdir)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
# def compute_systematics(outputdir, gvn_vs_mass_default, gvn_vs_mass, hchi2, hsignificance, trails, iFile=0):
#     print("ENTERED IN compute_systematics")
#     hvn, hsyst, hchi2_vs_trial, hsignificance_vs_trial = [], [], [], []
#     gref, gref_two, gsyst = [], [], []
#     hsyst_final = hchi2_trails0].Clone('hsyst_final')
#     hsyst_final.SetTitle('hsyst_final;#it{p}_{T} (GeV/#it{c});syst (%)')
#     hsyst_final.Reset()
#     canvas, leg = [], []

#     for i in range(0, gvn_vs_mass[0].GetN()): # loop over pt bins
#         canvas.append(TCanvas(f'syst_multitrail_pt{gvn_vs_mass[0].GetX()[i]}', f'syst_multitrail_pt{gvn_vs_mass[0].GetX()[i] - gvn_vs_mass[0].GetEXlow()[i]}_{gvn_vs_mass[0].GetX()[i] + gvn_vs_mass[0].GetEXhigh()[i]}', 800, 800))
#         canvas[-1].cd().Divide(2, 2)
#         hvn.append(TH1F('hvn', 'hvn;trial;#it{v}_{n}', len(gvn_vs_mass), 0, len(gvn_vs_mass)+1))
#         hsignificance_vs_trial.append(TH1F('hsignificance_vs_trial', 'hsignificance_vs_trial;trial;significance', len(gvn_vs_mass), 0, len(gvn_vs_mass)+1))
#         hchi2_vs_trial.append(TH1F('hchi2_vs_trial', 'hchi2_vs_trial;trial;#chi^{2}', len(gvn_vs_mass), 0, len(gvn_vs_mass)+1))
#         hsyst.append(TH1F('hsyst', 'hsyst;#it{v}_{n}^{trial} - #it{v}_{n}^{default};"', 200, -0.2, 0.2))
        
#         SetObjectStyle(hvn[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
#         SetObjectStyle(hsignificance_vs_trial[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
#         SetObjectStyle(hchi2_vs_trial[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
#         SetObjectStyle(hsyst[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
#         ipt = i+1
#         for j, _ in enumerate(gvn_vs_mass):    # loop over trials
#             if trails[j].gvnsimfits == None:
#                 print(f'No prompt v2 vs pt found for trial {j}: {trails[j].trail_path}')
#                 continue
#             significance = hsignificance[j].GetBinContent(ipt)
#             chi2 = hchi2[j].GetBinContent(ipt)

#             if len(trails[j].ryfiles) <= iFile:
#                 print(f'No ry file found for trial {j}: {trails[j].trail_path}')
#                 continue
#             # Skip chi2 higher than 5 and significance lower than 3
#             if (chi2 > 4 and chi2 != 0) or (significance < 6 or significance > 600):
#                 print(f'Skipping trial {j}: {trails[j].ryfiles[iFile]} for pt bin {ipt} with chi2 = {chi2} and significance = {significance}')
#                 continue
#             hchi2_vs_trial[-1].SetBinContent(j, chi2)
#             hchi2_vs_trial[-1].SetBinError(j, hchi2[j].GetBinError(ipt))
#             hsignificance_vs_trial[-1].SetBinContent(j, significance)
#             hsignificance_vs_trial[-1].SetBinError(j, hsignificance[j].GetBinError(ipt))
#             hvn[-1].SetBinContent(j, gvn_vs_mass[j].GetY()[i])
#             hvn[-1].SetBinError(j, gvn_vs_mass[j].GetEYlow()[i])
#             hsyst[-1].Fill(gvn_vs_mass[j].GetY()[i] - gvn_vs_mass_default.GetY()[i])
#         #input(f'pt bin {ipt} done. Number of trials skipped = {counter}')

#         # Compute systematic uncertainty
#         rms = hsyst[-1].GetRMS()
#         mean = hsyst[-1].GetMean()
#         syst = np.sqrt(rms**2 + mean**2)
#         print(f'pt bin {ipt}, mean = {mean}, rms = {rms}, syst = {syst}')
#         gsyst.append(TGraphAsymmErrors())
#         gsyst[-1].SetPoint(0, 0, hsyst[-1].GetMaximum()*0.5)
#         gsyst[-1].SetPointError(0, syst, syst,
#                                 hsyst[-1].GetMaximum()*0.5, hsyst[-1].GetMaximum()*0.5)
#         SetObjectStyle(gsyst[-1], markerstyle=20, markercolor=kOrange+2,
#                        markersize=1, linecolor=kOrange+2,
#                           linewidth=2, fillcolor=kOrange+2, fillstyle=3153,
#                           fillalpha=0.5, linestyle=9)

#         # Pad 1: vn vs trial
#         canvas[-1].cd(1).SetGrid()
#         # Define reference line
#         gref.append(TGraphAsymmErrors())
#         gref[-1].SetPoint(0, 0, gvn_vs_mass_default.GetY()[i])
#         gref[-1].SetPointError(0, 0, 0, gvn_vs_mass_default.GetEYlow()[i],
#                            gvn_vs_mass_default.GetEYhigh()[i])
#         gref[-1].SetPoint(1, len(gvn_vs_mass), gvn_vs_mass_default.GetY()[i])
#         gref[-1].SetPointError(1, 0, 0, gvn_vs_mass_default.GetEYlow()[i],
#                            gvn_vs_mass_default.GetEYhigh()[i])
#         SetObjectStyle(gref[-1], markerstyle=20, markercolor=kAzure+2,
#                        markersize=0, linecolor=kAzure+2,
#                        linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
#         hvn[-1].Draw('same')
#         gref[-1].Draw('c3 same')
#         canvas[-1].cd(2)

#         # Legend
#         leg.append(TLegend(0.6, 0.6, 0.9, 0.9))
#         leg[-1].SetBorderSize(0)
#         leg[-1].SetFillStyle(0)
#         leg[-1].SetTextSize(0.03)
#         leg[-1].SetHeader(f'{gvn_vs_mass[0].GetX()[i] - gvn_vs_mass[0].GetEXlow()[i]} < #it{{p}}_{{T}} < {gvn_vs_mass[0].GetX()[i] + gvn_vs_mass[0].GetEXhigh()[i]} GeV/#it{{c}}')
#         leg[-1].AddEntry(hvn[-1], f'vn vs trial ({hvn[-1].GetEntries()} trials)', 'p')
#         leg[-1].AddEntry(gsyst[-1], f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst:.3f}', 'f')
#         leg[-1].AddEntry(gref[-1], 'vn stat. unc.', 'f')

#         # Define reference vertical line at 1
#         gref_two.append(gref[-1].Clone())
#         gref_two[-1].SetPoint(0, 0, hsyst[-1].GetMaximum()*0.5)
#         gref_two[-1].SetPointError(0, gvn_vs_mass_default.GetEYlow()[i],
#                                    gvn_vs_mass_default.GetEYhigh()[i],
#                                    hsyst[-1].GetMaximum()*0.5, hsyst[-1].GetMaximum()*0.5)
#         SetObjectStyle(gref_two[-1], markerstyle=20, markercolor=kAzure+2,
#                        markersize=1, linecolor=kAzure+2,
#                        linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
#         hsyst[-1].GetYaxis().SetRangeUser(0, hsyst[-1].GetMaximum()*1.8)
#         hsyst[-1].GetXaxis().SetRangeUser(hsyst[-1].GetXaxis().GetXmin(), hsyst[-1].GetXaxis().GetXmax())
#         hsyst[-1].Draw('same')
#         gsyst[-1].Draw('2')
#         gref_two[-1].Draw('2')
#         gsyst[-1].Draw('2')
#         leg[-1].Draw()
#         # Pad 3: chi2 vs trial
#         canvas[-1].cd(3)
#         hchi2_vs_trial[-1].Draw('same')
#         # Pad 4: significance vs trial
#         canvas[-1].cd(4)
#         hsignificance_vs_trial[-1].Draw('same')
#         hsyst_final.SetBinContent(ipt, syst)
#         canvas[-1].Update()

#     canvsyst = TCanvas('canvsyst', 'canvsyst', 800, 800)
#     canvsyst.cd()
#     canvsyst.SetGrid()
#     SetObjectStyle(hsyst_final, markerstyle=20, markercolor=kAzure+2,
#                    markersize=1.,linecolor=kOrange+2,
#                    linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
#     hsyst_final.Draw('same e')

#     for icanv, canv in enumerate(canvas):
#         if icanv == 0:
#             suffix_pdf = '('
#         elif icanv == len(canvas) - 1:
#             suffix_pdf = ')'
#         else:
#             suffix_pdf = ''
#         if len(canvas) == 1:
#             suffix_pdf = ''
#         canv.SaveAs(f'{outputdir}/SystRy_{iFile}.pdf{suffix_pdf}')
#     canvsyst.SaveAs(f'{outputdir}/SystRy_vs_pt_{iFile}.pdf)')

#     outdir = os.path.join(outputdir, 'syst_multitrial')
#     if not os.path.exists(outdir):
#         os.makedirs(outdir)
#     outfilename = os.path.join(outdir, f'syst_ry_{iFile}.root')
#     outfile = TFile(outfilename, 'recreate')
#     for h in hvn:
#         h.Write()
#     for h in hsyst:
#         h.Write()
#     hsyst_final.Write()
#     outfile.Close()

# def compute_systematics_prompt(outputdir, default_trail, trails):

#     hvn_prompt, hsyst_prompt, hvn_fd, hsyst_fd = [], [], [], []
#     gref_prompt, gref_two_prompt, gsyst_prompt = [], [], []
#     gref_fd, gref_two_fd, gsyst_fd = [], [], []
#     hsyst_final_prompt = trails[0].hchi2s[0].Clone('hsyst_final_prompt')
#     hsyst_final_prompt_rel = trails[0].hchi2s[0].Clone('hsyst_final_prompt_rel')
#     hsyst_final_fd = trails[0].hchi2s[0].Clone('hsyst_final_fd')
#     hsyst_final_fd_rel = trails[0].hchi2s[0].Clone('hsyst_final_fd_rel')
#     hsyst_final_prompt.SetTitle('hsyst_final_prompt;#it{p}_{T} (GeV/#it{c});syst (abs)')
#     hsyst_final_prompt_rel.SetTitle('hsyst_final_prompt_rel;#it{p}_{T} (GeV/#it{c});syst (%)')
#     hsyst_final_fd.SetTitle('hsyst_final_fd;#it{p}_{T} (GeV/#it{c});syst (abs)')
#     hsyst_final_fd_rel.SetTitle('hsyst_final_fd_rel;#it{p}_{T} (GeV/#it{c});syst (%)')
#     hsyst_final_prompt.Reset()
#     hsyst_final_prompt_rel.Reset()
#     hsyst_final_fd.Reset()
#     hsyst_final_fd_rel.Reset()
#     canvas, leg_prompt, leg_fd = [], [], []
#     for i in range(0, default_trail.gvnsimfits[0].GetN()): # loop over pt bins
#         canvas.append(TCanvas(f'syst_multitrail_pt{default_trail.gvnsimfits[0].GetX()[i]}', f'syst_multitrail_pt{default_trail.gvnsimfits[0].GetX()[i] - default_trail.gvnsimfits[0].GetEXlow()[i]}_{default_trail.gvnsimfits[0].GetX()[i] + default_trail.gvnsimfits[0].GetEXhigh()[i]}', 800, 800))
#         canvas[-1].cd().Divide(2, 2)
#         hvn_prompt.append(TH1F('hvn_prompt', 'hvn_prompt;trial;prompt #it{v}_{2}', len(trails), 0, len(trails)+1))
#         hvn_fd.append(TH1F('hvn_fd', 'hvn_fd;trial;non-prompt #it{v}_{2}', len(trails), 0, len(trails)+1))
#         hsyst_prompt.append(TH1F('hsyst_prompt', 'hsyst_prompt;prompt #it{v}_{2}^{trial} - prompt #it{v}_{2}^{default};"', 200, -0.2, 0.2))
#         hsyst_fd.append(TH1F('hsyst_fd', 'hsyst_fd;non-prompt #it{v}_{2}^{trial} - non-prompt #it{v}_{2}^{default};"', 200, -0.2, 0.2))
        
#         SetObjectStyle(hvn_prompt[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
#         SetObjectStyle(hvn_fd[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
#         SetObjectStyle(hsyst_prompt[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
#         SetObjectStyle(hsyst_fd[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)

#         ipt = i+1
#         hvn_prompt, hsyst_prompt = process_trial(default_trail, trails, hvn_prompt, hsyst_prompt, ipt, prompt=True)
#         hvn_fd, hsyst_fd = process_trial(default_trail, trails, hvn_fd, hsyst_fd, ipt, prompt=False)
#         #input(f'pt bin {ipt} done. Number of trials skipped = {counter}')

#         # Compute systematic uncertainty
#         rms_prompt = hsyst_prompt[-1].GetRMS()
#         mean_prompt = hsyst_prompt[-1].GetMean()
#         syst_prompt = np.sqrt(rms_prompt**2 + mean_prompt**2)
#         print(f'pt bin {ipt}, mean = {mean_prompt}, rms = {rms_prompt}, syst = {syst_prompt}')
#         gsyst_prompt.append(TGraphAsymmErrors())
#         gsyst_prompt[-1].SetPoint(0, 0, hsyst_prompt[-1].GetMaximum()*0.5)
#         gsyst_prompt[-1].SetPointError(0, syst_prompt, syst_prompt,
#                                 hsyst_prompt[-1].GetMaximum()*0.5, hsyst_prompt[-1].GetMaximum()*0.5)
#         SetObjectStyle(gsyst_prompt[-1], markerstyle=20, markercolor=kOrange+2, markersize=1, linecolor=kOrange+2,
#                        linewidth=2, fillcolor=kOrange+2, fillstyle=3153, fillalpha=0.5, linestyle=9)

#         rms_fd = hsyst_fd[-1].GetRMS()
#         mean_fd = hsyst_fd[-1].GetMean()
#         syst_fd = np.sqrt(rms_fd**2 + mean_fd**2)
#         print(f'pt bin {ipt}, mean = {mean_fd}, rms = {rms_fd}, syst = {syst_fd}')
#         gsyst_fd.append(TGraphAsymmErrors())
#         gsyst_fd[-1].SetPoint(0, 0, hsyst_fd[-1].GetMaximum()*0.5)
#         gsyst_fd[-1].SetPointError(0, syst_fd, syst_fd,
#                                 hsyst_fd[-1].GetMaximum()*0.5, hsyst_fd[-1].GetMaximum()*0.5)
#         SetObjectStyle(gsyst_fd[-1], markerstyle=20, markercolor=kOrange+2, markersize=1, linecolor=kOrange+2,
#                        linewidth=2, fillcolor=kOrange+2, fillstyle=3153, fillalpha=0.5, linestyle=9)

#         # Pad 1: prompt v2 vs trial
#         canvas[-1].cd(1).SetGrid()
#         # Define reference line
#         gref_prompt.append(TGraphAsymmErrors())
#         gref_prompt[-1].SetPoint(0, 0, default_trail.hv2vsptprompt.GetBinContent(ipt))
#         gref_prompt[-1].SetPointError(0, 0, 0, default_trail.hv2vsptprompt.GetBinError(ipt),
#                            default_trail.hv2vsptprompt.GetBinError(ipt))
#         gref_prompt[-1].SetPoint(1, len(trails), default_trail.hv2vsptprompt.GetBinContent(ipt))
#         gref_prompt[-1].SetPointError(1, 0, 0, default_trail.hv2vsptprompt.GetBinError(ipt),
#                            default_trail.hv2vsptprompt.GetBinError(ipt))
#         SetObjectStyle(gref_prompt[-1], markerstyle=20, markercolor=kAzure+2,
#                        markersize=0, linecolor=kAzure+2,
#                        linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
#         hvn_prompt[-1].Draw('same')
#         gref_prompt[-1].Draw('c3 same')
        
#         # Pad 2: non-prompt v2 vs trial
#         canvas[-1].cd(2)
#         # Define reference line
#         gref_fd.append(TGraphAsymmErrors())
#         gref_fd[-1].SetPoint(0, 0, default_trail.hv2vsptfd.GetBinContent(ipt))
#         gref_fd[-1].SetPointError(0, 0, 0, default_trail.hv2vsptfd.GetBinError(ipt),
#                            default_trail.hv2vsptfd.GetBinError(ipt))
#         gref_fd[-1].SetPoint(1, len(trails), default_trail.hv2vsptfd.GetBinContent(ipt))
#         gref_fd[-1].SetPointError(1, 0, 0, default_trail.hv2vsptfd.GetBinError(ipt),
#                            default_trail.hv2vsptfd.GetBinError(ipt))
#         SetObjectStyle(gref_fd[-1], markerstyle=20, markercolor=kAzure+2,
#                        markersize=0, linecolor=kAzure+2,
#                        linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
#         hvn_fd[-1].Draw('same')
#         gref_fd[-1].Draw('c3 same')
        
#         # Pad 3: prompt v2 systematic uncertainty
#         canvas[-1].cd(3)
#         leg_prompt.append(TLegend(0.6, 0.6, 0.9, 0.9))
#         leg_prompt[-1].SetBorderSize(0)
#         leg_prompt[-1].SetFillStyle(0)
#         leg_prompt[-1].SetTextSize(0.03)
#         leg_prompt[-1].SetHeader(f'{default_trail.gvnsimfits[0].GetX()[i] - default_trail.gvnsimfits[0].GetEXlow()[i]} < #it{{p}}_{{T}} < {default_trail.gvnsimfits[0].GetX()[i] + default_trail.gvnsimfits[0].GetEXhigh()[i]} GeV/#it{{c}}')
#         leg_prompt[-1].AddEntry(hvn_prompt[-1], f'prompt v2 vs trial ({hvn_prompt[-1].GetEntries()} trials)', 'p')
#         leg_prompt[-1].AddEntry(gsyst_prompt[-1], f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst_prompt:.3f}', 'f')
#         leg_prompt[-1].AddEntry(gref_prompt[-1], 'prompt v2 stat. unc.', 'f')
        
#         gref_two_prompt.append(gref_prompt[-1].Clone())
#         gref_two_prompt[-1].SetPoint(0, 0, hsyst_prompt[-1].GetMaximum()*0.5)
#         gref_two_prompt[-1].SetPointError(0, default_trail.hv2vsptprompt.GetBinError(ipt),
#                                    default_trail.hv2vsptprompt.GetBinError(ipt),
#                                    hsyst_prompt[-1].GetMaximum()*0.5, hsyst_prompt[-1].GetMaximum()*0.5)
#         SetObjectStyle(gref_two_prompt[-1], markerstyle=20, markercolor=kAzure+2,
#                        markersize=1, linecolor=kAzure+2,
#                        linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        
#         hsyst_prompt[-1].GetYaxis().SetRangeUser(0, hsyst_prompt[-1].GetMaximum()*1.8)
#         hsyst_prompt[-1].GetXaxis().SetRangeUser(hsyst_prompt[-1].GetXaxis().GetXmin(), hsyst_prompt[-1].GetXaxis().GetXmax())
#         hsyst_final_prompt.SetBinContent(ipt, syst_prompt)
#         hsyst_final_prompt_rel.SetBinContent(ipt, syst_prompt/default_trail.hv2vsptprompt.GetBinContent(ipt))
#         hsyst_prompt[-1].Draw('same')
#         gsyst_prompt[-1].Draw('2')
#         gref_prompt[-1].Draw('2')
#         gref_two_prompt[-1].Draw('2')
#         leg_prompt[-1].Draw()

#         # Pad 4: non-prompt v2 systematic uncertainty
#         canvas[-1].cd(4)
#         leg_fd.append(TLegend(0.6, 0.6, 0.9, 0.9))
#         leg_fd[-1].SetBorderSize(0)
#         leg_fd[-1].SetFillStyle(0)
#         leg_fd[-1].SetTextSize(0.03)
#         leg_fd[-1].SetHeader(f'{default_trail.gvnsimfits[0].GetX()[i] - default_trail.gvnsimfits[0].GetEXlow()[i]} < #it{{p}}_{{T}} < {default_trail.gvnsimfits[0].GetX()[i] + default_trail.gvnsimfits[0].GetEXhigh()[i]} GeV/#it{{c}}')
#         leg_fd[-1].AddEntry(hvn_fd[-1], f'non-prompt v2 vs trial ({hvn_fd[-1].GetEntries()} trials)', 'p')
#         leg_fd[-1].AddEntry(gsyst_fd[-1], f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst_fd:.3f}', 'f')
#         leg_fd[-1].AddEntry(gref_fd[-1], 'non-prompt v2 stat. unc.', 'f')
        
#         gref_two_fd.append(gref_fd[-1].Clone())
#         gref_two_fd[-1].SetPoint(0, 0, hsyst_fd[-1].GetMaximum()*0.5)
#         gref_two_fd[-1].SetPointError(0, default_trail.hv2vsptfd.GetBinError(ipt),
#                                    default_trail.hv2vsptfd.GetBinError(ipt),
#                                    hsyst_fd[-1].GetMaximum()*0.5, hsyst_fd[-1].GetMaximum()*0.5)
#         SetObjectStyle(gref_two_fd[-1], markerstyle=20, markercolor=kAzure+2,
#                        markersize=1, linecolor=kAzure+2,
#                        linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        
#         hsyst_fd[-1].GetYaxis().SetRangeUser(0, hsyst_fd[-1].GetMaximum()*1.8)
#         hsyst_fd[-1].GetXaxis().SetRangeUser(hsyst_fd[-1].GetXaxis().GetXmin(), hsyst_fd[-1].GetXaxis().GetXmax())
#         hsyst_final_fd.SetBinContent(ipt, syst_fd)
#         hsyst_final_fd_rel.SetBinContent(ipt, syst_fd/default_trail.hv2vsptfd.GetBinContent(ipt))
#         hsyst_fd[-1].Draw('same')
#         gsyst_fd[-1].Draw('2')
#         gref_fd[-1].Draw('2')
#         gref_two_fd[-1].Draw('2')
#         leg_fd[-1].Draw()
#         canvas[-1].Update()

#     canvsyst_prompt = TCanvas('canvsyst_prompt', 'canvsyst_prompt', 800, 800)
#     canvsyst_prompt.cd()
#     canvsyst_prompt.SetGrid()
#     SetObjectStyle(hsyst_final_prompt, markerstyle=20, markercolor=kAzure+2,
#                    markersize=1.,linecolor=kOrange+2,
#                    linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
#     hsyst_final_prompt.Draw('same e')
#     SetObjectStyle(hsyst_final_prompt_rel, markerstyle=20, markercolor=kAzure+3,
#                    markersize=1.,linecolor=kOrange+3,
#                    linewidth=2, fillcolor=kOrange+3, fillstyle=3135, fillalpha=0.5)
#     # hsyst_final_prompt_rel.Draw('same e')
    
#     canvsyst_fd = TCanvas('canvsyst_fd', 'canvsyst_fd', 800, 800)
#     canvsyst_fd.cd()
#     canvsyst_fd.SetGrid()
#     SetObjectStyle(hsyst_final_fd, markerstyle=20, markercolor=kAzure+2,
#                    markersize=1.,linecolor=kOrange+2,
#                    linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
#     hsyst_final_fd.Draw('same e')
#     SetObjectStyle(hsyst_final_fd_rel, markerstyle=20, markercolor=kAzure+3,
#                    markersize=1.,linecolor=kOrange+3,
#                    linewidth=2, fillcolor=kOrange+3, fillstyle=3135, fillalpha=0.5)
#     # hsyst_final_fd_rel.Draw('same e')

#     for icanv, canv in enumerate(canvas):
#         if icanv == 0:
#             suffix_pdf = '('
#         elif icanv == len(canvas) - 1:
#             suffix_pdf = ')'
#         else:
#             suffix_pdf = ''
#         if len(canvas) == 1:
#             suffix_pdf = ''
#         canv.SaveAs(f'{outputdir}/Systv2.pdf{suffix_pdf}')
#     canvsyst_prompt.SaveAs(f'{outputdir}/SystPromptv2_vs_pt.pdf)')
#     canvsyst_fd.SaveAs(f'{outputdir}/SystFDv2_vs_pt.pdf)')

#     outdir = os.path.join(outputdir, 'syst_multitrial')
#     if not os.path.exists(outdir):
#         os.makedirs(outdir)
#     outfilename = os.path.join(outdir, 'syst_v2.root')
#     outfile = TFile(outfilename, 'recreate')
#     for iHist, (h_prompt, h_fd) in enumerate(zip(hvn_prompt, hvn_fd)):
#         h_prompt.Write(f'hvn_prompt_{iHist}')
#         h_fd.Write(f'hvn_fd_{iHist}')
#     for iHist, (h_prompt, h_fd) in enumerate(zip(hsyst_prompt, hsyst_fd)):
#         h_prompt.Write(f'hsyst_prompt_{iHist}')
#         h_fd.Write(f'hsyst_fd_{iHist}')
#     hsyst_final_prompt.Write()
#     hsyst_final_prompt_rel.Write()
#     hsyst_final_fd.Write()
#     hsyst_final_fd_rel.Write()
#     outfile.Close()