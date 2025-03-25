import ctypes
from ROOT import TFile, TH1, TF1, TCanvas, gStyle, TLegend, kRed, kBlue, kGreen, kMagenta, kCyan, kOrange
import argparse

parser = argparse.ArgumentParser(description="Fit psi2 distribution and save the plot.")
parser.add_argument("outdir", type=str, help="Output directory to save the plot")
parser.add_argument("cent", type=str, help="cent string")
parser.add_argument("input_trig_file", type=str, help="Path to file with hSparseEp sparse")
args = parser.parse_args()
outdir = args.outdir

# Harmonic number
harmonic = 2  

# File and histogram paths
sparse_trig_path = "hf-task-flow-charm-hadrons/ep/hSparseEp"
infile = TFile.Open(args.input_trig_file, 'READ')

# Check if file is open
if not infile or infile.IsZombie():
    print(f"Error: Cannot open file {args.input_trig_file}")
    exit()

# Retrieve sparse histogram
sparse_trig = infile.Get(sparse_trig_path)

if not sparse_trig:
    print(f"Error: Cannot retrieve histogram {sparse_trig_path}")
    exit()

print(f"Total entries in sparse_trig: {sparse_trig.GetEntries()}")

# Project psi2 angle
psi2_axis = 1
hist_psi2 = sparse_trig.Projection(psi2_axis)
outfile = TFile.Open(f"{outdir}/SystMiscal_{args.cent}.root", "RECREATE")
outfile.cd()
hist_psi2.Write()
hist_psi2.Rebin(4)
hist_psi2.SetStats(0)

# Fit functions for cos and sin components
func_cos_psi2 = TF1("fCosPsi2", "[0] * (1 + 2*[1]*TMath::Cos(2*x))", -3.14/2, 3.14/2)
func_sin_psi2 = TF1("fSinPsi2", "[0] * (1 + 2*[1]*TMath::Sin(2*x))", -3.14/2, 3.14/2)

# Perform fits
hist_psi2.Fit(func_cos_psi2, "SMR+")
avg_cos_2psi = func_cos_psi2.GetParameter(1)

hist_psi2.Fit(func_sin_psi2, "SMR+")
avg_sin_2psi = func_sin_psi2.GetParameter(1)

print(f"⟨cos(2ψ₂)⟩ = {avg_cos_2psi}")
print(f"⟨sin(2ψ₂)⟩ = {avg_sin_2psi}")

# Draw results in TCanvas
cFit = TCanvas("cFit", "Psi2 Fits", 600, 600)
cFit.cd()
cFit.SetLeftMargin(0.15)
gStyle.SetOptFit(1)  # Display fit statistics

hist_psi2.GetXaxis().SetRangeUser(-3.14/2, 3.14/2)
print(f"hist_psi2.GetMinimum(): {hist_psi2.GetMinimum()}")
print(f"hist_psi2.GetMinimum()*0.99999: {hist_psi2.GetMinimum()*0.99999}")
hist_psi2.GetYaxis().SetRangeUser(hist_psi2.GetMaximum()*0.9, hist_psi2.GetMaximum()*1.1)
hist_psi2.SetTitle(f"Fit to #psi_{{2}} Distribution -- cent{args.cent}")
hist_psi2.GetXaxis().SetTitle("#psi_{2} (rad)")
hist_psi2.GetYaxis().SetTitle("Counts")

hist_psi2.SetLineColor(1)  # Red
hist_psi2.Draw('pe')
func_cos_psi2.SetLineColor(2)  # Red
func_sin_psi2.SetLineColor(4)  # Blue
func_cos_psi2.Draw("same")
func_sin_psi2.Draw("same")

# Create legend
legend = TLegend(0.25, 0.65, 0.6, 0.85)  # Adjust position (x1, y1, x2, y2)
legend.SetBorderSize(0)  # No border
legend.SetFillStyle(0)  # Transparent background
legend.SetTextSize(0.04)  # Adjust text size
legend.AddEntry(func_cos_psi2, f"#LTcos(2#psi_{{2}})#GT = {func_cos_psi2.GetParameter(1):.5f} +/- {func_cos_psi2.GetParError(1):.5f}", "l")
legend.AddEntry(func_sin_psi2, f"#LTsin(2#psi_{{2}})#GT = {func_sin_psi2.GetParameter(1):.5f} +/- {func_sin_psi2.GetParError(1):.5f}", "l")
legend.Draw()

outfile.cd()
cFit.Write()

# Get the avg sin and cos of phi2
file_cos = TFile.Open(f"/home/mdicosta/FlowDplus/FinalResults/CheckTrigDistros/cutvar_{args.cent}_cos/ry/raw_yields_{args.cent}_cos_00.root", "r")
gAvgCosPhi2 = file_cos.Get("gvnSimFit")
file_sin = TFile.Open(f"/home/mdicosta/FlowDplus/FinalResults/CheckTrigDistros/cutvar_{args.cent}_sin/ry/raw_yields_{args.cent}_sin_00.root", "r")
gAvgSinPhi2 = file_sin.Get("gvnSimFit")

outfile.cd()
cos_x_values = []
cos_y_values = []
n_points = gAvgCosPhi2.GetN()
for i in range(n_points):
    x, y = ctypes.c_double(0), ctypes.c_double(0)
    gAvgCosPhi2.GetPoint(i, x, y)
    cos_x_values.append(float(x.value))
    cos_y_values.append(float(y.value))

print("[Cos] X values:", cos_x_values)
print("[Cos] Y values:", cos_y_values)

sin_x_values = []
sin_y_values = []
n_points = gAvgSinPhi2.GetN()
for i in range(n_points):
    x, y = ctypes.c_double(0), ctypes.c_double(0)
    gAvgSinPhi2.GetPoint(i, x, y)
    sin_x_values.append(float(x.value))
    sin_y_values.append(float(y.value))
    
print("[Sin] X values:", sin_x_values)
print("[Sin] Y values:", sin_y_values)

avg_cos_psi = func_cos_psi2.GetParameter(1)
avg_sin_psi = func_sin_psi2.GetParameter(1)

cos_delta_phi_psi = [cos_phi*avg_cos_psi + sin_phi*avg_sin_psi for cos_phi, sin_phi in zip(cos_y_values, sin_y_values)]
gAvgDeltaPhiPsi = gAvgCosPhi2.Clone("gAvgDeltaPhiPsi")
gAvgDeltaPhiPsi.SetName("gAvgDeltaPhiPsi")
gAvgDeltaPhiPsi.SetTitle("Average Delta Phi Psi")
for i in range(n_points):
    gAvgDeltaPhiPsi.SetPoint(i, cos_x_values[i], cos_delta_phi_psi[i])
    gAvgDeltaPhiPsi.SetPointError(i, gAvgCosPhi2.GetErrorX(i), gAvgCosPhi2.GetErrorX(i), 0, 0)  # Set Y errors to 0

inverted_cos_delta_phi_psi = [cos_phi*avg_sin_psi + sin_phi*avg_cos_psi for cos_phi, sin_phi in zip(cos_y_values, sin_y_values)]
gSinCosAvgDeltaPhiPsi = gAvgCosPhi2.Clone("gSinCosAvgDeltaPhiPsi")
gSinCosAvgDeltaPhiPsi.SetName("gSinCosAvgDeltaPhiPsi")
gSinCosAvgDeltaPhiPsi.SetTitle("Average Delta Phi Psi")
for i in range(n_points):
    gSinCosAvgDeltaPhiPsi.SetPoint(i, cos_x_values[i], inverted_cos_delta_phi_psi[i])
    gSinCosAvgDeltaPhiPsi.SetPointError(i, gAvgCosPhi2.GetErrorX(i), gAvgCosPhi2.GetErrorX(i), 0, 0)  # Set Y errors to 0

# Create canvas to draw
gAvgDeltaPhiPsi.Write()
gSinCosAvgDeltaPhiPsi.Write()
outfile.Close()

gAvgDeltaPhiPsi.GetYaxis().SetRangeUser(-0.00005, 0.0003)
gAvgDeltaPhiPsi.SetMarkerStyle(20)
gAvgDeltaPhiPsi.SetMarkerColor(2)
cFit.SaveAs(f"{outdir}/SinCosFits_{args.cent}.pdf")
canvas = TCanvas("canvas", "Updated Graph", 600, 600)
gAvgDeltaPhiPsi.Draw("ALP")
canvas.SaveAs(f"{outdir}/FinalSyst_{args.cent}.pdf")

print(f"cos_delta_phi_psi: {cos_delta_phi_psi}")
