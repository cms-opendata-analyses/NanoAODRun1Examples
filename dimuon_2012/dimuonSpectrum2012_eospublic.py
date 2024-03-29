# this example is a modified version of the one on 
# http://opendata.cern.ch/record/12342

import ROOT

# Enable multi-threading
# The default here is set to a single thread. You can choose the number of threads based on your system.
ROOT.ROOT.EnableImplicitMT()

# Create dataframe from NanoAODRun1 files  
df = ROOT.RDataFrame("Events", {"root://eospublic.cern.ch//eos/opendata/cms/derived-data/NanoAODRun1/01-Jul-22/Run2012B_DoubleMuParked_merged.root", "root://eospublic.cern.ch//eos/opendata/cms/derived-data/NanoAODRun1/01-Jul-22/Run2012C_DoubleMuParked_merged.root"})

# Select events with at least two muons
df_2mu = df.Filter("nMuon > 1", "Events with at least two muons")

# Select events with two muons of opposite charge and pT> 3 GeV 
df_os = df_2mu.Filter("(Muon_charge[0] != Muon_charge[1]) && (Muon_pt[0] > 3 && Muon_pt[1] > 3)", "Muons with opposite charge")

# Compute invariant mass of the dimuon system
# The following code just-in-time compiles the C++ function to compute
# the invariant mass, so that the function can be called in the Define node of
# the ROOT dataframe.
ROOT.gInterpreter.Declare(
"""
using namespace ROOT::VecOps;
float computeInvariantMass(RVec<float>& pt, RVec<float>& eta, RVec<float>& phi, RVec<float>& mass) {
    ROOT::Math::PtEtaPhiMVector m1(pt[0], eta[0], phi[0], mass[0]);
    ROOT::Math::PtEtaPhiMVector m2(pt[1], eta[1], phi[1], mass[1]);
    return (m1 + m2).mass();
}
""")
df_mass = df_os.Define("Dimuon_mass", "computeInvariantMass(Muon_pt, Muon_eta, Muon_phi, Muon_mass)")

# Book histogram of dimuon mass spectrum
bins = 30000 # Number of bins in the histogram
low = 0.25 # Lower edge of the histogram
up = 300.0 # Upper edge of the histogram
hist = df_mass.Histo1D(ROOT.RDF.TH1DModel("", "", bins, low, up), "Dimuon_mass")

# Request cut-flow report
report = df_mass.Report()

# Create canvas for plotting
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetTextFont(42)
c = ROOT.TCanvas("c", "", 800, 700)
c.SetLogx()
c.SetLogy()

# Draw histogram
hist.GetXaxis().SetTitle("m_{#mu#mu} (GeV)")
hist.GetXaxis().SetTitleSize(0.04)
hist.GetYaxis().SetTitle("N_{Events}")
hist.GetYaxis().SetTitleSize(0.04)
hist.Draw()

# Draw labels
label = ROOT.TLatex()
label.SetTextAlign(22)
label.DrawLatex(0.55, 4.7e4, "#eta")
label.DrawLatex(0.77, 8.0e4, "#rho,#omega")
label.DrawLatex(1.20, 6.0e4, "#phi")
label.DrawLatex(4.40, 4.0e5, "J/#psi")
label.DrawLatex(4.60, 2.0e4, "#psi'")
label.DrawLatex(12.0, 2.7e4, "Y(1,2,3S)")
label.DrawLatex(91.0, 1.8e4, "Z")
label.SetNDC(True)
label.SetTextAlign(11)
label.SetTextSize(0.04)
label.DrawLatex(0.10, 0.92, "#bf{CMS Open Data, NanoAODRun1}")
label.DrawLatex(0.10, 0.96, "#bf{2012 Run B + Run C with muon events >= 2 and pT-cut}")
label.SetTextAlign(31)
label.DrawLatex(0.95, 0.92, "#sqrt{s} = 8 TeV, L_{int} = 11.6 fb^{-1}")

# Save plot
c.SaveAs("dimuonSpectrum2012_py_eospublic.pdf")

# Print cut-flow report
report.Print()
