#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"

// this example is a modified version of the one on 
// http://opendata.cern.ch/record/12342

using namespace ROOT::VecOps;

// Compute the invariant mass of two muon four-vectors
float computeInvariantMass(RVec<float>& pt, RVec<float>& eta, RVec<float>& phi, RVec<float>& mass) {
    ROOT::Math::PtEtaPhiMVector m1(pt[0], eta[0], phi[0], mass[0]);
    ROOT::Math::PtEtaPhiMVector m2(pt[1], eta[1], phi[1], mass[1]);
    return (m1 + m2).mass();
}


void dimuonSpectrum2012_publicchain() {
    // Enable multi-threading
    // The default here is set to a single thread. You can choose the number of threads based on your system.
    ROOT::EnableImplicitMT();

    // Create dataframe from chain of NanoAODRun1 files on eospublic
    ROOT::RDataFrame df("Events", {"root://eospublic.cern.ch//eos/opendata/cms/derived-data/NanoAODRun1/01-Jul-22/Run2012B_DoubleMuParked/*.root", "root://eospublic.cern.ch//eos/opendata/cms/derived-data/NanoAODRun1/01-Jul-22/Run2012C_DoubleMuParked/*.root"});
    // RDataFrame interfaces to TTree and TChain. The "Events" part makes sure that within the root file, the data frame is taken from within the "Events" folder. 

    // Select events with at least two muons
    auto df_2mu = df.Filter("nMuon >= 2", "Events with two or more muons");
    // This line filters the "nMuon" branch of the "Events" tree to only select events with two or more muons"    

    // Select events with two muons of opposite charge and pt>3 GeV
    auto df_os = df_2mu.Filter("Muon_charge[0] != Muon_charge[1] && (Muon_pt[0] > 3 && Muon_pt[1] > 3)", "Muons with opposite charge");
    // This line filters the "nMuon" branch of the "Events" tree further to only select events with two muons of opposite charge and above pt threshold

    // Compute invariant mass of the dimuon system
    auto df_mass = df_os.Define("Dimuon_mass", computeInvariantMass,
                                {"Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"});
    // This line creates a new column "Dimuon_mass" with values from the computeInvariantMass float function. This takes pt as "Muon_pt", eta as "Muon_eta", phi as "Muon_phi" and mass as "Muon_mass".

    // Book histogram of dimuon mass spectrum
    const auto bins = 30000; // Number of bins in the histogram
    const auto low = 0.25; // Lower edge of the histogram
    const auto up = 300.0; // Upper edge of the histogram
    auto hist = df_mass.Histo1D({"", "", bins, low, up}, "Dimuon_mass");

    // Request cut-flow report
    auto report = df_mass.Report();
    // Obtains statistics on how many entries have been accepted and rejected by the filters. The method returns a ROOT::RDF::RCutFlowReport instance which can be queried programmatically to get information about the effects of the individual cuts. 

    // Create canvas for plotting
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    auto c = new TCanvas("c", "", 800, 700);
    c->SetLogx(); 
    c->SetLogy(); //logs both x and y

    // Draw histogram
    hist->GetXaxis()->SetTitle("m_{#mu#mu} (GeV)");
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitle("N_{Events}");
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->DrawClone();

    // Draw labels
    TLatex label;
    label.SetTextAlign(22);
    label.DrawLatex(0.55, 4.7e4, "#eta");
    label.DrawLatex(0.77, 8.0e4, "#rho,#omega");
    label.DrawLatex(1.20, 6.0e4, "#phi");
    label.DrawLatex(4.40, 4.0e5, "J/#psi");
    label.DrawLatex(4.60, 2.0e4, "#psi'");
    label.DrawLatex(12.0, 2.7e4, "Y(1,2,3S)");
    label.DrawLatex(91.0, 1.8e4, "Z");
    label.SetNDC(true);
    label.SetTextAlign(11);
    label.SetTextSize(0.04);
    label.DrawLatex(0.10, 0.92, "#bf{CMS Open Data, NanoAODrun1}");
    label.DrawLatex(0.10, 0.96, "#bf{2012 Run B + Run C with >= 2 muons and pT-cut}");
    label.SetTextAlign(31);
    label.DrawLatex(0.95, 0.92, "#sqrt{s} = 8 TeV, L_{int} = 11.6 fb^{-1}");

    // Save plot
    c->SaveAs("dimuonSpectrum2012_C_publicchain.pdf");

    // Print cut-flow report
    report->Print();
}


int main() {
    dimuonSpectrum2012_publicchain();
}

