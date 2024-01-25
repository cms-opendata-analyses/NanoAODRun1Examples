#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "Math/Vector4Dfwd.h"
#include "Math/Vector4D.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TStyle.h"
#include <iostream>
#include <chrono>
#include <unistd.h>


#define nThreads 12
#define nBins 620
#define x1 -0.4
#define x2 2.7

using namespace ROOT::VecOps;

void Dimuon2011_eospublic_RDF2() {

  // modified by Anuranjan Sarkar from Dimuon2011_eospublic_RDF, to approximately reproduce arXiv:1609.02366 Figure 68

    auto start = chrono::steady_clock::now();
    // Enable multi-threading
    // The default here is set to a single thread. You can choose the number of threads based on your system.
    ROOT::EnableImplicitMT(nThreads);

    // Run over double muon sample, for high pT double muon
    //ROOT::RDataFrame df_DoubleMu ("Events", "/nfs/dust/cms/user/geiser/eosdata/Run2011A_DoubleMu_merged.root");
    ROOT::RDataFrame df_DoubleMu ("Events", "root://eospublic.cern.ch//eos/opendata/cms/derived-data/NanoAODRun1/01-Jul-22/Run2011A_DoubleMu_merged.root");
    auto filter1 = df_DoubleMu.Filter("run < 170000", "Run number")
        .Filter("Trig_DoubleMuThresh > 12", "Dimuon threshold")
        .Define("Dimu_mass_cut", // name
            []( // lambda expression
            uint nDimu,
            const RVec<int>& Dimu_charge,
            const RVec<float>& Dimu_mass,
            const RVec<int>& Dimu_t1muIdx,
            const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt,
            const RVec<bool>& Muon_mediumId
            ) -> RVec<float>
            { // returned type is RVec<float>
                RVec<float> Dimu_mass_cut; // lambda function's body 
                for (size_t i = 0; i < nDimu; i++) {

                    if  (
                            (Dimu_charge[i] == 0) &&
                            (Muon_pt[Dimu_t1muIdx[i]] > 6.) &&
                            (Muon_pt[Dimu_t2muIdx[i]] > 6.) &&
                            (Muon_mediumId[Dimu_t1muIdx[i]]) &&
                            (Muon_mediumId[Dimu_t2muIdx[i]])
                        )

                        Dimu_mass_cut.push_back(Dimu_mass[i]);
                    }
                return Dimu_mass_cut;
                },
            {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId"}); // columns 

    auto h_dimulog1_ptr = filter1.Define("hist1data", "log10(Dimu_mass_cut)")
                                 .Define("hist1weight", "2./log(10.)/Dimu_mass_cut")
                                 .Histo1D({"h_dimulog1", "h_dimulog1", nBins, x1, 2.7}, "hist1data", "hist1weight");
    auto report1 = filter1.Report();

    // run over Muonia sample
    TChain* chain = new TChain("Events");
    //chain->Add("/nfs/dust/cms/user/yangq2/eosdata/Run2011A_MuOnia_merged.root");
    chain->Add("root://eospublic.cern.ch//eos/opendata/cms/derived-data/NanoAODRun1/01-Jul-22/Run2011A_MuOnia_merged.root");

    ROOT::RDataFrame df_MuOnia(*chain);

    // this one is historic, will not be used

    auto filter4 = df_MuOnia.Filter("run < 170000", "Run number")
        .Filter("!(Alsoon_DoubleMu && Trig_DoubleMuThresh > 12)", "Dimuon threshold & sample overlap")
        .Filter("Trig_JpsiThresh != 0", "J/Psi threshold")
        .Filter("!HLT_DoubleMu4_LowMass_Displaced && !HLT_DoubleMu4p5_LowMass_Displaced && !HLT_DoubleMu5_LowMass_Displaced && !HLT_Dimuon6p5_LowMass_Displaced && !HLT_Dimuon7_LowMass_Displaced && !HLT_DoubleMu4_Jpsi_Displaced && !HLT_DoubleMu5_Jpsi_Displaced && !HLT_Dimuon6p5_Jpsi_Displaced && !HLT_Dimuon7_Jpsi_Displaced && !HLT_Mu5_L2Mu2", "HLT")
        .Define("Dimu_mass_cut", [](
            uint nDimu,
            const RVec<int>& Dimu_charge,
            const RVec<float>& Dimu_mass,
            const RVec<int>& Dimu_t1muIdx,
            const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt,
            const RVec<bool>& Muon_mediumId
            ) -> RVec<float>
            {
                RVec<float> Dimu_mass_cut;
                for (size_t i = 0; i < nDimu; i++) {

                    if  (
                            (Dimu_charge[i] == 0) &&
                            (Muon_pt[Dimu_t1muIdx[i]] > 3) &&
                            (Muon_pt[Dimu_t2muIdx[i]] > 3) &&
                            (Muon_mediumId[Dimu_t1muIdx[i]]) &&
                            (Muon_mediumId[Dimu_t2muIdx[i]]) &&
                            (Dimu_mass[i] > 2)
                        )

                        Dimu_mass_cut.push_back(Dimu_mass[i]);  
                    }
                return Dimu_mass_cut;
                },
            {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId"});

    auto h_dimulog4_ptr = filter4.Define("hist4data", "log10(Dimu_mass_cut)")
        .Define("hist4weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog4", "h_dimulog4", nBins, x1, x2}, "hist4data", "hist4weight");
    auto report4 = filter4.Report();

    // Filling a histogram with low mass displaced dimuon sample  

    auto filter14 = df_MuOnia.Filter("run < 170000", "Run number")
        .Filter("!(Alsoon_DoubleMu && Trig_DoubleMuThresh > 12)", "Dimuon threshold & sample overlap")
    // HLT_DoubleMu5_LowMass_Displaced, HLT_DoubleMu4p5_LowMass_Displaced, HLT_DoubleMu4_LowMass_Displaced are not used since they don't have any contribution to the plot
        .Filter("HLT_Dimuon7_LowMass_Displaced || HLT_Dimuon6p5_LowMass_Displaced", "HLT")
        .Define("Dimu_mass_cut", [](
            uint nDimu,
            const RVec<int>& Dimu_charge,
            const RVec<float>& Dimu_mass,
            const RVec<int>& Dimu_t1muIdx,
            const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt,
            const RVec<bool>& Muon_mediumId
            ) -> RVec<float>
            {
                RVec<float> Dimu_mass_cut;
                for (size_t i = 0; i < nDimu; i++) {

                    if  (
                            (Dimu_charge[i] == 0) &&
                            (Muon_pt[Dimu_t1muIdx[i]] > 3) &&
                            (Muon_pt[Dimu_t2muIdx[i]] > 3) &&
                            (Muon_mediumId[Dimu_t1muIdx[i]]) &&
                            (Muon_mediumId[Dimu_t2muIdx[i]]) &&
                            (Dimu_mass[i] > 2) &&
			    //adding an auxiliary higher mass cut
			    (Dimu_mass[i] < 5.1)  
                        )

                        Dimu_mass_cut.push_back(Dimu_mass[i]);  
                    }
                return Dimu_mass_cut;
                },
            {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId"});

    auto h_dimulog14_ptr = filter14.Define("hist14data", "log10(Dimu_mass_cut)")
        .Define("hist14weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog14", "h_dimulog14", nBins, x1, x2}, "hist14data", "hist14weight");
    auto report14 = filter14.Report();

    // 2: for low pT double muon, vetos 1

    auto filter2 = df_MuOnia.Filter("run < 170000", "Run number")
        .Filter("!(Alsoon_DoubleMu && Trig_DoubleMuThresh > 12)", "Dimuon threshold and sample overlap")
        .Filter("HLT_DoubleMu3_Quarkonium", "HLT")
        .Define("Dimu_mass_cut", [](
            uint nDimu,
            const RVec<int>& Dimu_charge,
            const RVec<float>& Dimu_mass,
            const RVec<int>& Dimu_t1muIdx,
            const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt,
            const RVec<bool>& Muon_mediumId
        ) -> RVec<float> {
            RVec<float> Dimu_mass_cut;
                for (size_t i = 0; i < nDimu; i++) {

                    if  (
                            (Dimu_charge[i] == 0) &&
                            (Muon_pt[Dimu_t1muIdx[i]] > 2) &&
                            (Muon_pt[Dimu_t2muIdx[i]] > 2) &&
                            (Muon_mediumId[Dimu_t1muIdx[i]]) &&
                            (Muon_mediumId[Dimu_t2muIdx[i]]) &&
                            (Dimu_mass[i] > 2)
                        )

                        Dimu_mass_cut.push_back(Dimu_mass[i]);
                    }
                return Dimu_mass_cut;
        },
        {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId"});
    auto h_dimulog2_ptr = filter2.Define("hist2data", "log10(Dimu_mass_cut)")
        .Define("hist2weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog2", "h_dimulog2", nBins, x1, x2}, "hist2data", "hist2weight");
    auto report2 = filter2.Report();


    // 6: for Upsilon (includes 2, vetos 1)

    auto filter6 = df_MuOnia.Filter("run < 170000", "Run number")
        .Filter("!(Alsoon_DoubleMu && Trig_DoubleMuThresh > 12)", "Dimuon threshold and sample overlap")
        .Define("Dimu_mass_cut", [](
            uint nDimu,
            const RVec<int>& Dimu_charge,
            const RVec<float>& Dimu_mass,
            const RVec<int>& Dimu_t1muIdx,
            const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt,
            const RVec<bool>& Muon_mediumId,
            const bool& HLT_DoubleMu3_Quarkonium,
            const bool& HLT_Dimuon0_Upsilon,
            const bool& HLT_Dimuon0_Barrel_Upsilon,
            const bool& HLT_DoubleMu3_Upsilon,
            const bool& HLT_Dimuon5_Upsilon_Barrel,
            const bool& HLT_Dimuon7_Upsilon_Barrel
        ) -> RVec<float> {
            RVec<float> Dimu_mass_cut;
                for (size_t i = 0; i < nDimu; i++) {

                    if  (
                            (HLT_DoubleMu3_Quarkonium ||
                                (
                                    (HLT_Dimuon0_Upsilon || HLT_Dimuon0_Barrel_Upsilon || HLT_DoubleMu3_Upsilon ||
                                        HLT_Dimuon5_Upsilon_Barrel || HLT_Dimuon7_Upsilon_Barrel) &&

                                    (Dimu_mass[i] > 7.) &&
                                    (Dimu_mass[i] < 14.))
                            ) &&

                            (Dimu_mass[i] > 2.) &&
                            (Dimu_charge[i] == 0) &&
                            (Muon_pt[Dimu_t1muIdx[i]] > 2.) &&
                            (Muon_pt[Dimu_t2muIdx[i]] > 2.) &&
                            (Muon_mediumId[Dimu_t1muIdx[i]]) &&
                            (Muon_mediumId[Dimu_t2muIdx[i]])
                        )

                        Dimu_mass_cut.push_back(Dimu_mass[i]);
                    }
                return Dimu_mass_cut;
        },
        {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId", "HLT_DoubleMu3_Quarkonium",
        "HLT_Dimuon0_Upsilon", "HLT_Dimuon0_Barrel_Upsilon", "HLT_DoubleMu3_Upsilon", "HLT_Dimuon5_Upsilon_Barrel", "HLT_Dimuon7_Upsilon_Barrel"});

    auto h_dimulog6_ptr = filter6.Define("hist6data", "log10(Dimu_mass_cut)")
        .Define("hist6weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog6", "h_dimulog6", nBins, x1, x2}, "hist6data", "hist6weight");
    auto report6 = filter6.Report();

    // 7: for Bs (includes 2, vetos 1)

    auto filter7 = df_MuOnia.Filter("run < 170000", "Run number")
        .Filter("!(Alsoon_DoubleMu && Trig_DoubleMuThresh>12)", "Dimuon threshold and sample overlap")
        .Define("Dimu_mass_cut", [](
            uint nDimu,
            const RVec<int>& Dimu_charge,
            const RVec<float>& Dimu_mass,
            const RVec<int>& Dimu_t1muIdx,
            const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt, 
            const RVec<bool>& Muon_mediumId,
            const bool& HLT_DoubleMu3_Quarkonium,
            const bool& HLT_Dimuon6_Bs,
            const bool& HLT_Dimuon4_Bs_Barrel,
            const bool& HLT_DoubleMu4_Dimuon6_Bs,
            const bool& HLT_DoubleMu4_Dimuon4_Bs_Barrel,
            const bool& HLT_DoubleMu3_Bs,
            const bool& HLT_DoubleMu2_Bs
        ) -> RVec<float> {
            RVec<float> Dimu_mass_cut;
                for (size_t i = 0; i < nDimu; i++) {

                    if  (
                            (
                                (HLT_DoubleMu3_Quarkonium && Muon_pt[Dimu_t1muIdx[i]] > 2. && Muon_pt[Dimu_t2muIdx[i]] > 2.) ||
                                (
                                    (HLT_Dimuon6_Bs || HLT_Dimuon4_Bs_Barrel || HLT_DoubleMu4_Dimuon6_Bs ||
                                        HLT_DoubleMu4_Dimuon4_Bs_Barrel || HLT_DoubleMu3_Bs || HLT_DoubleMu2_Bs) &&
                                    (Dimu_mass[i] > 4.) &&
                                    (Dimu_mass[i] < 7.)
                                )
                            ) &&

                            (Dimu_mass[i] > 2.) &&
                            (Dimu_charge[i] == 0) &&
                            (Muon_pt[Dimu_t1muIdx[i]] > 2.) &&
                            (Muon_pt[Dimu_t2muIdx[i]] > 2.) &&
                            (Muon_mediumId[Dimu_t1muIdx[i]]) &&
                            (Muon_mediumId[Dimu_t2muIdx[i]])
                        )

                        Dimu_mass_cut.push_back(Dimu_mass[i]);
                    }
                return Dimu_mass_cut;
        },
        {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId", "HLT_DoubleMu3_Quarkonium",
        "HLT_Dimuon6_Bs", "HLT_Dimuon4_Bs_Barrel", "HLT_DoubleMu4_Dimuon6_Bs", "HLT_DoubleMu4_Dimuon4_Bs_Barrel", "HLT_DoubleMu3_Bs", "HLT_DoubleMu2_Bs"});


    auto h_dimulog7_ptr = filter7.Define("hist7data", "log10(Dimu_mass_cut)")
        .Define("hist7weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog7", "h_dimulog7", nBins, x1, x2}, "hist7data", "hist7weight");
    auto report7 = filter7.Report();

    // 8: J/psi, includes 2, vetos 14, vetos 1 

    auto filter8 = df_MuOnia.Filter("run < 170000", "Run number")
        .Filter("!(Alsoon_DoubleMu && Trig_DoubleMuThresh > 12)", "Dimuon threshold and sample overlap")
        .Define("Dimu_mass_cut", [](
            uint nDimu,
            const RVec<int>& Dimu_charge,
            const RVec<float>& Dimu_mass,
            const RVec<int>& Dimu_t1muIdx,
            const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt,
            const RVec<bool>& Muon_mediumId,
            const bool& HLT_DoubleMu3_Quarkonium,
            const bool& HLT_Dimuon6p5_Jpsi,
            const bool& HLT_Dimuon6p5_Barrel_Jpsi,
            const bool& HLT_Dimuon0_Jpsi,
            const bool& HLT_Dimuon13_Jpsi_Barrel,
            const bool& HLT_Dimuon10_Jpsi_Barrel,
            const bool& HLT_Dimuon7_LowMass_Displaced,
            const bool& HLT_Dimuon6p5_LowMass_Displaced 
        ) -> RVec<float> {
            RVec<float> Dimu_mass_cut;
                for (size_t i = 0; i < nDimu; i++) {

                    if  (
                            (
                                (HLT_DoubleMu3_Quarkonium && Muon_pt[Dimu_t1muIdx[i]] > 3. && Muon_pt[Dimu_t2muIdx[i]] > 3.) ||
                                (
                                    (HLT_Dimuon6p5_Jpsi || HLT_Dimuon6p5_Barrel_Jpsi) &&
                                    (Dimu_mass[i] > 2.8) &&
                                    (Dimu_mass[i] < 4.3)) ||
                                (
                                    (HLT_Dimuon0_Jpsi || HLT_Dimuon13_Jpsi_Barrel || HLT_Dimuon10_Jpsi_Barrel) &&
                                    (Dimu_mass[i] > 2.8) &&
                                    (Dimu_mass[i] < 3.4))
                            ) &&

                            (Dimu_mass[i] > 2.) &&
                            (Dimu_charge[i] == 0) &&
                            (Muon_pt[Dimu_t1muIdx[i]] > 1.5) &&
                            (Muon_pt[Dimu_t2muIdx[i]] > 1.5) &&
                            (Muon_mediumId[Dimu_t1muIdx[i]]) &&
                            (Muon_mediumId[Dimu_t2muIdx[i]]) &&
			    // add veto for displaced
                            !((HLT_Dimuon7_LowMass_Displaced || HLT_Dimuon6p5_LowMass_Displaced) && 
                              (Muon_pt[Dimu_t1muIdx[i]] > 3) &&
                              (Muon_pt[Dimu_t2muIdx[i]] > 3))
                        )

                        Dimu_mass_cut.push_back(Dimu_mass[i]);
                    }
                return Dimu_mass_cut;
    },
    {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId", "HLT_DoubleMu3_Quarkonium",
	"HLT_Dimuon6p5_Jpsi", "HLT_Dimuon6p5_Barrel_Jpsi", "HLT_Dimuon0_Jpsi", "HLT_Dimuon13_Jpsi_Barrel", "HLT_Dimuon10_Jpsi_Barrel", "HLT_Dimuon7_LowMass_Displaced", "HLT_Dimuon6p5_LowMass_Displaced"});

    auto h_dimulog8_ptr = filter8.Define("hist8data", "log10(Dimu_mass_cut)")
        .Define("hist8weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog8", "h_dimulog8", nBins, x1, x2}, "hist8data", "hist8weight");

    auto report8 = filter8.Report();

    // 12: psiprime, includes 2, includes 8, vetos 1

    auto filter12 = df_MuOnia.Filter("run < 170000", "Run number")
        .Filter("!(Alsoon_DoubleMu && Trig_DoubleMuThresh > 12)", "Dimuon threshold and sample overlap")
        .Define("Dimu_mass_cut", [](
            uint nDimu,
            const RVec<int>& Dimu_charge,
            const RVec<float>& Dimu_mass,
            const RVec<int>& Dimu_t1muIdx,
            const RVec<int>& Dimu_t2muIdx,
            const RVec<float>& Muon_pt,
            const RVec<bool>& Muon_mediumId,
            const bool& HLT_DoubleMu3_Quarkonium,
            const bool& HLT_Dimuon6p5_Jpsi,
            const bool& HLT_Dimuon6p5_Barrel_Jpsi,
            const bool& HLT_Dimuon0_Jpsi,
            const bool& HLT_Dimuon13_Jpsi_Barrel,
            const bool& HLT_Dimuon10_Jpsi_Barrel,
            const bool& HLT_Dimuon11_PsiPrime,
            const bool& HLT_Dimuon9_PsiPrime,
            const bool& HLT_Dimuon7_PsiPrime
        ) -> RVec<float> {
            RVec<float> Dimu_mass_cut;
                for (size_t i = 0; i < nDimu; i++) {

                    if  (
                            (
                                (HLT_DoubleMu3_Quarkonium && Muon_pt[Dimu_t1muIdx[i]] > 3. && Muon_pt[Dimu_t2muIdx[i]] > 3.) ||
                                (
                                    (HLT_Dimuon6p5_Jpsi || HLT_Dimuon6p5_Barrel_Jpsi) &&
                                    (Dimu_mass[i] > 2.5) &&
                                    (Dimu_mass[i] < 4.3)) ||
                                (
                                    (HLT_Dimuon0_Jpsi || HLT_Dimuon13_Jpsi_Barrel || HLT_Dimuon10_Jpsi_Barrel) &&
                                    (Dimu_mass[i] > 2.8) &&
                                    (Dimu_mass[i] < 3.4)) ||
                                (
                                    (HLT_Dimuon11_PsiPrime || HLT_Dimuon9_PsiPrime || HLT_Dimuon7_PsiPrime) &&
                                    (Dimu_mass[i] > 3.4) &&
                                    (Dimu_mass[i] < 4.3))
                            ) &&

                            (Dimu_mass[i] > 2.) &&
                            (Dimu_charge[i] == 0) &&
                            (Muon_pt[Dimu_t1muIdx[i]] > 1.5) &&
                            (Muon_pt[Dimu_t2muIdx[i]] > 1.5) &&
                            (Muon_mediumId[Dimu_t1muIdx[i]]) &&
                            (Muon_mediumId[Dimu_t2muIdx[i]])
                        )

                        Dimu_mass_cut.push_back(Dimu_mass[i]);
                    }
                return Dimu_mass_cut;
    },
    {"nDimu", "Dimu_charge", "Dimu_mass", "Dimu_t1muIdx", "Dimu_t2muIdx", "Muon_pt", "Muon_mediumId", "HLT_DoubleMu3_Quarkonium",
        "HLT_Dimuon6p5_Jpsi", "HLT_Dimuon6p5_Barrel_Jpsi", "HLT_Dimuon0_Jpsi", "HLT_Dimuon13_Jpsi_Barrel", "HLT_Dimuon10_Jpsi_Barrel",
        "HLT_Dimuon11_PsiPrime", "HLT_Dimuon9_PsiPrime", "HLT_Dimuon7_PsiPrime"});

    auto h_dimulog12_ptr = filter12.Define("hist12data", "log10(Dimu_mass_cut)")
        .Define("hist12weight", "2./log(10.)/Dimu_mass_cut")
        .Histo1D({"h_dimulog12", "h_dimulog12", nBins, x1, x2}, "hist12data", "hist12weight");
    auto report12 = filter12.Report();

    // Event loops is run here (once each)
    TH1D h_dimulog1  = *h_dimulog1_ptr;
    TH1D h_dimulog4  = *h_dimulog4_ptr;
    TH1D h_dimulog2  = *h_dimulog2_ptr;
    TH1D h_dimulog6  = *h_dimulog6_ptr;
    TH1D h_dimulog7  = *h_dimulog7_ptr;
    TH1D h_dimulog8  = *h_dimulog8_ptr;
    TH1D h_dimulog12 = *h_dimulog12_ptr;
    TH1D h_dimulog14 = *h_dimulog14_ptr;  // histogram with low mass displaced dimuon sample 

    h_dimulog1.SetDirectory(0);

    TH1D* h_dimulog1_  = &h_dimulog1;
    TH1D* h_dimulog2_  = &h_dimulog2;
    TH1D* h_dimulog3_  = new TH1D("h_dimulog3_", "h_dimulog3_", nBins, x1, x2);
    TH1D* h_dimulog4_  = &h_dimulog4;
    TH1D* h_dimulog5_  = new TH1D("h_dimulog5_", "h_dimulog5_", nBins, x1, x2);
    TH1D* h_dimulog6_  = &h_dimulog6;
    TH1D* h_dimulog7_  = &h_dimulog7;
    TH1D* h_dimulog8_  = &h_dimulog8;
    TH1D* h_dimulog9_  = new TH1D("h_dimulog9_", "h_dimulog9_", nBins, x1, x2);
    TH1D* h_dimulog10_ = new TH1D("h_dimulog10_", "h_dimulog10_", nBins, x1, x2);
    TH1D* h_dimulog11_ = new TH1D("h_dimulog11_", "h_dimulog11_", nBins, x1, x2);
    TH1D* h_dimulog12_ = &h_dimulog12;
    TH1D* h_dimulog13_ = new TH1D("h_dimulog13_", "h_dimulog13_", nBins, x1, x2);
    TH1D* h_dimulog14_ = &h_dimulog14;
    TH1D* h_dimulog15_ = new TH1D("h_dimulog15_", "h_dimulog15_", nBins, x1, x2);    // histogram with low mass displaced dimuon sample 


    h_dimulog3_->Add(h_dimulog1_, h_dimulog2_, 1, 1);     //  3 = 1 + (2-1)
    h_dimulog5_->Add(h_dimulog1_, h_dimulog4_, 1, 1);     //  5 = 1 + (4-1)   (not used)
    h_dimulog9_->Add(h_dimulog1_, h_dimulog6_, 1, 1);     //  9 = 1 + (6-1) 
    h_dimulog10_->Add(h_dimulog1_, h_dimulog7_, 1, 1);    // 10 = 1 + (7-1)
    h_dimulog11_->Add(h_dimulog1_, h_dimulog8_, 1, 1);    
    h_dimulog11_->Add(h_dimulog11_, h_dimulog14_, 1, 1);  // 11 = 1 + (8-2-14) + (14-1)
    h_dimulog13_->Add(h_dimulog1_, h_dimulog12_, 1, 1);   // 13 = 1 + (12-1)
    h_dimulog15_->Add(h_dimulog1_, h_dimulog14_, 1, 1);   // histogram with low mass displaced dimuon sample  
    h_dimulog15_->Add(h_dimulog15_, h_dimulog2_, 1, 1);   // Adding the dark grey contribution to the yellow part (1+(2-1)+14) 
    

    TCanvas *c = new TCanvas("c1", "c1", 1);
    c->SetLogy();
    gStyle->SetOptStat(0); // remove box
    gStyle->SetOptTitle(0); //modified to remove title 

    // Histogram 9 = Histogram 1 + Histogram 6    high pT + Upsilon (green)
    h_dimulog9_->GetXaxis()->SetTitle("Dimuon mass [Log_{10}(m/GeV/c^{2})]");  // x axis label
    h_dimulog9_->GetXaxis()->SetTitleOffset(1.35);  // Offset the x axis title
    h_dimulog9_->GetYaxis()->SetTitle("Number of Events/10 MeV");
    h_dimulog9_->SetMinimum(0.02);
    h_dimulog9_->SetMaximum(3.E6);
    h_dimulog9_->SetFillColor(8); // green
    h_dimulog9_->Draw("hist same");


    // Histogram 10 = Histogram 1 + Histogram 7    high pT + Bs  (light blue)
    h_dimulog10_->SetFillColor(7); // blue_green  // color changed from 29 to 7
    h_dimulog10_->Draw("hist same");

    // Histogram 13 = Histogram 1 + Histogram 12   high pT + psiprime  (dark blue)
    h_dimulog13_->SetFillColor(4); // dark blue  // color changed from 9 to 4 
    h_dimulog13_->Draw("hist same");

    // Histogram 11 = Histogram 1 + Histogram 8 + 14   high pT + J/psi (red)
    h_dimulog11_->SetFillColor(2); // red
    h_dimulog11_->Draw("hist same");

    // Histogram 15 = Histogram 1 + 2 + Histogram 14   // modified to make a new histogram with displaced dimuon sample 
    h_dimulog15_->SetFillColor(kOrange-4); // Yellow
    h_dimulog15_->Draw("hist same");

    // Histogram 3 = Histogram 1 + Histogram 2     high pT + low pT (dark grey)
    h_dimulog3_->SetFillColor(39); // dark grey  // color changed from 38 to 39
    h_dimulog3_->Draw("hist same");

    // Histogram 1                                 high pT (light grey)
    TH1D* h_dimulog1_2 = (TH1D*)h_dimulog1_->Clone();
    h_dimulog1_2->SetFillColor(17); // light grey  //color changed from 18 to 17
    h_dimulog1_2->Draw("hist same");


    gPad->RedrawAxis();

    // Adding legends
    TLegend *leg = new TLegend(0.62, 0.62, 0.89, 0.89);
    leg->AddEntry(h_dimulog13_, "#Psi'", "f");
    leg->AddEntry(h_dimulog11_, "J/#Psi", "f");
    leg->AddEntry(h_dimulog10_, "B_{S}#rightarrow#mu^{+}#mu^{-}", "f");
    leg->AddEntry(h_dimulog9_, "#Upsilon", "f");
    leg->AddEntry(h_dimulog15_, "low mass displaced #mu^{+}#mu^{-} ", "f");
    leg->AddEntry(h_dimulog3_, "low p_{T} double muon", "f");
    leg->AddEntry(h_dimulog1_2, "high p_{T} double muon", "f");
    leg->SetBorderSize(0);
    leg->Draw();


    // Adding label for information about plot
    TLatex *label1 = new TLatex();
    label1->SetNDC();
    label1->SetTextSize(0.0305);
    label1->DrawLatex(0.130, 0.86, "2011 Run, L = 1.2 fb^{-1}");
    label1->DrawLatex(0.130, 0.82, "CMS Open Data");
    label1->DrawLatex(0.130, 0.78, "#sqrt{s} = 7 TeV");

    // Adding labels for the mass peaks
    TLatex *label2 = new TLatex();
    label2->SetTextAlign(22);
    label2->SetTextSize(0.04);
    label2->DrawLatex(-0.12, 4.0e3, "#omega");
    label2->DrawLatex(0.02, 4.2e3, "#phi");
    label2->DrawLatex(0.6, 7.0e5, "J/#psi");
    label2->DrawLatex(0.62, 9.0e4, "#psi'");
    label2->DrawLatex(0.73, 2.8e4, "B_{s}");
    label2->DrawLatex(1.01, 8.0e4, "#Upsilon");
    label2->DrawLatex(1.96, 1.3e3, "Z");

    c->Print("Dimuon2011_eospublic_RDF2.png");
    c->SaveAs("Dimuon2011_eospublic_RDF2.pdf");

    TFile Dimuon2011("Dimuon2011_eospublic_RDF2.root", "RECREATE"); 
    h_dimulog1_->Write();
    h_dimulog2_->Write();
    h_dimulog3_->Write();
    h_dimulog4_->Write();
    h_dimulog5_->Write();
    h_dimulog6_->Write();
    h_dimulog7_->Write();
    h_dimulog8_->Write();
    h_dimulog9_->Write();
    h_dimulog10_->Write();
    h_dimulog11_->Write();
    h_dimulog12_->Write();
    h_dimulog13_->Write();
    h_dimulog14_->Write();  
    h_dimulog15_->Write();  
    

    report1->Print();
    report2->Print();
    report4->Print();
    report6->Print();
    report7->Print();
    report8->Print();
    report12->Print();
    report14->Print();

    auto end = chrono::steady_clock::now();
    std::cout << "Elapsed time in seconds: "
        << chrono::duration_cast<chrono::seconds>(end - start).count()
        << " sec" << std::endl;
 
    std::cout << "Elapsed time in minutes: "
        << chrono::duration_cast<chrono::minutes>(end - start).count()
        << " min" << std::endl;

    std::cout << "Number of threads: " 
        << nThreads << std::endl;

    
}


int main() {

    Dimuon2011_eospublic_RDF2();

    return 0;
}
