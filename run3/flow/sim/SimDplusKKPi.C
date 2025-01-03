#include "TH1F.h"
#include "TTree.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TCanvas.h"

#include <array>
#include <map>
#include <string>
#include <vector>

#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom3.h>

#include "Pythia8/Pythia.h"
using namespace Pythia8;

void SimDplusKKPi(int nevents=10000, int seed = 1233) { 

    /*  
        Script that computes the invariant mass spectrum obtained
        from D+->PhiPi->KKPi and D+->KK(892)->KKPi with the K 
        misidentified as a pion

        Simulated 200M events with seed from 1 to 1000
    */
    Pythia8::Pythia pythia;
    
    pythia.readString(Form("Random:seed %d", seed));
    pythia.readString("Random:setSeed = on");
    pythia.readString("Tune:pp = 14");
    pythia.readString("SoftQCD:nonDiffractive = on");
    pythia.readString("411:onMode = off");
    pythia.readString("333:onMode = off");
    pythia.readString("313:onMode = off");
    pythia.readString("411:onIfMatch = 333 211");     // D+ -> Phi Pi+
    pythia.readString("333:onIfMatch = 321 -321");    // Phi -> K+ K-
    pythia.readString("411:onIfMatch = 321 -313");    // D+ -> K+ K0(892)bar
    pythia.readString("313:onIfMatch = 321 -211");    // K0(892) -> K+ Pi-

    int DplusPDG = 411;
    pythia.init();

    TH1F * hInvMass = new TH1F("hInvMass", ";M(#pi K#pi) (GeV/c^{2});Counts", 1500, 1.3, 2.5);
    TH1F * hPt = new TH1F("hPt", ";p_{T} (GeV/c);Counts", 200, 0, 50);
    TH1F * hDecayFlag = new TH1F("hDecayFlag", ";Decay flag;Counts", 3, 0, 3);
    
    // Create a ROOT Tree
    double invMass;
    double pt;
    int decayFlag; // 1 for Phi, 2 for K(892)
    TTree tree("DplusKKPiMisidKPi","DplusKKPiMisidKPi");
    tree.Branch("mass",&invMass,"M(#pi K#pi)/D");    
    tree.Branch("pt",&pt,"p_{T}/D");    
    tree.Branch("decayflag",&decayFlag,"Decay Flag/I");    

    double MPi = 0.13957;
    for (int iEvent = 0; iEvent < nevents; ++iEvent) {
        if (!pythia.next()) {
            continue;
        }
        for (int iPart = 2; iPart < pythia.event.size(); iPart++) {
            Particle part = pythia.event.at(iPart);
            if(abs(part.id()) == DplusPDG) {
                Particle KPlus, KMinus, PiPlus;
                decayFlag = 0;
                std::vector<int> allDaughtersIdxs = part.daughterListRecursive();
                for(int iFirstPart=0; iFirstPart<allDaughtersIdxs.size(); iFirstPart++) {
                    Particle daug = pythia.event.at(allDaughtersIdxs[iFirstPart]); 
                    if(abs(pythia.event.at(allDaughtersIdxs[iFirstPart]).id()) == 333) {
                        decayFlag=1;
                    }
                    if(abs(pythia.event.at(allDaughtersIdxs[iFirstPart]).id()) == 313) {
                        decayFlag=2;
                    }
                    if(part.id() == DplusPDG) {  // particle
                        if(pythia.event.at(allDaughtersIdxs[iFirstPart]).id() == 321)  KPlus=daug;
                        if(pythia.event.at(allDaughtersIdxs[iFirstPart]).id() == -321) KMinus=daug;
                        if(pythia.event.at(allDaughtersIdxs[iFirstPart]).id() == 211)  PiPlus=daug;
                    } else {                  // antiparticle
                        if(pythia.event.at(allDaughtersIdxs[iFirstPart]).id() == -321) KPlus=daug;
                        if(pythia.event.at(allDaughtersIdxs[iFirstPart]).id() == 321)  KMinus=daug;
                        if(pythia.event.at(allDaughtersIdxs[iFirstPart]).id() == -211) PiPlus=daug;
                    }
                }

                // The KPlus can be wrongly identified as a PiPlus, this changes its energy
                double EKPlus = TMath::Sqrt(MPi*MPi + KPlus.px()*KPlus.px() + 
                                            KPlus.py()*KPlus.py() + KPlus.pz()*KPlus.pz());
                double totE = EKPlus + KMinus.e() + PiPlus.e();
                double totPx = KPlus.px() + KMinus.px() + PiPlus.px();
                double totPy = KPlus.py() + KMinus.py() + PiPlus.py();
                double totPz = KPlus.pz() + KMinus.pz() + PiPlus.pz();
                TLorentzVector wrongRecoDPlus(totPx, totPy, totPz, totE);
                invMass = wrongRecoDPlus.M();
                pt = wrongRecoDPlus.Pt();
                hInvMass->Fill(invMass);
                hPt->Fill(pt);
                hDecayFlag->Fill(decayFlag);
                tree.Fill();
            }
        }
    }
    
    TFile oFile(Form("./outputDplus/DplusToKKPiWrongRecoPiKPi_%i.root", seed), "recreate"); 
    oFile.cd();
    hInvMass->Write();
    hPt->Write();
    hDecayFlag->Write();
    tree.Write();
    oFile.Close();      
    
}