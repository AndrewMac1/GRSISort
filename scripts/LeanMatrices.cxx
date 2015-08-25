
//g++ exAnalysis.C -std=c++0x -I$GRSISYS/include -L$GRSISYS/libraries -lAnalysisTreeBuilder -lGriffin -lSceptar -lDescant -lPaces -lGRSIDetector -lTGRSIFit -lTigress -lSharc -lCSM -lTriFoil -lTGRSIint -lGRSILoop -lMidasFormat -lGRSIRootIO -lDataParser -lGRSIFormat -lMidasFormat -lXMLParser -lXMLIO -lProof -lGuiHtml `grsi-config --cflags --libs` `root-config --cflags --libs`  -lTreePlayer -lGROOT -lX11 -lXpm -lSpectrum
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cstdio>
#include <sys/stat.h>

#include "TROOT.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TVirtualIndex.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TGRSIRunInfo.h"
#include "TGRSISortInfo.h"
#include "Globals.h"

#ifndef __CINT__ 
#include "TGriffin.h"
#include "TSceptar.h"
#endif

//This code is an example of how to write an analysis script to analyse an analysis tree
//The compiled script works like this
//
//  1. Starts in the main function by finding the analysis tree, setting up input and output files
//  2. Calls the exAnalysis function
//  3. exAnalysis creates 1D and 2D histograms and adds them to a list
//  4. Some loops over event "packets" decides which histograms should be filled and fills them.
//  5. The list of now filled histograms is returned to the main function
//  6. The list is written (meaning all of the histograms are written) to the output root file
//  7. Papers are published, theses are granted, high-fives are made
//
/////////////////////////////////////////////////////////////////////////////////////////

//This function gets run if running interpretively
//Not recommended for the analysis scripts
#ifdef __CINT__ 
void exAnalysis() {
   if(!AnalysisTree) {
      printf("No analysis tree found!\n");
      return;
   }
   //coinc window = 0-20, bg window 40-60, 6000 bins from 0. to 6000. (default is 4000)
   TList *list = exAnalysis(AnalysisTree, 0.);

   TFile *outfile = new TFile("output.root","recreate");
   list->Write();
}
#endif

std::vector<std::pair<double,int> > AngleCombinations(double distance = 110., bool folding = false) {
   std::vector<std::pair<double,int> > result;
   std::vector<std::pair<double,int> > grouped_result;
   std::vector<double> angle;
   for(int firstDet = 1; firstDet <= 12; ++firstDet) {
      for(int firstCry = 0; firstCry < 4; ++firstCry) {
         for(int secondDet = 1; secondDet <= 12; ++secondDet) {
            for(int secondCry = 0; secondCry < 4; ++secondCry) {
               if(firstDet == secondDet && firstCry == secondCry) {
                  continue;
               }
               angle.push_back(TGriffin::GetPosition(firstDet, firstCry, distance).Angle(TGriffin::GetPosition(secondDet, secondCry, distance))*180./TMath::Pi());
              // if(folding && angle.back() > 90.) {
              //                 //  angle.back() = 180. - angle.back();
              //                               // }
              //                                             // std::cout<<firstDet<<", "<<firstCry<<" - "<<secondDet<<", "<<secondCry<<": "<<angle.back()<<std::endl;
               }
            }
        }
   }

   std::sort(angle.begin(),angle.end());
   size_t r;
   for(size_t a = 0; a < angle.size(); ++a) {
      for(r = 0; r < result.size(); ++r) {
         if(angle[a] >= result[r].first-0.001 && angle[a] <= result[r].first+0.001) {
            (result[r].second)++;
            break;
         }
      }
      if(result.size() == 0 || r == result.size()) {
         result.push_back(std::make_pair(angle[a],1));
      }
   }

   if(folding) {//if we fold we also want to group
     std::vector<std::pair<double,int> > groupedResult;
     for(size_t i = 0, r = 0; i < result.size(); ++i) {
       switch(i) {
       case 0:
       case 1:
         groupedResult.push_back(result[i]);
         break;
       case 2:
       case 4:
       case 6:
         if(i+1 >= result.size()) {
           std::cerr<<"Error!"<<std::endl;
         }
         groupedResult.push_back(std::make_pair((result[i].first+result[i+1].first)/2.,result[i].second+result[i+1].second));
         ++i;
         break;
       default:
         groupedResult.push_back(std::make_pair((result[i].first+result[i+1].first+result[i+2].first)/3.,result[i].second+result[i+1].second+result[i+2].second));
         i+=2;
         break;
       }
     }
     return groupedResult;
   }

                 return result;

}





TList *exAnalysis(TTree* tree, long maxEntries = 0, TStopwatch* w = NULL) {

   ///////////////////////////////////// SETUP ///////////////////////////////////////
   //Histogram paramaters
   Double_t low = 0;
   Double_t high = 10000;
   Double_t nofBins = 10000;

   //Coincidence Parameters
   Double_t ggTlow = 0.;   //Times are in 10's of ns
   Double_t ggThigh = 60.;
   Double_t gbTlow =  0.;
   Double_t gbThigh = 100.;

   Double_t ggBGlow = 100.;
   Double_t ggBGhigh = 175.;
   Double_t gbBGlow = -160.;
   Double_t gbBGhigh = 0.;
   Double_t ggBGScale = (ggThigh - ggTlow)/(ggBGhigh - ggBGlow);
   Double_t gbBGScale = (gbThigh - gbTlow)/(gbBGhigh - gbBGlow);

   if(w == NULL) {
      w = new TStopwatch;
      w->Start();
   }
   TList* list = new TList;

   const size_t MEM_SIZE = (size_t)1024*(size_t)1024*(size_t)1024*(size_t)8; // 8 GB

   std::vector<std::pair<double,int> > angleCombinations = AngleCombinations(110., false);
   std::cout<<"got "<<angleCombinations.size()<<" angles"<<std::endl;
   for(auto ang = angleCombinations.begin(); ang != angleCombinations.end(); ang++) {
      std::cout<<(*ang).first<<" degree: "<<(*ang).second<<" combinations"<<std::endl;
   }

   double* xBins = new double[1501];
   double* yBins = new double[1501];
   for(int i = 0; i <= 1500; ++i) {
      xBins[i] = (double) i;
      yBins[i] = (double) i;
   }

   double* zBins = new double[angleCombinations.size()+1];
   for(int i = 0; i<=52; i++){
        zBins[i] = (double) i;
   }
   //We create some spectra and then add it to the list
   TH1F* gammaSingles = new TH1F("gammaSingles","#gamma singles;energy[keV]",nofBins, low, high); list->Add(gammaSingles);
   TH1F* gammaSinglesB = new TH1F("gammaSinglesB","#beta #gamma;energy[keV]",nofBins, low, high); list->Add(gammaSinglesB);
   TH1F* gammaSinglesBt = new TH1F("gammaSinglesBt","#beta #gamma t-rand-corr; energy[keV]",nofBins, low, high); list->Add(gammaSinglesBt);
   TH1F* ggTimeDiff = new TH1F("ggTimeDiff", "#gamma-#gamma time difference", 300,0,300); list->Add(ggTimeDiff);
   TH1F* gbTimeDiff = new TH1F("gbTimeDiff", "#gamma-#beta time difference", 2000,-1000,1000); list->Add(gbTimeDiff); 
   TH1F* gtimestamp = new TH1F("gtimestamp", "#gamma time stamp", 10000,0,1000); list->Add(gtimestamp);
   TH1F* btimestamp = new TH1F("btimestamp", "#beta time stamp", 10000,0,1000); list->Add(btimestamp);
   TH2F* gbEnergyvsgTime = new TH2F("gbEnergyvsgTime", "#gamma #beta coincident: #gamma timestamp vs. #gamma energy; Time [s]; Energy [keV]", 1000,0,1000, nofBins, low, high); list->Add(gbEnergyvsgTime);
   TH2F* gbEnergyvsbTime = new TH2F("gbEnergyvsbTime", "#gamma #beta coincident: #beta timestamp vs. #gamma energy; Time [s]; Energy [keV]", 1000,0,1000, nofBins, low, high); list->Add(gbEnergyvsbTime);
   TH2F* ggmatrix = new TH2F("ggmatrix","#gamma-#gamma matrix",nofBins, low, high,nofBins, low, high); list->Add(ggmatrix);
   TH2F* ggmatrixt = new TH2F("ggmatrixt","#gamma-#gamma matrix t-corr",nofBins,low,high,nofBins,low,high); list->Add(ggmatrixt);
   TH2F* gammaSinglesB_hp = new TH2F("gammaSinglesB_hp", "#gamma-#beta vs. SC channel", nofBins,low,high,20,1,21); list->Add(gammaSinglesB_hp);
   TH2F* ggbmatrix = new TH2F("ggbmatrix","#gamma-#gamma-#beta matrix", nofBins, low, high, nofBins, low, high); list->Add(ggbmatrix);
   TH2F* ggbmatrixt = new TH2F("ggbmatrixt","#gamma-#gamma-#beta matrix t-corr", nofBins, low, high, nofBins, low, high); list->Add(ggbmatrixt);
   TH2F* ggEcosTheta1408 = new TH2F("ggEcosTheta1408","#gamma-#gamma Energy cosTheta index ", 5200, low, 2600, 52, low, 51); list->Add(ggEcosTheta1408);
   TH2F* ggEcosTheta1112 = new TH2F("ggEcosTheta1112","#gamma-#gamma Energy cosTheta index ", 5200, low, 2600, 52, low, 51); list->Add(ggEcosTheta1112);
   TH2F* ggEcosTheta964 = new TH2F("ggEcosTheta964","#gamma-#gamma Energy cosTheta index ", 5200, low, 2600, 52, low, 51); list->Add(ggEcosTheta964);
   TH2F* ggEcosTheta778 = new TH2F("ggEcosTheta778","#gamma-#gamma Energy cosTheta index ", 5200, low, 2600, 52, low, 51); list->Add(ggEcosTheta778);
   TH2F* ggEcosTheta1089 = new TH2F("ggEcosTheta1089","#gamma-#gamma Energy cosTheta index ", 5200, low, 2600, 52, low, 51); list->Add(ggEcosTheta1089);
   TH2F* ggEcosTheta411 = new TH2F("ggEcosTheta411","#gamma-#gamma Energy cosTheta index ", 5200, low, 2600, 52, low, 51); list->Add(ggEcosTheta411);
   TH2F *ggEcosThetaIso = new TH2F("ggEcosThetaIso","#gamma-#gamma Energy cosTheta Isotropic", 5200, low, 2600, 52, low, 51); list->Add(ggEcosThetaIso);
   TH2F* ggmatrixIndex[52];
   TH2F* ggmatrixIndexIso[52];
   for(int i = 0; i <51; i++){
        ggmatrixIndex[i] = new TH2F(Form("ggmatrixIndex_%d",i),Form("#gamma-#gamma Energy cosTheta index %d",i), 5200, low, 2600, 5200, low, 2600); list->Add(ggmatrixIndex[i]);
        ggmatrixIndexIso[i] = new TH2F(Form("ggmatrixIndexIso_%d",i),Form("#gamma-#gamma Energy cosTheta Isotropic index %d",i), 5200, low, 2600, 5200, low, 2600); list->Add(ggmatrixIndexIso[i]);
   }
   TH3F* angCorr_coinc_Binned = new TH3F("angCorr_coinc_Binned","angular correlation cube;energy [keV];energy [keV];angle [^{o}]", 1500, xBins, 1500, yBins, angleCombinations.size(), zBins); list->Add(angCorr_coinc_Binned);
   TH2F* grifscep_hp = new TH2F("grifscep_hp","Sceptar vs Griffin hit pattern",64,0,64,20,0,20); list->Add(grifscep_hp);
   TH2F* gbTimevsg = new TH2F("gbTimevsg","#gamma energy vs. #gamma-#beta timing",300,-150,150,nofBins,low,high); list->Add(gbTimevsg); 
   TH2F* gammaSinglesCyc = new TH2F("gammaSinglesCyc", "Cycle time vs. #gamma energy", 25000.,0.,25000., 4500,low,4500); list->Add(gammaSinglesCyc);
   TH2F* gammaSinglesBCyc = new TH2F("gammaSinglesBCyc", "Cycle time vs. #beta coinc #gamma energy", 25000.,0.,25000.,4500,low,4500); list->Add(gammaSinglesBCyc);
   ///////////////////////////////////// PROCESSING /////////////////////////////////////

   //set up branches
   //Each branch can hold multiple hits
   //ie TGriffin grif holds 3 gamma rays on a triples event 
   TGriffin* grif = 0;
   TSceptar* scep = 0;
   tree->SetBranchAddress("TGriffin", &grif); //We assume we always have a Griffin

   bool gotSceptar;
   if(tree->FindBranch("TSceptar") == 0) {   //We check to see if we have a Scepter branch in the analysis tree
      gotSceptar = false;
   } else {
      tree->SetBranchAddress("TSceptar", &scep);
      gotSceptar = true;
   }

   //tree->LoadBaskets(MEM_SIZE);   

   long entries = tree->GetEntries();

   //These are the indices of the two hits being compared
   int one;
   int two;
   int three;
   int index = 0;
   std::cout<<std::fixed<<std::setprecision(1); //This just make outputs not look terrible
   int entry;
   size_t angIndex;
   if(maxEntries == 0 || maxEntries > tree->GetEntries()) {
      maxEntries = tree->GetEntries();
   }
   for(entry = 1; entry < maxEntries; ++entry) { //Only loop over the set number of entries
                                                 //I'm starting at entry 1 because of the weird high stamp of 4.
      bool iso_flag = false;
      if(entry > 2)tree->GetEntry(entry-2);
        if(grif->GetMultiplicity()>1)
                iso_flag = true;

      TGriffin* grif2 = (TGriffin*)(grif->Clone()); 
     tree->GetEntry(entry);
      
      //loop over the gammas in the event packet
      //grif is the variable which points to the current TGriffin
      for(one = 0; one < (int) grif->GetMultiplicity(); ++one) {
         //We want to put every gamma ray in this event into the singles
         gammaSingles->Fill(grif->GetGriffinHit(one)->GetEnergy()); 
              for(three = 0; three < (int) grif2->GetMultiplicity(); ++three) {
                   double ang2 = grif->GetGriffinHit(one)->GetPosition().Angle(grif2->GetGriffinHit(three)->GetPosition())*180./TMath::Pi();
                         for(int i = 0; i < 52; i++){
                                if(ang2 >= angleCombinations[i].first - 0.0005 && ang2 <= angleCombinations[i].first + 0.0005){
                                        ggmatrixIndexIso[i]->Fill(grif->GetGriffinHit(one)->GetEnergy(),grif2->GetGriffinHit(three)->GetEnergy());
                                        if(grif->GetGriffinHit(one)->GetEnergy() > 1405.5 && grif->GetGriffinHit(one)->GetEnergy() < 1411.5){
                                                 if(iso_flag && grif2->GetMultiplicity() > 0.){
                                                         ggEcosThetaIso->Fill(grif2->GetGriffinHit(three)->GetEnergy(), i);
                                                }//TMath::Cos(TMath::DegToRad()*ang));
                                                 }
                                              //   if(grif2->GetGriffinHit(three)->GetEnergy() < 1334.5 && grif2->GetGriffinHit(three)->GetEnergy() > 1330.5){
                                                //        if(iso_flag && grif2->GetMultiplicity() > 0.){
                                                  //              ggEcosThetaIso->Fill(grif->GetGriffinHit(one)->GetEnergy(), i);
                                               // }//TMath::Cos(TMath::DegToRad()*ang));
                                                  //      }
                                                        break;
                                                   }
                                        }
              }         
         gtimestamp->Fill(grif->GetGriffinHit(one)->GetTime()/100000000.);
         Long_t time = (Long_t)(grif->GetHit(one)->GetTime());
    //     time = time%2268500000L;
         time = time%957500000L;
         gammaSinglesCyc->Fill(grif->GetCycleTimeInMilliSeconds(grif->GetHit(one)->GetTime()/10.), grif->GetHit(one)->GetEnergy());
         //We now want to loop over any other gammas in this packet
         for(two = 0; two < (int) grif->GetMultiplicity(); ++two) {
            if(two == one){ //If we are looking at the same gamma we don't want to call it a coincidence
               continue;
            }
	    double ang = grif->GetGriffinHit(one)->GetPosition().Angle(grif->GetGriffinHit(two)->GetPosition())*180./TMath::Pi();
            for(angIndex = 0; angIndex < angleCombinations.size(); ++angIndex) {
                if(TMath::Abs(ang - angleCombinations[angIndex].first) < 1.0) {
                  break;
               }
            }
            if(ggTlow <= TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) && TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) < ggThigh) {
               for(int i = 0; i < 52; i++){
                   if(ang >= angleCombinations[i].first - 0.0005 && ang <= angleCombinations[i].first + 0.0005){
                     index = i;
                     angCorr_coinc_Binned->Fill(grif->GetGriffinHit(one)->GetEnergy(),grif->GetGriffinHit(two)->GetEnergy(),(double) index,1.);
                        break;
                   }
               }
                for(int i = 0; i < 52; i++){
                   if((ang >= (angleCombinations[i].first - 0.0005)) && (ang <= (angleCombinations[i].first + 0.0005))){
                        ggmatrixIndex[i]->Fill(grif->GetGriffinHit(one)->GetEnergy(),grif->GetGriffinHit(two)->GetEnergy());
                        if(grif->GetGriffinHit(one)->GetEnergy() > 1405.5 && grif->GetGriffinHit(one)->GetEnergy() < 1411.5){
                                ggEcosTheta1408->Fill(grif->GetGriffinHit(two)->GetEnergy(), i);//TMath::Cos(TMath::DegToRad()*ang));
                        }
                        if(grif->GetGriffinHit(one)->GetEnergy() > 1108.5 && grif->GetGriffinHit(one)->GetEnergy() < 1115.5){
                                ggEcosTheta1112->Fill(grif->GetGriffinHit(two)->GetEnergy(), i);//TMath::Cos(TMath::DegToRad()*ang));
                        }
                        if(grif->GetGriffinHit(one)->GetEnergy() > 960.5 && grif->GetGriffinHit(one)->GetEnergy() < 967.5){
                                ggEcosTheta964->Fill(grif->GetGriffinHit(two)->GetEnergy(), i);//TMath::Cos(TMath::DegToRad()*ang));
                        }
                        if(grif->GetGriffinHit(one)->GetEnergy() > 774.5 && grif->GetGriffinHit(one)->GetEnergy() < 782.5){
                                ggEcosTheta778->Fill(grif->GetGriffinHit(two)->GetEnergy(), i);//TMath::Cos(TMath::DegToRad()*ang));
                        }
                        if(grif->GetGriffinHit(one)->GetEnergy() > 1085.5 && grif->GetGriffinHit(one)->GetEnergy() < 1093.5){
                                ggEcosTheta1089->Fill(grif->GetGriffinHit(two)->GetEnergy(), i);//TMath::Cos(TMath::DegToRad()*ang));
                        }
                        if(grif->GetGriffinHit(one)->GetEnergy() > 407.5 && grif->GetGriffinHit(one)->GetEnergy() < 414.5){
                                ggEcosTheta411->Fill(grif->GetGriffinHit(two)->GetEnergy(), i);//TMath::Cos(TMath::DegToRad()*ang));
			 }
                      //  if(grif->GetGriffinHit(two)->GetEnergy() < 1334.5 && grif->GetGriffinHit(two)->GetEnergy() > 1330.5){
                         //       ggEcosTheta->Fill(grif->GetGriffinHit(one)->GetEnergy(), i);//TMath::Cos(TMath::DegToRad()*ang));

                      //  }
                        break;

                      }
                    }
               }
            //Check to see if the two gammas are close enough in time
            ggTimeDiff->Fill(TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()));
            if(ggTlow <= TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) && TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) < ggThigh) { 
               //If they are close enough in time, fill the gamma-gamma matrix. This will be symmetric because we are doing a double loop over gammas
               ggmatrix->Fill(grif->GetGriffinHit(one)->GetEnergy(), grif->GetGriffinHit(two)->GetEnergy());
            }
            if(ggBGlow <= TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) && TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) < ggBGhigh) { 
               //If they are not close enough in time, fill the time-random gamma-gamma matrix. This will be symmetric because we are doing a double loop over gammas
               ggmatrixt->Fill(grif->GetGriffinHit(one)->GetEnergy(), grif->GetGriffinHit(two)->GetEnergy());
            }
         }
      }
      

      //Now we make beta gamma coincident matrices
      if(gotSceptar && scep->GetMultiplicity() > 0) {
        //We do an outside loop on gammas so that we can break on the betas if we see a beta in coincidence (we don't want to bin twice just because we have two betas)
         for(int b = 0; b < scep->GetMultiplicity(); ++b) {
            btimestamp->Fill(scep->GetHit(b)->GetTime()/1e8);
            for(one = 0; one < (int) grif->GetMultiplicity(); ++one) {
               //Be careful about time ordering!!!! betas and gammas are not symmetric out of the DAQ
               //Fill the time diffrence spectra
               gbTimeDiff->Fill(grif->GetHit(one)->GetTime()-scep->GetHit(b)->GetTime());
               gbTimevsg->Fill(grif->GetHit(one)->GetTime()-scep->GetHit(b)->GetTime(),grif->GetHit(one)->GetEnergy());
               if((gbTlow <= grif->GetHit(one)->GetTime()-scep->GetHit(b)->GetTime()) && (grif->GetHit(one)->GetTime()-scep->GetHit(b)->GetTime() <= gbThigh)) {
                  gbEnergyvsbTime->Fill(scep->GetHit(b)->GetTime(),grif->GetHit(one)->GetEnergy());
               }
               if((gbTlow <= grif->GetHit(one)->GetTime()-scep->GetHit(b)->GetTime()) && (grif->GetHit(one)->GetTime()-scep->GetHit(b)->GetTime() <= gbThigh)) {
                  //Plots a gamma energy spectrum in coincidence with a beta
                  gbEnergyvsgTime->Fill(grif->GetGriffinHit(one)->GetTime()/1e8, grif->GetGriffinHit(one)->GetEnergy());
                  gammaSinglesB->Fill(grif->GetGriffinHit(one)->GetEnergy());
                  gammaSinglesB_hp->Fill(grif->GetGriffinHit(one)->GetEnergy(),scep->GetSceptarHit(b)->GetDetector());
                  grifscep_hp->Fill(grif->GetGriffinHit(one)->GetArrayNumber(),scep->GetSceptarHit(b)->GetDetector());
                  gammaSinglesBCyc->Fill(grif->GetCycleTimeInMilliSeconds(grif->GetHit(one)->GetTime())/100.,grif->GetHit(one)->GetEnergy());
                  //Now we want to loop over gamma rays if they are in coincidence.
                  if(grif->GetMultiplicity() > 1){
                     for(two = 0; two < (int) grif->GetMultiplicity(); ++two) {
                       if(two == one){ //If we are looking at the same gamma we don't want to call it a coincidence
                         continue;
                       }
                     
                        if(ggTlow <= TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) && TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) < ggThigh) { 
                           //If they are close enough in time, fill the gamma-gamma-beta matrix. This will be symmetric because we are doing a double loop over gammas
                           ggbmatrix->Fill(grif->GetGriffinHit(one)->GetEnergy(), grif->GetGriffinHit(two)->GetEnergy());
                        }
                        if(ggBGlow <= TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) && TMath::Abs(grif->GetGriffinHit(two)->GetTime()-grif->GetGriffinHit(one)->GetTime()) < ggBGhigh) { 
                           //If they are not close enough in time, fill the gamma-gamma-beta time-random matrix. This will be symmetric because we are doing a double loop over gammas
                           ggbmatrixt->Fill(grif->GetGriffinHit(one)->GetEnergy(), grif->GetGriffinHit(two)->GetEnergy());
                        }
                     }                   
                  }
               }
               if((gbBGlow <= grif->GetHit(one)->GetTime()-scep->GetHit(b)->GetTime()) && (grif->GetHit(one)->GetTime()-scep->GetHit(b)->GetTime() <= gbBGhigh)) {
                  gammaSinglesBt->Fill(grif->GetHit(one)->GetEnergy());                 
               }
            }
         }
      }
      if((entry%10000) == 1){
         printf("Completed %d of %d \r",entry,maxEntries);
      }

   }
   ggmatrixt->Scale(-ggBGScale);
   ggmatrixt->Add(ggmatrix);

   ggbmatrixt->Scale(-ggBGScale);
   ggbmatrixt->Add(ggbmatrix);

   gammaSinglesBt->Scale(-gbBGScale);
   gammaSinglesBt->Add(gammaSinglesB);

   list->Sort(); //Sorts the list alphabetically
   std::cout << "creating histograms done after " << w->RealTime() << " seconds" << std::endl;
   w->Continue();
   return list;

}











//This function gets run if running in compiled mode
#ifndef __CINT__ 
int main(int argc, char **argv) {
   if(argc != 4 && argc != 3 && argc != 2) {
      printf("try again (usage: %s <analysis tree file> <optional: output file> <max entries>).\n",argv[0]);
      return 0;
   }

   //We use a stopwatch so that we can watch progress
   TStopwatch w;
   w.Start();

   TFile* file = new TFile(argv[1]);

   if(file == NULL) {
      printf("Failed to open file '%s'!\n",argv[1]);
      return 1;
   }
   if(!file->IsOpen()) {
      printf("Failed to open file '%s'!\n",argv[1]);
      return 1;
   }
   printf("Sorting file:" DBLUE " %s" RESET_COLOR"\n",file->GetName());
   TGRSIRunInfo* runinfo  = (TGRSIRunInfo*)file->Get("TGRSIRunInfo");

   TTree* tree = (TTree*) file->Get("AnalysisTree");
   TChannel::ReadCalFromTree(tree);
   if(tree == NULL) {
      printf("Failed to find analysis tree in file '%s'!\n",argv[1]);
      return 1;
   }
   //Get the TGRSIRunInfo from the analysis Tree.
   
   TList *list;//We return a list because we fill a bunch of TH1's and shove them into this list.
   TFile * outfile;
   if(argc<3)
   {
      if(!runinfo){
         printf("Could not find run info, please provide output file name\n");
         return 0;
      }
      int runnumber = runinfo->RunNumber();
      int subrunnumber = runinfo->SubRunNumber();
      outfile = new TFile(Form("Leanmatrix%05d_%03d.root",runnumber,subrunnumber),"recreate");
   }
   else
   {
      outfile = new TFile(argv[2],"recreate");
   }

   std::cout << argv[0] << ": starting Analysis after " << w.RealTime() << " seconds" << std::endl;
   w.Continue();
   if(argc < 4) {
      list = exAnalysis(tree,0, &w);
   } else {
      int entries = atoi(argv[3]);
      std::cout<<"Limiting processing of analysis tree to "<<entries<<" entries!"<<std::endl;
      list = exAnalysis(tree, entries, &w);
   }
   printf("Writing to File: " DYELLOW "%s" RESET_COLOR"\n",outfile->GetName());
   list->Write();
   //Write the run info into the tree as well if there is run info in the Analysis Tree
   TGRSISortList *sortinfolist = new TGRSISortList;
   if(runinfo){
      TGRSISortInfo *info = new TGRSISortInfo(runinfo);
      sortinfolist->AddSortInfo(info);
      sortinfolist->Write("TGRSISortList",TObject::kSingleKey);
   }
   outfile->Close();

   std::cout << argv[0] << " done after " << w.RealTime() << " seconds" << std::endl << std::endl;


   return 0;
}

#endif
