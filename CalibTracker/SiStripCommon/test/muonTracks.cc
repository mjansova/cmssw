//./muonTracks test_shallowTrackCRUZET_2017VRNewClusteringTimingMuonTimingTestsZOffsetPlusCTvzPlusEta.root  muonTracksInfo 
//

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream> 
#include <cmath>
#include <exception>
#include <ctime>
#include <cmath>

#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TCanvas.h"
#include "TFitResult.h"
#include "TPaveText.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TObject.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TStyle.h"
#include "TROOT.h"
#include "THashList.h"
#include "TApplication.h"
#include "TGraph.h"
#include "TMath.h"
#include "TLegend.h"
#include "TProfile.h"

using namespace std;

int main(int argc, char *argv[]){

 gStyle->SetOptStat(0);
 gROOT->ForceStyle();
 
 TH1::SetDefaultSumw2();
    if(argc != 3)
        throw std::runtime_error("Bad number of arguments!");
    
    TString file1 = argv[1];
    TString dir = argv[2];

    float narrowness = 0.5;


/*    if(subDet == "TIB")
        sdId = 3;
    else if(subDet == "TID")
        sdId = 4;
    else if(subDet == "TOB")
        sdId = 5;
    else if(subDet == "TEC")
        sdId = 6;
    else
        throw std::runtime_error("Wrong partition entered");
*/
    TFile* f1 = NULL;
    TTree* t1 = NULL;
    f1 = TFile::Open(file1); 
    if(f1==NULL)
        throw std::runtime_error("File 1 address not set");
    t1 = dynamic_cast<TTree*>(f1->Get("testTree/tree"));
    if(t1==NULL)
        throw std::runtime_error("Tree 1 address not set");


       vector<float>*  muonsDTMuontrackInOutTop = 0;
       vector<float>*  muonsDTMuontrackInOutErrTop = 0;
       vector<float>*  muonsDTMuontrackOutInTop = 0;
       vector<float>*  muonsDTMuontrackOutInErrTop = 0;
       vector<float>*  innerXtop = 0;
       vector<float>*  innerYtop = 0;
       vector<float>*  innerVZtop = 0;
       vector<float>*  outerXtop = 0;
       vector<float>*  outerYtop = 0;
       vector<float>*  outerZtop = 0;
       vector<float>*  outerEtatop = 0;

       vector<float>  submuonsDTMuontrackInOutTop;
       vector<float>  submuonsDTMuontrackInOutErrTop;
       vector<float>  submuonsDTMuontrackOutInTop ;
       vector<float>  submuonsDTMuontrackOutInErrTop;
       vector<float>  subinnerXtop;
       vector<float>  subinnerYtop;
       vector<float>  subinnerVZtop ;
       vector<float>  subouterXtop;
       vector<float>  subouterYtop;
       vector<float>  subouterZtop;
       vector<float>  subouterEtatop;




       t1->SetBranchAddress("muonsDTMuontrackInOutTop",  &muonsDTMuontrackInOutTop );
       t1->SetBranchAddress("muonsDTMuontrackInOutErrTop",  &muonsDTMuontrackInOutErrTop );
       t1->SetBranchAddress("muonsDTMuontrackOutInTop",  &muonsDTMuontrackOutInTop );
       t1->SetBranchAddress("muonsDTMuontrackOutInErrTop",  &muonsDTMuontrackOutInErrTop );
       t1->SetBranchAddress("innerXtop",  &innerXtop );
       t1->SetBranchAddress("innerYtop",  &innerYtop );
       t1->SetBranchAddress("innerVZtop",  &innerVZtop );
       t1->SetBranchAddress("outerXtop",  &outerXtop );
       t1->SetBranchAddress("outerYtop",  &outerYtop );
       t1->SetBranchAddress("outerZtop",  &outerZtop );
       t1->SetBranchAddress("outerEtatop",  &outerEtatop );

   //data always first
    

   uint32_t evCount=0;
   
   //cout << "in here a" << endl;
   Int_t nentries = (Int_t)t1->GetEntries();

   //cout << "in here b" << endl;
   ///fill variables from tree 1
   for (Int_t e=0; e<nentries; e++) 
   {
       t1->GetEntry(e);
          
           //per cluster
           uint32_t up = muonsDTMuontrackInOutTop->size();
           for(uint32_t k=0; k<up;k++)
           {
               if( true )
               {
                   if( true ) //@MJ@ TODO trigger
                   {
                       submuonsDTMuontrackInOutTop.push_back(muonsDTMuontrackInOutTop->at(k));
                       submuonsDTMuontrackInOutErrTop.push_back(muonsDTMuontrackInOutErrTop->at(k));
                       submuonsDTMuontrackOutInTop.push_back(muonsDTMuontrackOutInTop->at(k));
                       submuonsDTMuontrackOutInErrTop.push_back(muonsDTMuontrackOutInErrTop->at(k));
                       subinnerXtop.push_back(innerXtop->at(k));
                       subinnerYtop.push_back(innerYtop->at(k));
                       subinnerVZtop.push_back(innerVZtop->at(k));
                       subouterXtop.push_back(outerXtop->at(k));
                       subouterYtop.push_back(outerYtop->at(k));
                       subouterZtop.push_back(outerZtop->at(k));
                       subouterEtatop.push_back(outerEtatop->at(k));
                   }
               }
           }
           //perstrip

   }


   //}maybe




       TH2F* OutInVSinnerZ = new TH2F("OutInVSinnerVZ","OutInVSinnerVZ", 100, -100, 100, 200, -300, 300);         
       TH2F* OutInErrVSinnerZ = new TH2F("OutInErrVSinnerVZ","OutInErrVSinnerVZ", 40, 0, 40, 200, -300, 300);         

       TH1F*  OutInCenter=  new TH1F("OutInCenter", "OutInCenter" , 200, -100, 100) ;
       TH1F*  OutInEdge=  new TH1F("OutInEdge", "OutInEdge" , 200, -100, 100) ;
       TH1F*  OutInErrCenter=  new TH1F("OutInErrCenter", "OutInErrCenter" , 40, 0, 20) ;
       TH1F*  OutInErrEdge=  new TH1F("OutInErrEdge", "OutErrInEdge" , 40, 0, 20) ;


       for(uint32_t t=0; t<submuonsDTMuontrackOutInTop.size(); t++)
       {
           OutInVSinnerZ->Fill(  submuonsDTMuontrackOutInTop.at(t) , subinnerVZtop.at(t));
           OutInErrVSinnerZ->Fill(  submuonsDTMuontrackOutInErrTop.at(t) , subinnerVZtop.at(t));
           if(subinnerVZtop.at(t) > -50 && subinnerVZtop.at(t) < 50)
           {
               OutInCenter->Fill(submuonsDTMuontrackOutInTop.at(t)) ;
               OutInErrCenter->Fill(submuonsDTMuontrackOutInErrTop.at(t)) ;
           }
           if(subinnerVZtop.at(t) < -50 || subinnerVZtop.at(t) > 50)
           {
               OutInEdge->Fill(submuonsDTMuontrackOutInTop.at(t)) ;
               OutInErrEdge->Fill(submuonsDTMuontrackOutInErrTop.at(t)) ;
           }
       }
       cout<< "correlation OutInVSinnerZ" << OutInVSinnerZ->GetCorrelationFactor() << endl;
       cout<< "correlation OutInErrVSinnerZ" << OutInErrVSinnerZ->GetCorrelationFactor() << endl;

//////////////////////filling
       //directionVSsector->SetMarkerStyle(kFullCircle); 
       //directionVSsector->SetMarkerColor(kBlack); 

       //directionVSsector->SetTitle(""); 


       //directionVSsector->SetMaximum(1.5*  directionVSsector->GetMaximum());


       TCanvas c37("OutInVSinnerVZ","OutInVSinnerVZ");
       OutInVSinnerZ->Draw("colz"); 
       c37.SaveAs(("output/"+(string)dir+"/"+  "OutInVSinnerVZ.root").c_str());
       c37.SaveAs(("output/"+(string)dir+"/"+  "OutInVSinnerVZ.eps").c_str());

       TCanvas c38("OutInErrVSinnerVZ","OutInErrVSinnerVZ");
       OutInErrVSinnerZ->Draw("colz"); 
       c38.SaveAs(("output/"+(string)dir+"/"+  "OutInErrVSinnerVZ.root").c_str());
       c38.SaveAs(("output/"+(string)dir+"/"+  "OutInErrVSinnerVZ.eps").c_str());

       TCanvas c39("OutInVsCenterEdge","OutInVsCenterEdge");
       OutInCenter->SetLineColor(kBlue);
       OutInEdge->SetTitle("");
       OutInEdge->GetXaxis()->SetTitle("OutIn time");
       OutInEdge->SetLineColor(kRed);
       OutInEdge->DrawNormalized("hist");
       OutInCenter->DrawNormalized("same hist");
       c39.SaveAs(("output/"+(string)dir+"/"+  "OutInVsCenterEdge.root").c_str());
       c39.SaveAs(("output/"+(string)dir+"/"+  "OutInVsCenterEdge.eps").c_str());

       TCanvas c40("OutInErrVsCenterEdge","OutInErrVsCenterEdge");
       OutInErrCenter->SetLineColor(kBlue);
       OutInErrEdge->SetLineColor(kRed);
       OutInErrCenter->DrawNormalized("hist");
       OutInErrEdge->DrawNormalized("same hist");
       c40.SaveAs(("output/"+(string)dir+"/"+  "OutInErrVsCenterEdge.root").c_str());
       c40.SaveAs(("output/"+(string)dir+"/"+  "OutInErrVsCenterEdge.eps").c_str());


    //cout << "in here 7"  << endl;

   return 0;
}
