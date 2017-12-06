//./muonHits test_shallowTrackCRUZET_2017VRNewClusteringTimingTestOnlyGlobalTracks.root  muHits 
//./muonHits test_shallowTrackCRUZET_2017VRNewClusteringTimingTestOnlyGlobalTracksiTimeTest.root  muHits 
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


       vector<float>*  muonsDTMuontrackDirection= 0;
       vector<float>*  muonsDTMuontrackSector= 0;
       vector<float>*  muonsDTMuontrackInOut= 0;
       vector<float>*  muonsDTMuontrackOutIn= 0;


       vector<float> submuonsDTMuontrackDirection;
       vector<float> submuonsDTMuontrackSector;
       vector<float> submuonsDTMuontrackInOut;
       vector<float> submuonsDTMuontrackOutIn;


       t1->SetBranchAddress("muonsDTMuontrackDirection",  &muonsDTMuontrackDirection );
       t1->SetBranchAddress("muonsDTMuontrackSector",  &muonsDTMuontrackSector );
       t1->SetBranchAddress("muonsDTMuontrackInOut",  &muonsDTMuontrackInOut );
       t1->SetBranchAddress("muonsDTMuontrackOutIn",  &muonsDTMuontrackOutIn );

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
           uint32_t up = muonsDTMuontrackDirection->size();
           for(uint32_t k=0; k<up;k++)
           {
               if( true )
               {
                   if( true ) //@MJ@ TODO trigger
                   {
                       submuonsDTMuontrackDirection.push_back(muonsDTMuontrackDirection->at(k));
                       submuonsDTMuontrackSector.push_back(muonsDTMuontrackSector->at(k));
                       submuonsDTMuontrackInOut.push_back(muonsDTMuontrackInOut->at(k));
                       submuonsDTMuontrackOutIn.push_back(muonsDTMuontrackOutIn->at(k));
                   }
               }
           }
           //perstrip

   }


   //}maybe




       TH2F* directionVSsector = new TH2F("directionVSsector","directionVSsector", 20, 0, 20, 4, -2, 2);         
       TH2F* InOutVSsector = new TH2F("InOutVSsector","InOutVSsector", 200, -100, 100, 20, 0, 20);         
       TH2F* OutInVSsector = new TH2F("InOutVSsector","InOutVSsector", 200, -100, 100, 20, 0, 20);         


       for(uint32_t t=0; t<submuonsDTMuontrackSector.size(); t++)
       {
           directionVSsector->Fill( submuonsDTMuontrackSector.at(t) , submuonsDTMuontrackDirection.at(t));
           if(submuonsDTMuontrackDirection.at(t)> 0)
               InOutVSsector->Fill( submuonsDTMuontrackInOut.at(t) , submuonsDTMuontrackSector.at(t)  );
           if(submuonsDTMuontrackDirection.at(t)< 0)
               OutInVSsector->Fill( submuonsDTMuontrackOutIn.at(t) , submuonsDTMuontrackSector.at(t)  );
       }
           cout<< "correlation " << directionVSsector->GetCorrelationFactor() << endl;

//////////////////////filling
       //directionVSsector->SetMarkerStyle(kFullCircle); 
       //directionVSsector->SetMarkerColor(kBlack); 

       //directionVSsector->SetTitle(""); 


       //directionVSsector->SetMaximum(1.5*  directionVSsector->GetMaximum());

       TCanvas c36("directionVSsector","directionVSsector");
       directionVSsector->Draw("colz"); 
       c36.SaveAs(("output/"+(string)dir+"/"+  "directionVSsector.root").c_str());
       c36.SaveAs(("output/"+(string)dir+"/"+  "directionVSsector.eps").c_str());

       TCanvas c37("InOutVSsector","InOutVSsector");
       InOutVSsector->Draw("colz"); 
       c37.SaveAs(("output/"+(string)dir+"/"+  "InOutVSsector.root").c_str());
       c37.SaveAs(("output/"+(string)dir+"/"+  "InOutVSsector.eps").c_str());

       TCanvas c37b("OutInVSsector","OutInVSsector");
       OutInVSsector->Draw("colz"); 
       c37b.SaveAs(("output/"+(string)dir+"/"+  "OutInVSsector.root").c_str());
       c37b.SaveAs(("output/"+(string)dir+"/"+  "OutInVSsector.eps").c_str());

    //cout << "in here 7"  << endl;

   return 0;
}
