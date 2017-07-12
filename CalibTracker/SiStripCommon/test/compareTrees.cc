//./compareTrees onTrackClustercosmicsData0T.root onTrackClustercosmicsMC0T.root TEC yes 0TdataMC > logsMean/log0TdataMcTEC
//./compareTrees onTrackClustercosmicsData38T.root onTrackClustercosmicsMC38T.root TEC yes 38TdataMC > logsMean/log38TdataMcTEC
//./compareTrees onTrackClustercosmicsData0T.root onTrackClustercosmicsData38T.root TIB yes 0Tto38Tdata > logsMean/log0Tto38TdataTIB
//./compareTrees test_shallowTrackClusterPUData.root test_shallowTrackClusterPUMCupgrade2017.root TID yes dataMCPU2017 > logsMean/logdataMCPU2017TID
//./compareTrees onTrackCluster284078NoPUData.root test_shallowTrackClusterNoPUMCupgrade2017.root TID yes 284078dataMC > logsMean/log284078dataMCTID
//./compareTrees test_shallowTrackAndClusterFullInfodataPU.root test_shallowTrackAndClusterFullInfoMCPUupgrade2017.root TID yes PUdataMCdataMCfullInfoUp2017 > logsMean/logPUdataMCdataMCfullInfoUp2017TID
//
//./compareTrees onTrackCluster284078NoPUData.root test_shallowTrackClusterNoPUMCrun2_mcTAG.root TID yes 284078data2016MC > logsMean/log284078data2016MCTID
//./compareTrees test_shallowTrackClusterPUData.root test_shallowTrackClusterPUMCrun2_mc.root TID yes data2016MCPU > logsMean/logdata2016MCPUTID
//./compareTrees onTrackCluster284078NoPUData.root test_shallowTrackClusterPUData.root TID yes noPUtoPUdata > logsMean/noPUtoPUdataTID
//
//./compareTrees test_shallowTrackAndClusterFullInfoNoPU2016Data.root test_shallowTrackAndClusterFullInfoNoPU2016MC.root TID yes noPU2016dataToMCFullInfo > logsMean/noPU2016dataToMCFullInfoTID
//./compareTrees noPU2016Data_v2_fullInfo.root noPU2016MC_v2_fullInfo.root TOB yes noPU2016dataToMCFullInfo > logsMean/noPU2016dataToMCFullInfoTOB
//./compareTrees test_shallowTrackAndClusterFullInfoNoPU2016Data.root test_shallowTrackAndClusterFullInfo.rootLA_0.022step3.root TID yes noPU2016dataToMCFullInfoFixedLA > logsMean/noPU2016dataToMCFullInfoFixedLATID
//./compareTrees test_shallowTrackAndClusterFullInfoNoPU2016Data.root test_shallowTrackAndClusterFullInfo.rootLA_0.001step3.root TID yes noPU2016dataToMCFullInfoFixedLAVerySmall > logsMean/noPU2016dataToMCFullInfoFixedLAVerySmallTID

#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream> 
#include <cmath>
#include <exception>
#include <ctime>

#include "TNtuple.h"
#include "TROOT.h"
#include "TString.h"
#include "TChain.h"
#include "TH1F.h"
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

using namespace std;

int main(int argc, char *argv[]){

    gStyle->SetOptStat(0);
gROOT->ForceStyle();

    if(argc != 6)
        throw std::runtime_error("Bad number of arguments!");
    
    TString file1 = argv[1];
    TString file2 = argv[2];
    TString subDet = argv[3];
    TString norm = argv[4];
    TString dir = argv[5];

    uint16_t sdId = 0;

    if(subDet == "TIB")
        sdId = 3;
    else if(subDet == "TID")
        sdId = 4;
    else if(subDet == "TOB")
        sdId = 5;
    else if(subDet == "TEC")
        sdId = 6;
    else
        throw std::runtime_error("Wrong partition entered");

    TFile* f1 = NULL;
    TTree* t1 = NULL;
    f1 = TFile::Open(file1); 
    TFile* f2 = NULL;
    TTree* t2 = NULL;
    f2 = TFile::Open(file2); 
    if(f1==NULL)
        throw std::runtime_error("File 1 address not set");
    t1 = dynamic_cast<TTree*>(f1->Get("testTree/tree"));
    if(t1==NULL)
        throw std::runtime_error("Tree 1 address not set");

    if(f2==NULL)
        throw std::runtime_error("File 2 address not set");
    t2 = dynamic_cast<TTree*>(f2->Get("testTree/tree"));
    if(t2==NULL)
        throw std::runtime_error("Tree 2 address not set");

    TObjArray *branchList = t1->GetListOfBranches();
    TObjArray *branchListT2 = t2->GetListOfBranches();

    if( branchList->GetEntries() != branchListT2->GetEntries() ) //@MJ@ TODO could be more sophiticated, liek comparing elements
        throw std::runtime_error("Trees have different content");

    vector<string> observables;
    for(int i = 0; i < branchList->GetEntries(); ++i) 
    { 
            observables.push_back(branchList->At(i)->GetName());
            //cout << branchList->At(i)->GetName() << endl; 
    }

   if(observables.size() == 0)
       throw std::runtime_error("The vector of observables not filled");

 
   vector< vector<float>* > var;
   vector< vector<float>* > varT2;
   vector< vector<float> > varTot;
   vector< vector<float> > varTotT2;
   var.resize(observables.size());
   varTot.resize(observables.size());
   varT2.resize(observables.size());
   varTotT2.resize(observables.size());

   uint32_t posSdId = 0;
   uint32_t puId = 0;

   cout << "data hist size " << var.size() << " mc hist size " << varT2.size() << endl;

   for(uint32_t v=0;v<var.size();v++)
   {
       if(observables.at(v) == "clustersubdetid")
           posSdId = v; 
       if(observables.at(v) == "PU")
          puId = v; 
       t1->SetBranchAddress(observables.at(v).c_str(), & (var.at(v)) );
       t2->SetBranchAddress(observables.at(v).c_str(), & (varT2.at(v)) );
   }
   //data always first
    

   uint32_t evCount=0;
   
   //cout << "in here a" << endl;
   Int_t nentries = (Int_t)t1->GetEntries();
   Int_t nentriesT2 = (Int_t)t2->GetEntries();

   ///fill variables from tree 1
   for (Int_t e=0; e<nentries; e++) 
   {
       t1->GetEntry(e);

      // float PU = var.at(puId).at(0);
       //if(PU>=5)
         //  continue;

       evCount++;

   //cout << "in here b" << endl;
       for(uint32_t j=0; j<var.size();j++)
       {
           //cout << "variable " << j << " address " << & (var.at(j)) << " size " << var.at(j)->size() << endl; 
          
           uint32_t up = var.at(posSdId)->size() < var.at(j)->size() ? var.at(posSdId)->size() :var.at(j)->size();
           for(uint32_t k=0; k<up;k++)
           {
               if(observables.at(j) == "clusternumber")
                   varTot.at(j).push_back(var.at(j)->at(k));

               if(var.at(posSdId)->at(k) == sdId)
               {
                   varTot.at(j).push_back(var.at(j)->at(k));
               }
   //cout << "in here d" << endl;
           }
       }

   }

   //fill variables from tree 2
   //for (Int_t e=0; e<58; e++) //test_shallowTrackClusterRunHMedPU.root  test_shallowTrackCluster900MC.root
   for (Int_t e=0; e<nentriesT2; e++) 
   {
       t2->GetEntry(e);

       for(uint32_t j=0; j<varT2.size();j++)
       {
   //cout << "in here e" << endl;
           
           uint32_t up = varT2.at(posSdId)->size() < varT2.at(j)->size() ? varT2.at(posSdId)->size() :varT2.at(j)->size();
           for(uint32_t k=0; k<up;k++)
           {
               if(observables.at(j) == "clusternumber")
                   varTotT2.at(j).push_back(varT2.at(j)->at(k));

               if(varT2.at(posSdId)->at(k) == sdId)
               {
                   varTotT2.at(j).push_back(varT2.at(j)->at(k));
               }
   //cout << "in here f" << endl;
           }
       }

   }

   cout << "in here g" << endl;
   vector<TH1F*> hists;
   vector<TH1F*> histsT2;
   //cout << "in here h" << endl;
   hists.resize(varTot.size());
   histsT2.resize(varTotT2.size());
   if(hists.size() != histsT2.size())
       throw std::runtime_error("The histogras do not have same size, something went terribly wrong!");

   for(uint32_t h=0; h<varTot.size(); h++) //@MJ@ TODO in case of empty vector do something!!!
   {
       if(varTot.at(h).size()==0) //@MJ@ TODO do it better
       {
           cout << "continuing 1 " << endl;
           continue;
       }
       if(varTotT2.at(h).size()==0) //@MJ@ TODO do it better
       {
           cout << "continuing 2 " << endl;
           continue;
       }

       sort(varTot.at(h).begin(), varTot.at(h).end());
       sort(varTotT2.at(h).begin(), varTotT2.at(h).end());

    //cout << "in here 1"  << endl;
       //float start = varTot.at(h).at(0) < varTotT2.at(h).at(0) ? varTot.at(h).at(0) : varTotT2.at(h).at(0);
       float start = varTot.at(h).at(0);
    //cout << "in here 1.1"  << endl;
       float size1 = varTot.at(h).size() -1;
    //cout << "in here 1.2"  << endl;
       float size2 = varTotT2.at(h).size() -1;
    //cout << "in here 1.3"  << endl;
       float end = varTot.at(h).at(size1);//
       //float end = varTot.at(h).at(size1) > varTotT2.at(h).at(size2) ? varTot.at(h).at(size1) : varTotT2.at(h).at(size2);
    //cout << "in here 1.4"  << endl;
       hists.at(h) = new TH1F(observables.at(h).c_str(), observables.at(h).c_str(), 1000, start, end );
       //hists.at(h) = new TH1F(observables.at(h).c_str(), observables.at(h).c_str(), 10, 0, 10 );
    //cout << "in here 1.5"  << endl;
       histsT2.at(h) = new TH1F((observables.at(h)+"T2").c_str(), (observables.at(h)+"T2").c_str(), 1000, start, end );
       //histsT2.at(h) = new TH1F((observables.at(h)+"T2").c_str(), (observables.at(h)+"T2").c_str(), 10, 0, 10 );

    cout << "in here 2"  << endl;
       for(uint32_t m = 0; m<varTot.at(h).size(); m++)
       {
           hists.at(h)->Fill(varTot.at(h).at(m));
       }
    cout << "in here 3"  << endl;
       for(uint32_t m = 0; m<varTotT2.at(h).size(); m++)
       {
           histsT2.at(h)->Fill(varTotT2.at(h).at(m));
       }
    cout << "in here 4"  << endl;
       TCanvas c(observables.at(h).c_str(),observables.at(h).c_str());
       TPad pad1(("pad1"+observables.at(h)).c_str(), "The pad 80% of the height",0.0,0.2,1.0,1.0);
       TPad pad2(("pad2"+observables.at(h)).c_str(), "The pad 20% of the height",0.0,0.0,1.0,0.2);
       pad1.Draw();
       pad2.Draw();

       pad1.cd();
 
    //cout << "in here 5"  << endl;

       //hists.at(h)->SetLineColor(kBlack);
       hists.at(h)->SetMarkerStyle(kFullCircle);
       hists.at(h)->SetMarkerColor(kBlack);
       hists.at(h)->GetXaxis()->SetTitle(observables.at(h).c_str());
       hists.at(h)->GetYaxis()->SetTitle("entries");
       float maxm = hists.at(h)->GetMaximum();
       hists.at(h)->SetMaximum(1.3*maxm);
       hists.at(h)->SetTitle("");
       histsT2.at(h)->SetTitle("");
       //hists.at(h)->Draw("*H");
       hists.at(h)->Draw("P");

       double n = hists.at(h)->GetEntries();
       //cout << "norm 1 " << hists.at(h)->GetNormFactor() << " netries " << n << end;
       double n2 = histsT2.at(h)->GetEntries();
       //double n = 80654/5948; 
       histsT2.at(h)->SetLineColor(kRed);
       if(norm == "yes")
       {
           //histsT2.at(h)->Scale(n);
           //histsT2.at(h)->Draw("same");
           //histsT2.at(h)->DrawNormalized("same", n);
           histsT2.at(h)->DrawNormalized("same", n);
       }
       else
           histsT2.at(h)->Draw("same");
  
       cout << "hist " << observables.at(h) << " mean 1 " << hists.at(h)->GetMean() << "  mean2 " <<  histsT2.at(h)->GetMean()<< " mean 1 err " << hists.at(h)->GetMeanError() << "  mean2 err " <<  histsT2.at(h)->GetMeanError() << endl;
    //cout << "in here 6"  << endl;
       pad2.cd();
       TH1F* histD = (TH1F*) histsT2.at(h)->Clone();
       histD->Scale(n/n2);
       histD->Divide(hists.at(h));
       histD->SetLineColor(kBlue);
       histD->SetMaximum(2);
       histD->SetMinimum(0);
       histD->Draw("");
       cout << "hist " << h << " drawn" << endl; 

       c.Update();
       c.SaveAs(("output/"+(string)dir+"/"+observables.at(h) + (string)subDet + ".root").c_str());
       //c.SaveAs(("output/"+(string)dir+"/"+observables.at(h) + (string)subDet + ".eps").c_str());
       delete histD;
   }

    //cout << "in here 7"  << endl;

   return 0;
}
