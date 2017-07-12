//./compareData test_shallowCluster910RunHVtx.root test_shallowCluster910RunEVtx.root TOB no
//./compareData test_shallowCluster910RunHVtx.root test_shallowCluster910RunBVtx.root TOB no
//./compareData test_shallowCluster910RunBVtx.root test_shallowCluster910RunEVtx.root TOB no


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

    if(argc != 5)
        throw std::runtime_error("Bad number of arguments!");
    
    TString file1 = argv[1];
    TString file2 = argv[2];
    TString subDet = argv[3];
    TString norm = argv[4];

    uint16_t sdId = 0;
    uint16_t layerId = 0;

    if(subDet == "TIB")
        sdId = 3;
    else if(subDet == "TID")
        sdId = 4;
    else if(subDet == "TOB")
    {
        sdId = 5;
        layerId = 0;
    }
    else if(subDet == "TOBl1")
    {
        sdId = 5;
        layerId = 1;
    }
    else if(subDet == "TOBlm3")
    {
        sdId = 5;
        layerId = 3;
    }
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
   //vector< vector<float> > wvarTot;
   vector< vector<float> > wvarTotT2;
   var.resize(observables.size());
   varTot.resize(observables.size());
   varT2.resize(observables.size());
   varTotT2.resize(observables.size());
   //wvarTot.resize(observables.size());
   wvarTotT2.resize(observables.size());

   uint32_t posSdId = 0;
   uint32_t posLayerId = 0;
   uint32_t puId = 0;

   cout << "data hist size " << var.size() << " mc hist size " << varT2.size() << endl;

   for(uint32_t v=0;v<var.size();v++)
   {
       if(observables.at(v) == "clustersubdetid")
           posSdId = v; 
       if(observables.at(v) == "clusterlayerwheel")
           posLayerId = v; 
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
   //Int_t nentries = 500;
   //Int_t nentriesT2 = 500;

   
   TH1F* hPU1 = new TH1F("PU1", "PU1", 100, 0, 100 );
   TH1F* hPU2 = new TH1F("PU2", "PU2", 100, 0, 100 );
   TH1F* hPU;

   ///fill variables from tree 1
   for (Int_t e=0; e<nentries; e++) 
   {
       t1->GetEntry(e);
       hPU1->Fill(var.at(puId)->at(0));

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
               if(observables.at(j) == "clusternumber") //not correct for layers
               {
                   varTot.at(j).push_back(var.at(j)->at(k));
                   //wvarTot.at(j).push_back(var.at(puId)->at(0));
               }

               if(var.at(posSdId)->at(k) == sdId)
               {

                   if((sdId == 5 && layerId ==0) || sdId!=5 )
                   varTot.at(j).push_back(var.at(j)->at(k));

                   
                   if(sdId == 5 && layerId ==1)
                   {
                       if(var.at(posLayerId)->at(k) == layerId)
                           varTot.at(j).push_back(var.at(j)->at(k));
                   }
                   if(sdId == 5 && layerId>1)
                   {
                       if(var.at(posLayerId)->at(k) >= layerId)
                           varTot.at(j).push_back(var.at(j)->at(k));
                   }
                   //wvarTot.at(j).push_back(var.at(puId)->at(0));
               }
   //cout << "in here d" << endl;
           }
       }

   }

   //fill variables from tree 2
   //for (Int_t e=0; e<81; e++) 
   for (Int_t e=0; e<nentriesT2; e++) 
   {
       t2->GetEntry(e);
       hPU2->Fill(varT2.at(puId)->at(0));

       for(uint32_t j=0; j<varT2.size();j++)
       {
   //cout << "in here e" << endl;
           
           uint32_t up = varT2.at(posSdId)->size() < varT2.at(j)->size() ? varT2.at(posSdId)->size() :varT2.at(j)->size();
           for(uint32_t k=0; k<up;k++)
           {
               if(observables.at(j) == "clusternumber")
               {
                   varTotT2.at(j).push_back(varT2.at(j)->at(k));
                   wvarTotT2.at(j).push_back(varT2.at(puId)->at(0));
               }

               if(varT2.at(posSdId)->at(k) == sdId)
               {
                   if((sdId == 5 && layerId ==0) || sdId!=5)
                   {
                       varTotT2.at(j).push_back(varT2.at(j)->at(k));
                       wvarTotT2.at(j).push_back(varT2.at(puId)->at(0));
                   }

                   
                   if(sdId == 5 && layerId ==1)
                   {
                       if(varT2.at(posLayerId)->at(k) == layerId)
                       {
                       varTotT2.at(j).push_back(varT2.at(j)->at(k));
                       wvarTotT2.at(j).push_back(varT2.at(puId)->at(0));
                       }
                   }
                   if(sdId == 5 && layerId>1)
                   {
                       if(varT2.at(posLayerId)->at(k) >= layerId)
                       {
                       varTotT2.at(j).push_back(varT2.at(j)->at(k));
                       wvarTotT2.at(j).push_back(varT2.at(puId)->at(0));
                       }
                   }
               }
   //cout << "in here f" << endl;
           }
       }

   }

   hPU = (TH1F*) hPU1->Clone();
   hPU->Divide(hPU2);

   cout << "in here g" << endl;
   vector<TH1F*> hists;
   vector<TH1F*> histsT2;
   //ncout << "in here h" << endl;
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

       vector<float> varTots = varTot.at(h);
       vector<float> varTotT2s = varTotT2.at(h);

       sort(varTots.begin(), varTots.end());
       sort(varTotT2s.begin(), varTotT2s.end());

    //cout << "in here 1"  << endl;
       float start = varTots.at(0) < varTotT2s.at(0) ? varTots.at(0) : varTotT2s.at(0);
    //cout << "in here 1.1"  << endl;
       float size1 = varTots.size() -1;
    //cout << "in here 1.2"  << endl;
       float size2 = varTotT2s.size() -1;
    //cout << "in here 1.3"  << endl;
       float end = varTots.at(size1) > varTotT2s.at(size2) ? varTots.at(size1) : varTotT2s.at(size2);
    //cout << "in here 1.4"  << endl;
       hists.at(h) = new TH1F(observables.at(h).c_str(), observables.at(h).c_str(), 500, start, end );
       //hists.at(h) = new TH1F(observables.at(h).c_str(), observables.at(h).c_str(), 100, 0, 1000 );
    //cout << "in here 1.5"  << endl;
       histsT2.at(h) = new TH1F((observables.at(h)+"T2").c_str(), (observables.at(h)+"T2").c_str(), 500, start, end );
       //histsT2.at(h) = new TH1F((observables.at(h)+"T2").c_str(), (observables.at(h)+"T2").c_str(), 100, 0, 1000 );

    //cout << "in here 2"  << endl;
       for(uint32_t m = 0; m<varTot.at(h).size(); m++)
       {
           hists.at(h)->Fill(varTot.at(h).at(m));
       }
    //cout << "in here 3"  << endl;
       for(uint32_t m = 0; m<varTotT2.at(h).size(); m++)
       {
           uint32_t b = hPU->FindBin(wvarTotT2.at(h).at(m));
           float w = hPU->GetBinContent(b);
    //cout << "in here 4"  << endl;
           //cout << " PU2 " << wvarTotT2.at(h).at(m) << " weight "<< w << endl;
           histsT2.at(h)->Fill(varTotT2.at(h).at(m), w);
       }
    //cout << "in here 4"  << endl;
       TCanvas c(observables.at(h).c_str(),observables.at(h).c_str());
    //cout << "in here 5"  << endl;
       TH1F* histD = (TH1F*) hists.at(h)->Clone();

       //hists.at(h)->SetLineColor(kBlack);
       hists.at(h)->SetMarkerStyle(kFullCircle);
       hists.at(h)->SetMarkerColor(kBlack);
       hists.at(h)->GetXaxis()->SetTitle(observables.at(h).c_str());
       float maxm = hists.at(h)->GetMaximum();
       hists.at(h)->SetMaximum(4*maxm);
       //hists.at(h)->Draw("*H");
       hists.at(h)->Draw("P");

       double n = hists.at(h)->GetEntries();
       histsT2.at(h)->SetMarkerStyle(kFullCircle);
       histsT2.at(h)->SetMarkerColor(kRed);
       if(norm == "yes")
           histsT2.at(h)->DrawNormalized("same P", n);
       else
           histsT2.at(h)->Draw("same P");
   
       cout << "hist " << observables.at(h) << " mean 1 " << hists.at(h)->GetMean() << " mean2 " <<  histsT2.at(h)->GetMean() << endl;
    //cout << "in here 6"  << endl;
       histD->Divide(histsT2.at(h));
       histD->Draw("same");
       cout << "hist " << h << " drawn" << endl; 
       c.SaveAs(("output/"+observables.at(h) + (string)subDet + ".root").c_str());
       delete histD;
   }

    TCanvas ca("PUd","PUd");
    hPU->Draw();
    ca.SaveAs("hPU.root");
    TCanvas ca1("PUd1","PUd1");
    hPU1->Draw();
    ca1.SaveAs("hPU1.root");
    TCanvas ca2("PUd2","PUd2");
    hPU2->Draw();
    ca2.SaveAs("hPU2.root");
    cout << "in here 7"  << endl;

   return 0;
}
