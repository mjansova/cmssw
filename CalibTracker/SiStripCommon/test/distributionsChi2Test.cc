//./distributionsChi2Test noPU2016Data_v2_fullInfo.root test_shallowTrackAndClusterFullInfo.rootCT3_0.8401step3.root 0.07 0.01 TOB 0 5 yes CTtuning 
//
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

 //gStyle->SetOptStat(0);
 //gROOT->ForceStyle();
 TH1::SetDefaultSumw2();
    if(argc != 10)
    {
        cout << "argc " << argc<< endl;
        for(uint32_t a= 0; a<argc; a++)
        {
            cout << a << " argument value " << argv[a] << endl;
        }
        throw std::runtime_error("Bad number of arguments!");
    }
    
    TString file1 = argv[1];
    TString file2 = argv[2];
    TString spar2 = argv[3];
    TString spar3 = argv[4];
    TString subDet = argv[5];
    TString slayerCutDown = argv[6];
    TString slayerCutUp = argv[7];
    TString norm = argv[8];
    TString dir = argv[9];

    float par2 = spar2.Atof();
    float par3 = spar3.Atof();
    float par1 = 1- (2*par2) -(2*par3);
    TString spar1 = TString(to_string(par1));

    int layerCutDown = slayerCutDown.Atoi();
    int layerCutUp = slayerCutUp.Atoi();

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

       vector<float>* partition = 0;
       vector<float>* partition2 = 0;
       vector<float>* clustercharge = 0;
       vector<float>* clustercharge2 =0;
       vector<float>* clusterwidth = 0;
       vector<float>* clusterwidth2 = 0;
       vector<float>* clusterlayerwheel = 0;
       vector<float>* clusterlayerwheel2 = 0;

       vector<float> subpartition;
       vector<float> subpartition2;
       vector<float> subclustercharge;
       vector<float> subclustercharge2;
       vector<float> subclusterwidth;
       vector<float> subclusterwidth2;
       vector<float> subclusterlayerwheel;
       vector<float> subclusterlayerwheel2;

       t1->SetBranchAddress("clustersubdetid",  &partition );
       t1->SetBranchAddress("clustercharge",  &clustercharge );
       t1->SetBranchAddress("clusterwidth",  &clusterwidth );
       t1->SetBranchAddress("clusterlayerwheel",  &clusterlayerwheel );

       t2->SetBranchAddress("clustersubdetid",  &partition2 );
       t2->SetBranchAddress("clustercharge",  &clustercharge2 );
       t2->SetBranchAddress("clusterwidth",  &clusterwidth2 );
       t2->SetBranchAddress("clusterlayerwheel",  &clusterlayerwheel2 );

   //data always first
    

   uint32_t evCount=0;
   
   cout << "in here a" << endl;
   Int_t nentries = (Int_t)t1->GetEntries();
   Int_t nentriesT2 = (Int_t)t2->GetEntries();
   cout << "entries " << nentries << "enstries 2" << nentriesT2 << endl;

   cout << "in here b" << endl;
   ///fill variables from tree 1
   for (Int_t e=0; e<nentries; e++) 
   {
       t1->GetEntry(e);
          
           //per cluster
           uint32_t up = partition->size();
           for(uint32_t k=0; k<up;k++)
           {
               if( partition->at(k) == sdId )
               {
                   subpartition.push_back(partition->at(k));
                   subclustercharge.push_back(clustercharge->at(k));
                   subclusterwidth.push_back(clusterwidth->at(k));
                   subclusterlayerwheel.push_back(clusterlayerwheel->at(k));
               }
           }

   }

   //fill variables from tree 2
   //for (Int_t e=0; e<58; e++) //test_shallowTrackClusterRunHMedPU.root  test_shallowTrackCluster900MC.root
   for (Int_t e=0; e<nentriesT2; e++) 
   {
       t2->GetEntry(e);

          //per cluster 
           uint32_t up = partition2->size();
           for(uint32_t k=0; k<up;k++)
           {

               if( partition2->at(k) == sdId)
               {
                   subpartition2.push_back(partition2->at(k));
                   subclustercharge2.push_back(clustercharge2->at(k));
                   subclusterwidth2.push_back(clusterwidth2->at(k));
                   subclusterlayerwheel2.push_back(clusterlayerwheel2->at(k));
               }
           }
   }



       TH1F* chargeForAllWidthsData = new TH1F("chargeForAllWidths", "chargeForAllWidths" , 100, 0, 1000 );
       TH1F* chargeForAllWidthsMC = new TH1F("chargeForAllWidthsMC", "chargeForAllWidthsMC" , 100, 0, 1000 );
       TH1F* widthAllTo10Data = new TH1F("widthAllTo10Data", "widthAllTo10Data" , 9, 1, 10  );
       TH1F* widthAllTo10MC = new TH1F("widthAllTo10MC", "widthAllTo10MC" , 9, 1, 10 );
       TH1F* widthAllTo15Data = new TH1F("widthAllTo15Data", "widthAllTo15Data" , 14, 1, 15  );
       TH1F* widthAllTo15MC = new TH1F("widthAllTo15MC", "widthAllTo15MC" , 14, 1, 15 );
       TH1F* widthAllTo20Data = new TH1F("widthAllTo20Data", "widthAllTo20Data" , 19, 1, 20  );
       TH1F* widthAllTo20MC = new TH1F("widthAllTo20MC", "widthAllTo20MC" , 19, 1, 20 );
   
       cout << "in here 2"  << endl;
       for(uint32_t m = 0; m<subclusterwidth.size(); m++)
       {
           if(layerCutDown< subclusterlayerwheel.at(m) && subclusterlayerwheel.at(m) < layerCutUp)
           {
               widthAllTo10Data->Fill(subclusterwidth.at(m));
               widthAllTo15Data->Fill(subclusterwidth.at(m));
               widthAllTo20Data->Fill(subclusterwidth.at(m));
           }
         
       }
   
    cout << "in here 3"  << endl;
       for(uint32_t m = 0; m<subclusterwidth2.size(); m++)
       {
           if(layerCutDown< subclusterlayerwheel2.at(m) &&  subclusterlayerwheel2.at(m)< layerCutUp)
           {
               widthAllTo10MC->Fill(subclusterwidth2.at(m));
               widthAllTo15MC->Fill(subclusterwidth2.at(m));
               widthAllTo20MC->Fill(subclusterwidth2.at(m));
           }
       }

       for(uint32_t c=0; c<widthAllTo10MC->GetNbinsX(); c++)
       {
           cout << "MC bin " << c+1 << " n events " << widthAllTo10MC->GetBinContent(c+1) << endl;
       }

       TFile *f_chi2 = new TFile (("output/"+(string)dir+"/"+(string)dir+"chi2.root").c_str(), "UPDATE");
       TH2D* chi22D = NULL;
       chi22D = dynamic_cast<TH2D*>(f_chi2->Get("chi22D"));
       if(chi22D == NULL)
       {
           cout << "new chi2 histogram created " << endl;
           chi22D = new TH2D("chi22D", "chi22D", 1000, 0.01, 0.2, 1000, 0, 0.07);
       }

       
       TAxis* xaxis = chi22D->GetXaxis();
       TAxis *yaxis = chi22D->GetYaxis();
       int binx = xaxis->FindBin(par2);
       int biny = yaxis->FindBin(par3);

       chi22D->SetBinContent(binx, biny, widthAllTo10Data->Chi2Test(widthAllTo10MC, "CHI2"));
       chi22D->Write(chi22D->GetName(),TObject::kOverwrite);
       f_chi2->Close();

       TFile *f_Kolmogorov = new TFile (("output/"+(string)dir+"/"+(string)dir+"Kolmogorov.root").c_str(), "UPDATE");
       TH2D* kol2D = NULL;
       kol2D = dynamic_cast<TH2D*>(f_Kolmogorov->Get("kol2D"));
       if(kol2D == NULL)
       {
           cout << "new kol histogram created " << endl;
           kol2D = new TH2D("kol2D", "kol2D", 1000, 0.01, 0.2, 1000, 0, 0.07);
       }
       kol2D->SetBinContent(binx, biny, widthAllTo10Data->KolmogorovTest(widthAllTo10MC));
       kol2D->Write(kol2D->GetName(),TObject::kOverwrite);
       f_Kolmogorov->Close();

/*       cout << "chi2 widthAllTo10 " << widthAllTo10Data->Chi2Test(widthAllTo10MC, "CHI2") << endl;;
       cout << "Kolmogorov widthAllTo10 " << widthAllTo10Data->KolmogorovTest(widthAllTo10MC) << endl;;
       cout << "chi2 widthAllTo15 " << widthAllTo15Data->Chi2Test(widthAllTo15MC, "CHI2") << endl;;
       cout << "Kolmogorov widthAllTo15 " << widthAllTo15Data->KolmogorovTest(widthAllTo15MC) << endl;;
       cout << "chi2 widthAllTo20 " << widthAllTo20Data->Chi2Test(widthAllTo20MC, "CHI2") << endl;;
       cout << "Kolmogorov widthAllTo20 " << widthAllTo20Data->KolmogorovTest(widthAllTo20MC) << endl;;
*/

//////////////////////filling
       //chargeForAllWidthsData->SetMarkerStyle(kFullCircle);
       //chargeForAllWidthsMC->SetMarkerColor(kBlack);
       widthAllTo10Data->SetMarkerStyle(kFullCircle);
       widthAllTo10Data->SetMarkerColor(kBlack);


       //chargeForAllWidthsMC->SetTitle("");
       widthAllTo10Data->SetTitle("");

       //chargeForAllWidthsData->GetXaxis()->SetTitle("cluster charge");
       widthAllTo10Data->GetXaxis()->SetTitle("cluster width");

       //chargeForAllWidthsData->SetMaximum(1.5*  (chargeForAllWidthsData->GetMaximum()));
       widthAllTo10Data->SetMaximum(1.5*  (widthAllTo10Data->GetMaximum()));

       //chargeForAllWidthsMC->SetLineColor(kRed);
       widthAllTo10MC->SetLineColor(kRed);



       //can
       /*TCanvas c0("chargeForAllWidths","chargeForAllWidths");
       chargeForAllWidthsData->DrawNormalized("P");
       chargeForAllWidthsMC->DrawNormalized("same"); 
       c0.SaveAs(("output/"+(string)dir+"/"+ (string)slayerCutDown+ "lt" +(string)subDet +"mt"+(string)slayerCutUp +"_CT_"+ (string)spar1 + "_" + (string)spar2 + "_" +(string)spar3 + "_"+  "chargeForAllWidths.root").c_str());
       c0.SaveAs(("output/"+(string)dir+"/"+ (string)slayerCutDown+ "lt" +(string)subDet +"mt"+(string)slayerCutUp +"_CT_"+ (string)spar1 + "_" + (string)spar2 + "_" +(string)spar3 + "_"+  "chargeForAllWidths.eps").c_str());
*/

       TCanvas c1("chargeForAllWidths","chargeForAllWidths");
       TPad pad1("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0);
       TPad pad2("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.2);
       pad1.Draw();
       pad2.Draw();
       
       pad1.cd();
       widthAllTo10Data->DrawNormalized("P");
       widthAllTo10MC->DrawNormalized("same"); 
       pad2.cd();
       TH1F* widthAllTo10Div = (TH1F*) widthAllTo10Data->Clone();
       widthAllTo10Div->Divide(widthAllTo10MC);
       widthAllTo10Div->SetMaximum(2);
       widthAllTo10Div->SetMinimum(0);
       widthAllTo10Div->Draw("hist");
       c1.SaveAs(("output/"+(string)dir+"/"+ (string)slayerCutDown+ "lt" +(string)subDet +"mt"+(string)slayerCutUp +"_CT_"+ (string)spar1 + "_" + (string)spar2 + "_" +(string)spar3 + "_"+  "widthAllTo10.root").c_str());
       c1.SaveAs(("output/"+(string)dir+"/"+ (string)slayerCutDown+ "lt" +(string)subDet +"mt"+(string)slayerCutUp +"_CT_"+ (string)spar1 + "_" + (string)spar2 + "_" +(string)spar3 + "_"+  "widthAllTo10.eps").c_str());

      /////////////////////////////////////////////////////




    //cout << "in here 7"  << endl;

   return 0;
}
