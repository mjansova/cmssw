//./distributionsChi2TestGeneral test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC252Tuned0.805step3.root 252 TOB 0 5 yes CTtuning 
//./distributionsChi2TestGeneral test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC257Tuned0.805step3.root 257 TOB 0 5 yes CTtuning 
//./distributionsChi2TestGeneral test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC262Tuned0.805step3.root 262 TOB 0 5 yes CTtuning 
//./distributionsChi2TestGeneral test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC267Tuned0.805step3.root 267 TOB 0 5 yes CTtuning 
//
//-rw-r--r-- 1 mjansova cms 22213879  9 oct.  19:41 test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_min20_252.0step3.root
//-rw-r--r-- 1 mjansova cms 23048735  9 oct.  19:42 test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC252Tuned0.805step3.root
//-rw-r--r-- 1 mjansova cms 22153780  9 oct.  20:48 test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_min10_257.0step3.root
//-rw-r--r-- 1 mjansova cms 22965898  9 oct.  20:48 test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC257Tuned0.805step3.root
//-rw-r--r-- 1 mjansova cms 22079516  9 oct.  21:54 test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_pl10_262.0step3.root
//-rw-r--r-- 1 mjansova cms 22859092  9 oct.  21:55 test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC262Tuned0.805step3.root
//-rw-r--r-- 1 mjansova cms 21986281  9 oct.  22:59 test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_pl20_267.0step3.root
//-rw-r--r-- 1 mjansova cms 22786514  9 oct.  23:01 test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC267Tuned0.805step3.root


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
    if(argc != 9)
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
    TString subDet = argv[4];
    TString slayerCutDown = argv[5];
    TString slayerCutUp = argv[6];
    TString norm = argv[7];
    TString dir = argv[8];

    float par2 = spar2.Atof();

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



       TH1F* chargeData = new TH1F("chargeForAllWidths", "chargeForAllWidths" , 50, 100, 500 );
       TH1F* chargeMC = new TH1F("chargeForAllWidthsMC", "chargeForAllWidthsMC" , 50, 100, 500 );
       TH1F* chargeDataPart1 = new TH1F("chargeForAllWidthsPart1", "chargeForAllWidthsPart1" , 25, 100, 230 );
       TH1F* chargeMCPart1 = new TH1F("chargeForAllWidthsMCPart1", "chargeForAllWidthsMCPart1" , 25, 100, 230 );
       TH1F* chargeDataPart2 = new TH1F("chargeForAllWidthsPart2", "chargeForAllWidthsPart2" , 25, 230, 500 );
       TH1F* chargeMCPart2 = new TH1F("chargeForAllWidthsMCPart2", "chargeForAllWidthsMCPart2" , 25, 230, 500 );
       TH1F* widthData = new TH1F("widthAllTo10Data", "widthAllTo10Data" , 10, 0, 10  );
       TH1F* widthMC = new TH1F("widthAllTo10MC", "widthAllTo10MC" , 10, 0, 10 );
       TH1F* widthAllTo15Data = new TH1F("widthAllTo15Data", "widthAllTo15Data" , 14, 1, 15  );
       TH1F* widthAllTo15MC = new TH1F("widthAllTo15MC", "widthAllTo15MC" , 14, 1, 15 );
       TH1F* widthAllTo20Data = new TH1F("widthAllTo20Data", "widthAllTo20Data" , 19, 1, 20  );
       TH1F* widthAllTo20MC = new TH1F("widthAllTo20MC", "widthAllTo20MC" , 19, 1, 20 );
   
       cout << "in here 2"  << endl;
       for(uint32_t m = 0; m<subclusterwidth.size(); m++)
       {
           if(layerCutDown< subclusterlayerwheel.at(m) && subclusterlayerwheel.at(m) < layerCutUp)
           {
               chargeData->Fill(subclustercharge.at(m));
               chargeDataPart1->Fill(subclustercharge.at(m));
               chargeDataPart2->Fill(subclustercharge.at(m));
               widthData->Fill(subclusterwidth.at(m));
               widthAllTo15Data->Fill(subclusterwidth.at(m));
               widthAllTo20Data->Fill(subclusterwidth.at(m));
           }
         
       }
   
    cout << "in here 3"  << endl;
       for(uint32_t m = 0; m<subclusterwidth2.size(); m++)
       {
           if(layerCutDown< subclusterlayerwheel2.at(m) &&  subclusterlayerwheel2.at(m)< layerCutUp)
           {
               chargeMC->Fill(subclustercharge2.at(m));
               chargeMCPart1->Fill(subclustercharge2.at(m));
               chargeMCPart2->Fill(subclustercharge2.at(m));
               widthMC->Fill(subclusterwidth2.at(m));
               widthAllTo15MC->Fill(subclusterwidth2.at(m));
               widthAllTo20MC->Fill(subclusterwidth2.at(m));
           }
       }


       
       float Ndata=  chargeData->Integral();
       float NMC=  chargeMC->Integral();
       float scale = Ndata/NMC;
       chargeMC->Scale(scale);

       float NdataW=  widthData->Integral();
       float NMCW=  widthMC->Integral();
       float scaleW = NdataW/NMCW;
       widthMC->Scale(scaleW);

       /*for(uint32_t c=0; c<widthAllTo10MC->GetNbinsX(); c++)
       {
           cout << "MC bin " << c+1 << " n events " << widthAllTo10MC->GetBinContent(c+1) << endl;
       }*/

       TFile *f_chi2 = new TFile (("output/"+(string)dir+"/"+(string)dir+"chi2.root").c_str(), "UPDATE");
       TH1F* chi21D = NULL;
       chi21D = dynamic_cast<TH1F*>(f_chi2->Get("chi21D"));
       if(chi21D == NULL)
       {
           cout << "new chi2 histogram created " << endl;
           chi21D = new TH1F("chi21D", "chi21D", 1000, 190, 300);
       }

       /*TFile *f_chi2Part1 = new TFile (("output/"+(string)dir+"/"+(string)dir+"chi2Part1.root").c_str(), "UPDATE");
       TH1F* chi21DPart1 = NULL;
       chi21DPart1 = dynamic_cast<TH1F*>(f_chi2Part1->Get("chi21DPart1"));
       if(chi21DPart1 == NULL)
       {
           cout << "new chi2 histogram created " << endl;
           chi21DPart1 = new TH1F("chi21DPart1", "chi21DPart1", 1000, 190, 300);
       }

       TFile *f_chi2Part2 = new TFile (("output/"+(string)dir+"/"+(string)dir+"chi2Part2.root").c_str(), "UPDATE");
       TH1F* chi21DPart2 = NULL;
       chi21DPart2 = dynamic_cast<TH1F*>(f_chi2Part2->Get("chi21DPart2"));
       if(chi21DPart2 == NULL)
       {
           cout << "new chi2 histogram created " << endl;
           chi21DPart2 = new TH1F("chi21DPart2", "chi21DPart2", 1000, 190, 300);
       }*/

       
       TAxis* xaxis = chi21D->GetXaxis();
       int binx = xaxis->FindBin(par2);

       chi21D->SetBinContent(binx, chargeData->Chi2Test(chargeMC, "UW CHI2/NDF P"));
       chi21D->Write(chi21D->GetName(),TObject::kOverwrite);
       f_chi2->Close();

/*       chi21DPart1->SetBinContent(binx, chargeDataPart1->Chi2Test(chargeMCPart1, "UW CHI2/NDF"));
       chi21DPart1->Write(chi21DPart1->GetName(),TObject::kOverwrite);
       f_chi2Part1->Close();

       chi21DPart2->SetBinContent(binx, chargeDataPart2->Chi2Test(chargeMCPart2, "UW CHI2/NDF"));
       chi21DPart2->Write(chi21DPart2->GetName(),TObject::kOverwrite);
       f_chi2Part2->Close();
*/
       TFile *f_Kolmogorov = new TFile (("output/"+(string)dir+"/"+(string)dir+"Kolmogorov.root").c_str(), "UPDATE");
       TH1F* kol1D = NULL;
       kol1D = dynamic_cast<TH1F*>(f_Kolmogorov->Get("kol1D"));
       if(kol1D == NULL)
       {
           cout << "new kol histogram created " << endl;
           kol1D = new TH1F("kol1D", "kol1D", 1000,190,300 );
       }
       kol1D->SetBinContent(binx, chargeData->KolmogorovTest(chargeMC));
       kol1D->Write(kol1D->GetName(),TObject::kOverwrite);
       f_Kolmogorov->Close();

/*       cout << "chi2 widthAllTo10 " << widthAllTo10Data->Chi2Test(widthAllTo10MC, "CHI2") << endl;;
       cout << "Kolmogorov widthAllTo10 " << widthAllTo10Data->KolmogorovTest(widthAllTo10MC) << endl;;
       cout << "chi2 widthAllTo15 " << widthAllTo15Data->Chi2Test(widthAllTo15MC, "CHI2") << endl;;
       cout << "Kolmogorov widthAllTo15 " << widthAllTo15Data->KolmogorovTest(widthAllTo15MC) << endl;;
       cout << "chi2 widthAllTo20 " << widthAllTo20Data->Chi2Test(widthAllTo20MC, "CHI2") << endl;;
       cout << "Kolmogorov widthAllTo20 " << widthAllTo20Data->KolmogorovTest(widthAllTo20MC) << endl;;
*/

//////////////////////filling
       chargeData->SetMarkerStyle(kFullCircle);
       chargeMC->SetMarkerColor(kBlack);
       widthData->SetMarkerStyle(kFullCircle);
       widthData->SetMarkerColor(kBlack);


       chargeMC->SetTitle("");
       chargeData->SetTitle("");
       widthData->SetTitle("");


       chargeData->GetXaxis()->SetTitle("cluster charge");
       widthData->GetXaxis()->SetTitle("cluster width");

       chargeData->SetMaximum(1.5*  (chargeData->GetMaximum()));
       widthData->SetMaximum(1.5*  (widthData->GetMaximum()));

       chargeMC->SetLineColor(kRed);
       widthMC->SetLineColor(kRed);


       TCanvas c1("chargeForAllWidths","chargeForAllWidths");
       TPad pad1("pad1", "The pad 80% of the height",0.0,0.2,1.0,1.0);
       TPad pad2("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.2);
       pad1.Draw();
       pad2.Draw();
       
       pad1.cd();
       chargeData->Draw("P");
       chargeMC->Draw("same"); 
       pad2.cd();
       TH1F* chargeDataClone = (TH1F*) chargeData->Clone();
       chargeDataClone->Divide(chargeMC);
       chargeDataClone->SetMaximum(2);
       chargeDataClone->SetMinimum(0);
       chargeDataClone->Draw("e");

       TCanvas c2("widthForAllWidths","widthForAllWidths");
       TPad pad3("pad3", "The pad 80% of the height",0.0,0.2,1.0,1.0);
       TPad pad4("pad4", "The pad 20% of the height",0.0,0.0,1.0,0.2);
       pad3.Draw();
       pad4.Draw();
       
       pad3.cd();
       widthData->Draw("P");
       widthMC->Draw("same"); 
       pad4.cd();
       TH1F* widthDataClone = (TH1F*) widthData->Clone();
       widthDataClone->Divide(widthMC);
       widthDataClone->SetMaximum(2);
       widthDataClone->SetMinimum(0);
       widthDataClone->Draw("e");

       c1.SaveAs(("output/"+(string)dir+"/"+ (string)slayerCutDown+ "lt" +(string)subDet +"mt"+(string)slayerCutUp +"_eToADC_"+ (string)spar2  + "_"+  "charge.root").c_str());
       c1.SaveAs(("output/"+(string)dir+"/"+ (string)slayerCutDown+ "lt" +(string)subDet +"mt"+(string)slayerCutUp +"_eToADC_"+ (string)spar2  +"_" +  "charge.eps").c_str());
       c2.SaveAs(("output/"+(string)dir+"/"+ (string)slayerCutDown+ "lt" +(string)subDet +"mt"+(string)slayerCutUp +"_eToADC_"+ (string)spar2  + "_"+  "width.root").c_str());
       c2.SaveAs(("output/"+(string)dir+"/"+ (string)slayerCutDown+ "lt" +(string)subDet +"mt"+(string)slayerCutUp +"_eToADC_"+ (string)spar2  +"_" +  "width.eps").c_str());

      /////////////////////////////////////////////////////




    //cout << "in here 7"  << endl;

   return 0;
}
