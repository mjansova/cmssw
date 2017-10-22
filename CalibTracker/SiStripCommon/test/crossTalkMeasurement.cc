//./crossTalkMeasurement test_shallowTrackCRUZET_2017VRdigis.root crossTalkMeasurement //small stats
//./crossTalkMeasurement test_shallowTrackCRUZET_2017VRNewClusteringTimingTestFULL.root crossTalkMeasurement
//./crossTalkMeasurement test_shallowTrackCRUZET_2017VRNewClusteringTimingTest.root crossTalkMeasurement

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


       vector<float>* CTstripChargeSubdetid = 0;
       vector<float>* CTstripCharge = 0;
       vector<float>* CTstripChargeLayerwheel = 0;
       //vector<float>* CTstripChargeStripNr = 0;
       vector<float>* CTstripChargeTotWidth = 0;
       vector<float>* CTstripChargeTotCharge = 0;
       vector<float>* CTstripChargeLocalTrackPhi = 0;
       //vector<float>* CTstripChargeGlobalTrackPhi = 0;
       vector<float>* CTstripChargeLocalTrackTheta = 0;
       //vector<float>* CTstripChargeGlobalTrackTheta = 0;
       //vector<unsigned>* CTstripChargeDetid = 0;
       //vector<float>* CTstripChargeLocalX = 0;
       //vector<float>* CTstripChargeLocalY = 0;
       //vector<float>* CTstripChargetrackPt = 0;
       vector<float>* CTstripChargelocalpitch = 0;
       vector<float>* CTstripChargesensorThickness = 0;
       vector<float>* CTCmbtimeVtx = 0;
       vector<float>* CTCmbtimeVtxr = 0;
       vector<float>* CTCmbtimeVtxErr = 0;
       vector<float>* CTCmbtimeVtxrErr = 0;
       vector<float>* CTMuontrackDirection = 0;
       //vector<float>* CTstripChargeBdotY = 0;
       //unsigned int  nroftracks;
       //unsigned int  nrofevents;
       vector<float>* tsosrhglobalphi = 0;
       vector<uint32_t>* clusterdetid = 0;

       vector<float> subCTstripChargeSubdetid;
       vector<float> subCTstripCharge;
       vector<float> subCTstripChargeLayerwheel;
       //vector<float> subCTstripChargeStripNr;
       vector<float> subCTstripChargeTotWidth;
       vector<float> subCTstripChargeTotCharge;
       vector<float> subCTstripChargeLocalTrackPhi;
       //vector<float> subCTstripChargeGlobalTrackPhi;
       vector<float> subCTstripChargeLocalTrackTheta;
       //vector<float> subCTstripChargeGlobalTrackTheta;
       //vector<unsigned> subCTstripChargeDetid;
       //vector<float> subCTstripChargeLocalX;
       //vector<float> subCTstripChargeLocalY;
       //vector<float> subCTstripChargetrackPt;
       vector<float> subCTstripChargelocalpitch;
       vector<float> subCTstripChargesensorThickness;
       vector<float> subCTCmbtimeVtx;
       vector<float> subCTCmbtimeVtxr;
       vector<float> subCTCmbtimeVtxErr;
       vector<float> subCTCmbtimeVtxrErr;
       vector<float> subCTMuontrackDirection;
       //vector<float> subCTstripChargeBdotY;
       //vector<unsigned int>  subnroftracks;
       //vector<unsigned int>  subnrofevents;
       vector<float> subtsosrhglobalphi;
       vector<uint32_t> subclusterdetid;

       float bx;
       vector<float> bxPerStrip;

       t1->SetBranchAddress("CTstripChargeSubdetid",  &CTstripChargeSubdetid );
       t1->SetBranchAddress("CTstripCharge",  &CTstripCharge );
       t1->SetBranchAddress("CTstripChargeLayerwheel",  &CTstripChargeLayerwheel );
       //t1->SetBranchAddress("CTstripChargeStripNr",  &CTstripChargeStripNr );
       t1->SetBranchAddress("CTstripChargeTotWidth",  &CTstripChargeTotWidth );
       t1->SetBranchAddress("CTstripChargeTotCharge",  &CTstripChargeTotCharge );
       t1->SetBranchAddress("CTstripChargeLocalTrackPhi",  &CTstripChargeLocalTrackPhi );
       //t1->SetBranchAddress("CTstripChargeGlobalTrackPhi",  &CTstripChargeGlobalTrackPhi );
       t1->SetBranchAddress("CTstripChargeLocalTrackTheta",  &CTstripChargeLocalTrackTheta );
       //t1->SetBranchAddress("CTstripChargeGlobalTrackTheta",  &CTstripChargeGlobalTrackTheta );
       //t1->SetBranchAddress("CTstripChargeDetid",  &CTstripChargeDetid );
       //t1->SetBranchAddress("CTstripChargeLocalX",  &CTstripChargeLocalX );
       //t1->SetBranchAddress("CTstripChargeLocalY",  &CTstripChargeLocalY );
       //t1->SetBranchAddress("CTstripChargetrackPt",  &CTstripChargetrackPt );
       t1->SetBranchAddress("CTstripChargelocalpitch",  &CTstripChargelocalpitch );
       t1->SetBranchAddress("CTstripChargesensorThickness",  &CTstripChargesensorThickness );
       t1->SetBranchAddress("CTCmbtimeVtx",  &CTCmbtimeVtx );
       t1->SetBranchAddress("CTCmbtimeVtxr",  &CTCmbtimeVtxr );
       t1->SetBranchAddress("CTCmbtimeVtxErr",  &CTCmbtimeVtxErr );
       t1->SetBranchAddress("CTCmbtimeVtxrErr",  &CTCmbtimeVtxrErr );
       t1->SetBranchAddress("CTMuontrackDirection",  &CTMuontrackDirection );
       //t1->SetBranchAddress("CTstripChargeBdotY",  &CTstripChargeBdotY );
       //t1->SetBranchAddress("nroftracks",  &nroftracks );
       //t1->SetBranchAddress("nrofevents",  &nrofevents );
       t1->SetBranchAddress("tsosrhglobalphi",  &tsosrhglobalphi );
       t1->SetBranchAddress("clusterdetid",  &clusterdetid );
       //t1->SetBranchAddress("bx",  &bx );
   //data always first
    

   uint32_t evCount=0;
   vector<uint32_t> eventCount;
   
   //cout << "in here a" << endl;
   Int_t nentries = (Int_t)t1->GetEntries();

   //cout << "in here b" << endl;
   ///fill variables from tree 1
   for (Int_t e=0; e<nentries; e++) 
   {
       t1->GetEntry(e);
          
           //per cluster
           //subnroftracks.push_back(nroftracks);
           //subnrofevents.push_back(nrofevents);
           //perstrip
           uint32_t upStrip = CTstripChargeSubdetid->size();
           for(uint32_t k=0; k<upStrip;k++)
           {
               if( true)
               {
                   if(true ) //@MJ@ TODO trigger
                   {
                       subCTstripChargeSubdetid.push_back(CTstripChargeSubdetid->at(k));
                       subCTstripCharge.push_back(CTstripCharge->at(k));
                       subCTstripChargeLayerwheel.push_back(CTstripChargeLayerwheel->at(k));
                       //subCTstripChargeStripNr.push_back(CTstripChargeStripNr->at(k));
                       subCTstripChargeTotWidth.push_back(CTstripChargeTotWidth->at(k));
                       subCTstripChargeTotCharge.push_back(CTstripChargeTotCharge->at(k));
                       subCTstripChargeLocalTrackPhi.push_back(CTstripChargeLocalTrackPhi->at(k));
                       //subCTstripChargeGlobalTrackPhi.push_back(CTstripChargeGlobalTrackPhi->at(k));
                       subCTstripChargeLocalTrackTheta.push_back(CTstripChargeLocalTrackTheta->at(k));
                       //subCTstripChargeGlobalTrackTheta.push_back(CTstripChargeGlobalTrackTheta->at(k));
                       //subCTstripChargeDetid.push_back(CTstripChargeDetid->at(k));
                       //subCTstripChargeLocalX.push_back(CTstripChargeLocalX->at(k));
                       //subCTstripChargeLocalY.push_back(CTstripChargeLocalY->at(k));
                       //subCTstripChargetrackPt.push_back(CTstripChargetrackPt->at(k));
                       subCTstripChargelocalpitch.push_back(CTstripChargelocalpitch->at(k));
                       subCTstripChargesensorThickness.push_back(CTstripChargesensorThickness->at(k));
                       subCTCmbtimeVtx.push_back(CTCmbtimeVtx->at(k));
                       subCTCmbtimeVtxr.push_back(CTCmbtimeVtxr->at(k));
                       subCTCmbtimeVtxErr.push_back(CTCmbtimeVtxErr->at(k));
                       subCTCmbtimeVtxrErr.push_back(CTCmbtimeVtxrErr->at(k));
                       subCTMuontrackDirection.push_back(CTMuontrackDirection->at(k));
                       //subCTstripChargeBdotY.push_back(CTstripChargeBdotY->at(k));
                       //bxPerStrip.push_back(bx);
                       eventCount.push_back(e);
                   }
               }
           }
           for(uint32_t j=0; j<tsosrhglobalphi->size() ;j++)//@MJ@ TODO just temporaray
           {
               subtsosrhglobalphi.push_back(tsosrhglobalphi->at(j));
               subclusterdetid.push_back(clusterdetid->at(j));
           }
   }



       //F
       vector<TH1F*> narrowTrackSharing1Data;
       vector<TH1F*> narrowTrackSharing2Data;

       vector<TH1F*> timingErrorTop;
       vector<TH1F*> timingErrorTopInOut;
       vector<TCanvas*> timingErrorTopCan;
       vector<TCanvas*> timingErrorTopInOutCan;

       vector<TProfile*> widthTopInOut;
       vector<TProfile*> widthBottomInOut;
       vector<TProfile*> widthTopOutIn;
       vector<TProfile*> widthBottomOutIn;

       vector<TProfile*> chargeTopInOut;
       vector<TProfile*> chargeBottomInOut;
       vector<TProfile*> chargeTopOutIn;
       vector<TProfile*> chargeBottomOutIn;

       vector<TCanvas*> narrowTrackSharing1DataCan;
       vector<TCanvas*> narrowTrackSharing2DataCan;

       vector<TCanvas*> widthTopInOutCan;
       vector<TCanvas*> widthBottomInOutCan;
       vector<TCanvas*> widthTopOutInCan;
       vector<TCanvas*> widthBottomOutInCan;

       vector<TCanvas*> chargeTopInOutCan;
       vector<TCanvas*> chargeBottomInOutCan;
       vector<TCanvas*> chargeTopOutInCan;
       vector<TCanvas*> chargeBottomOutInCan;


       vector<TProfile*> widthTopInOutOnly;
       vector<TProfile*> widthBottomInOutOnly;
       vector<TProfile*> chargeTopInOutOnly;
       vector<TProfile*> chargeBottomInOutOnly;
       vector<TCanvas*> widthTopInOutOnlyCan;
       vector<TCanvas*> widthBottomInOutOnlyCan;
       vector<TCanvas*> chargeTopInOutOnlyCan;
       vector<TCanvas*> chargeBottomInOutOnlyCan;

       narrowTrackSharing1Data.resize(20, NULL);
       narrowTrackSharing2Data.resize(20, NULL);
       narrowTrackSharing1DataCan.resize(20, NULL);
       narrowTrackSharing2DataCan.resize(20, NULL);

       widthTopInOut.resize(20, NULL);
       widthBottomInOut.resize(20, NULL);
       widthTopInOutOnly.resize(20, NULL);
       widthBottomInOutOnly.resize(20, NULL);
       widthTopOutIn.resize(20, NULL);
       widthBottomOutIn.resize(20, NULL);

       chargeTopInOut.resize(20, NULL);
       chargeBottomInOut.resize(20, NULL);
       chargeTopInOutOnly.resize(20, NULL);
       chargeBottomInOutOnly.resize(20, NULL);
       chargeTopOutIn.resize(20, NULL);
       chargeBottomOutIn.resize(20, NULL);

       widthTopInOutCan.resize(20, NULL);
       widthBottomInOutCan.resize(20, NULL);
       widthTopInOutOnlyCan.resize(20, NULL);
       widthBottomInOutOnlyCan.resize(20, NULL);
       widthTopOutInCan.resize(20, NULL);
       widthBottomOutInCan.resize(20, NULL);

       chargeTopInOutCan.resize(20, NULL);
       chargeBottomInOutCan.resize(20, NULL);
       chargeTopInOutOnlyCan.resize(20, NULL);
       chargeBottomInOutOnlyCan.resize(20, NULL);
       chargeTopOutInCan.resize(20, NULL);
       chargeBottomOutInCan.resize(20, NULL);


       timingErrorTop.resize(20, NULL);
       timingErrorTopInOut.resize(20, NULL);
       timingErrorTopCan.resize(20, NULL);
       timingErrorTopInOutCan.resize(20, NULL);


       uint8_t TIDoffset = 3;
       uint8_t TECoffset = 6;

       vector<string> canvases;

       for(uint32_t m = 0; m<subCTstripCharge.size()/5; m++) //@MJ@ TODO nr of clusters in here, correct variables
       {
           int32_t clusterStart = m*5;


           int32_t partPos = -1;
           string  parName;
               if(subCTstripChargeSubdetid.at(clusterStart) == 3 ) //TIB
               {
                   if(subCTstripChargeLayerwheel.at(clusterStart)>0 && subCTstripChargeLayerwheel.at(clusterStart)<3)
                   {
                       partPos = 0;
                       parName = "TIB1";
                   }
                   else
                   {
                       partPos = 1;
                       parName = "TIB2";

                   }

               }
               if(subCTstripChargeSubdetid.at(clusterStart) == 5 ) //TOB
               {
                   if(subCTstripChargeLayerwheel.at(clusterStart)>0 && subCTstripChargeLayerwheel.at(clusterStart)<5)
                   {
                       partPos = 2;
                       parName = "TOB1";
                   }
                   else
                   {
                       partPos = 3;
                       parName = "TOB2";
                   }

               }
               if(subCTstripChargeSubdetid.at(clusterStart) == 4 ) //TID
               {
                       //cout << "TID filling bin " << TIDoffset + subCTstripChargeLayerwheel.at(clusterStart) << endl;
                       partPos = TIDoffset + subCTstripChargeLayerwheel.at(clusterStart);
                       parName = "TID"+to_string(partPos);
               }
               if(subCTstripChargeSubdetid.at(clusterStart) == 6 ) //TEC
               {
                       //cout << "TEC filling bin " << TECoffset + subCTstripChargeLayerwheel.at(clusterStart) << " for layer " << subCTstripChargeLayerwheel.at(clusterStart) << endl;
                       partPos = TECoffset + subCTstripChargeLayerwheel.at(clusterStart);
                       parName = "TEC"+to_string(partPos);
               }


           //is the cluster narrow?           
           if( abs(tan(subCTstripChargeLocalTrackTheta.at(clusterStart))*cos(subCTstripChargeLocalTrackPhi.at(clusterStart))) < narrowness*(subCTstripChargelocalpitch.at(clusterStart)/subCTstripChargesensorThickness.at(clusterStart)) )
           {
               {
                   if(narrowTrackSharing1Data.at(partPos) == NULL)
                   {
                       narrowTrackSharing1Data.at(partPos) = new TH1F( ("narrowTrackSharing1"+ parName).c_str() , ("narrowTrackSharing1"+parName).c_str() , 200, -1, 1 );
                       narrowTrackSharing1DataCan.at(partPos) = new TCanvas( ("cnarrowTrackSharing1"+ parName).c_str() , ("cnarrowTrackSharing1"+parName).c_str()  );
                       canvases.push_back(parName);
                   }
                   if(narrowTrackSharing2Data.at(partPos) == NULL)
                   {
                       narrowTrackSharing2Data.at(partPos) = new TH1F( ("narrowTrackSharing2"+ parName).c_str() , ("narrowTrackSharing2"+parName).c_str() , 200, -1, 1 );
                       narrowTrackSharing2DataCan.at(partPos) = new TCanvas( ("cnarrowTrackSharing2"+ parName).c_str() , ("cnarrowTrackSharing2"+parName).c_str());
                   }
                   //@MJ@ TODO create histogram if not done yet
               if(subCTstripCharge.at(clusterStart+1) != -333)
                   narrowTrackSharing1Data.at(partPos)->Fill((float) subCTstripCharge.at(clusterStart+1)/subCTstripCharge.at(clusterStart+2));
               if(subCTstripCharge.at(clusterStart+3) != -333)
                   narrowTrackSharing1Data.at(partPos)->Fill( (float) subCTstripCharge.at(clusterStart+3)/subCTstripCharge.at(clusterStart+2));
               if(subCTstripCharge.at(clusterStart) != -333)
                   narrowTrackSharing2Data.at(partPos)->Fill( (float) subCTstripCharge.at(clusterStart)/subCTstripCharge.at(clusterStart+2));
               if(subCTstripCharge.at(clusterStart+4) != -333)
                   narrowTrackSharing2Data.at(partPos)->Fill( (float) subCTstripCharge.at(clusterStart+4)/subCTstripCharge.at(clusterStart+2));


                   if(widthTopInOut.at(partPos) == NULL)
                   {
		       widthTopInOut.at(partPos) =  new TProfile(("widthTopInOut"+ parName).c_str(), ("widthTopInOut"+ parName).c_str() , 200, -100, 100, 0, 10 ) ;
		       widthBottomInOut.at(partPos) = new TProfile(("widthBottomInOut"+ parName).c_str(), ("widthBottomInOut"+ parName).c_str() , 200, -100, 100, 0, 10 ) ;
		       widthTopInOutOnly.at(partPos) =  new TProfile(("widthTopInOutOnly"+ parName).c_str(), ("widthTopInOutOnly"+ parName).c_str() , 200, -100, 100, 0, 10 ) ;
		       widthBottomInOutOnly.at(partPos) = new TProfile(("widthBottomInOutOnly"+ parName).c_str(), ("widthBottomInOutOnly"+ parName).c_str() , 200, -100, 100, 0, 10 ) ;
		       widthTopOutIn.at(partPos) =  new TProfile(("widthTopOutIn"+ parName).c_str(), ("widthTopOutIn"+ parName).c_str() , 200, -100, 100, 0, 10 ) ;
		       widthBottomOutIn.at(partPos) =  new TProfile(("widthBottomOutIn"+ parName).c_str(), ("widthBottomOutIn"+ parName).c_str() , 200, -100, 100, 0, 10 ) ;

		       chargeTopInOut.at(partPos) =  new TProfile(("chargeTopInOut"+ parName).c_str(), ("chargeTopInOut"+ parName).c_str() , 200, -100, 100, 0, 1000 ) ;
		       chargeBottomInOut.at(partPos) =  new TProfile(("chargeBottomInOut"+ parName).c_str(), ("chargeBottomInOut"+ parName).c_str() , 200, -100, 100, 0, 1000 ) ;
		       chargeTopInOutOnly.at(partPos) =  new TProfile(("chargeTopInOutOnly"+ parName).c_str(), ("chargeTopInOutOnly"+ parName).c_str() , 200, -100, 100, 0, 1000 ) ;
		       chargeBottomInOutOnly.at(partPos) =  new TProfile(("chargeBottomInOutOnly"+ parName).c_str(), ("chargeBottomInOutOnly"+ parName).c_str() , 200, -100, 100, 0, 1000 ) ;
		       chargeTopOutIn.at(partPos)=  new TProfile(("chargeTopOutIn"+ parName).c_str(), ("chargeTopOutIn"+ parName).c_str() , 200, -100, 100, 0, 1000 ) ;
		       chargeBottomOutIn.at(partPos) =  new TProfile(("chargeBottomOutIn"+ parName).c_str(), ("chargeBottomOutIn"+ parName).c_str() , 200, -100, 100, 0, 1000 ) ;

		       widthTopInOutCan.at(partPos) =   new TCanvas( ("widthTopInOut"+ parName).c_str(), ("widthTopInOut"+ parName).c_str()) ;
		       widthBottomInOutCan.at(partPos) = new TCanvas( ("widthBottomInOut"+ parName).c_str(), ("widthBottomInOut"+ parName).c_str()) ;
		       widthTopInOutOnlyCan.at(partPos) =   new TCanvas( ("widthTopInOutOnly"+ parName).c_str(), ("widthTopInOutOnly"+ parName).c_str()) ;
		       widthBottomInOutOnlyCan.at(partPos) = new TCanvas( ("widthBottomInOutOnly"+ parName).c_str(), ("widthBottomInOutOnly"+ parName).c_str()) ;
		       widthTopOutInCan.at(partPos) =   new TCanvas( ("widthTopOutIn"+ parName).c_str(), ("widthTopOutIn"+ parName).c_str()) ;
		       widthBottomOutInCan.at(partPos) =   new TCanvas( ("widthBottomOutIn"+ parName).c_str(), ("widthBottomOutIn"+ parName).c_str()) ;

		       chargeTopInOutCan.at(partPos) =   new TCanvas( ("chargeTopInOut"+ parName).c_str(), ("chargeTopInOut"+ parName).c_str()) ;
		       chargeBottomInOutCan.at(partPos) =  new TCanvas( ("chargeBottomInOut"+ parName).c_str(), ("chargeBottomInOut"+ parName).c_str()) ;
		       chargeTopInOutOnlyCan.at(partPos) =   new TCanvas( ("chargeTopInOutOnly"+ parName).c_str(), ("chargeTopInOutOnly"+ parName).c_str()) ;
		       chargeBottomInOutOnlyCan.at(partPos) =  new TCanvas( ("chargeBottomInOutOnly"+ parName).c_str(), ("chargeBottomInOutOnly"+ parName).c_str()) ;
		       chargeTopOutInCan.at(partPos)=   new TCanvas( ("chargeTopOutIn"+ parName).c_str(), ("chargeTopOutIn"+ parName).c_str()) ;
		       chargeBottomOutInCan.at(partPos) = new TCanvas( ("chargeBottomOutIn"+ parName).c_str(), ("chargeBottomOutIn"+ parName).c_str()) ;

                       timingErrorTop.at(partPos) = new TH1F( ("timingErrorTop"+ parName).c_str() , ("timingErrorTop"+parName).c_str() , 100, 0, 20 );
                       timingErrorTopInOut.at(partPos)  = new TH1F( ("timingErrorTopInOut"+ parName).c_str() , ("timingErrorTopInOut"+parName).c_str() , 100, 0, 20 );
                       timingErrorTopCan.at(partPos) = new TCanvas( ("timingErrorTop"+ parName).c_str() , ("timingErrorTop"+parName).c_str() ) ;
                       timingErrorTopInOutCan.at(partPos) = new TCanvas( ("timingErrorTopInOut"+ parName).c_str() , ("timingErrorTopInOut"+parName).c_str()  ) ;


                   }

                   if(true) //if(subclusterlayerwheel.at(m) == 3)
                   {
                       cout << "partition name " << parName << " event count " << eventCount.at(clusterStart) << "layer " << subCTstripChargeLayerwheel.at(clusterStart) << " detid " << subclusterdetid.at(m) << endl;
                       if(subtsosrhglobalphi.at(m)>0) //top
                       {
                           widthTopInOutOnly.at(partPos)->Fill(subCTCmbtimeVtx.at(clusterStart), subCTstripChargeTotWidth.at(clusterStart));
                           chargeTopInOutOnly.at(partPos)->Fill(subCTCmbtimeVtx.at(clusterStart), subCTstripChargeTotCharge.at(clusterStart));
                           if(subCTMuontrackDirection.at(clusterStart) > 0)
                           {
                               cout << "top In Out "  << endl;
                               widthTopInOut.at(partPos)->Fill(subCTCmbtimeVtx.at(clusterStart), subCTstripChargeTotWidth.at(clusterStart));
                               chargeTopInOut.at(partPos)->Fill(subCTCmbtimeVtx.at(clusterStart), subCTstripChargeTotCharge.at(clusterStart));
                           }
                           else if(subCTMuontrackDirection.at(clusterStart) < 0)
                           {
                               cout << "top Out In "  << endl;
                               widthTopOutIn.at(partPos)->Fill(subCTCmbtimeVtxr.at(clusterStart), subCTstripChargeTotWidth.at(clusterStart));
                               chargeTopOutIn.at(partPos)->Fill(subCTCmbtimeVtxr.at(clusterStart), subCTstripChargeTotCharge.at(clusterStart));
                               timingErrorTop.at(partPos)->Fill(subCTCmbtimeVtxrErr.at(clusterStart));
                               timingErrorTopInOut.at(partPos)->Fill(subCTCmbtimeVtxErr.at(clusterStart));
                           }
                           else
                               cout << "direction undefined" << endl;
                       }
                       else if(subtsosrhglobalphi.at(m)<0)
                       {
                           widthBottomInOutOnly.at(partPos)->Fill(subCTCmbtimeVtx.at(clusterStart), subCTstripChargeTotWidth.at(clusterStart));
                           chargeBottomInOutOnly.at(partPos)->Fill(subCTCmbtimeVtx.at(clusterStart), subCTstripChargeTotCharge.at(clusterStart));
                           if(subCTMuontrackDirection.at(clusterStart) > 0)
                           {
                               cout << "bottom In Out "  << endl;
                               widthBottomInOut.at(partPos)->Fill(subCTCmbtimeVtx.at(clusterStart), subCTstripChargeTotWidth.at(clusterStart));
                               chargeBottomInOut.at(partPos)->Fill(subCTCmbtimeVtx.at(clusterStart), subCTstripChargeTotCharge.at(clusterStart));
                           }
                           else if(subCTMuontrackDirection.at(clusterStart) < 0)
                           {
                               cout << "bottom Out In "  << endl;
                               widthBottomOutIn.at(partPos)->Fill(subCTCmbtimeVtxr.at(clusterStart), subCTstripChargeTotWidth.at(clusterStart));
                               chargeBottomOutIn.at(partPos)->Fill(subCTCmbtimeVtxr.at(clusterStart), subCTstripChargeTotCharge.at(clusterStart));
                           }
                           else
                               cout << "direction undefined" << endl;

                       }
                       else
                       {
                           cout << "weird" << endl;
                       }
                   }
               }
           }

       }

       TH1F* dummyCTplot = new TH1F("dummyCTplot", "dummyCTplot" , 5, 0, 5 );
       dummyCTplot->SetBinContent(1,1);
       dummyCTplot->SetBinContent(2,2);
       dummyCTplot->SetBinContent(3,3);
       dummyCTplot->SetBinContent(4,2);
       dummyCTplot->SetBinContent(5,1);
       dummyCTplot->SetFillStyle(1001);
       dummyCTplot->SetFillColor(kGray+2);
       dummyCTplot->SetTitle("");
       dummyCTplot->SetMaximum(1.5*(dummyCTplot->GetMaximum()));
       
       

//////////////////////filling
//     
//
       for(uint32_t c= 0; c<canvases.size(); c++ )
       {
           if(narrowTrackSharing1Data.at(c)==NULL || narrowTrackSharing2Data.at(c) == NULL)
               break;

       narrowTrackSharing1Data.at(c)->SetMarkerStyle(kFullCircle);
       narrowTrackSharing1Data.at(c)->SetMarkerColor(kBlack);
       narrowTrackSharing2Data.at(c)->SetMarkerStyle(kFullCircle);
       narrowTrackSharing2Data.at(c)->SetMarkerColor(kBlack);

       narrowTrackSharing1Data.at(c)->SetTitle("");
       narrowTrackSharing2Data.at(c)->SetTitle("");

       narrowTrackSharing1Data.at(c)->GetXaxis()->SetTitle("#eta #pm 1");
       narrowTrackSharing2Data.at(c)->GetXaxis()->SetTitle("#eta #pm 2");

       narrowTrackSharing1Data.at(c)->SetMaximum(2.0*  narrowTrackSharing1Data.at(c)->GetMaximum());
       narrowTrackSharing2Data.at(c)->SetMaximum(1.5*  narrowTrackSharing2Data.at(c)->GetMaximum());


       narrowTrackSharing1DataCan.at(c)->cd();
       narrowTrackSharing1Data.at(c)->DrawNormalized(""); 
       narrowTrackSharing1DataCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ narrowTrackSharing1Data.at(c)->GetName()+ ".root").c_str());
       narrowTrackSharing1DataCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ narrowTrackSharing1Data.at(c)->GetName()+ ".eps").c_str());


       narrowTrackSharing2DataCan.at(c)->cd();
       narrowTrackSharing2Data.at(c)->DrawNormalized(""); 
       narrowTrackSharing2DataCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ narrowTrackSharing2Data.at(c)->GetName() + ".root").c_str());
       narrowTrackSharing2DataCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ narrowTrackSharing2Data.at(c)->GetName()+ ".eps").c_str());
       //c37.SaveAs(("output/"+(string)dir+"/"+  "narrowTrackSharing1.eps").c_str());
       //
       //
       
       
       widthTopInOutCan.at(c)->cd();
       widthTopInOut.at(c)->Draw("");
       widthTopInOutCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthTopInOut.at(c)->GetName()+ ".eps").c_str());
       widthTopInOutCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthTopInOut.at(c)->GetName()+ ".root").c_str());

       widthTopInOutOnlyCan.at(c)->cd();
       widthTopInOutOnly.at(c)->Draw("");
       widthTopInOutOnlyCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthTopInOutOnly.at(c)->GetName()+ ".eps").c_str());
       widthTopInOutOnlyCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthTopInOutOnly.at(c)->GetName()+ ".root").c_str());

       widthTopOutInCan.at(c)->cd();
       widthTopOutIn.at(c)->Draw("");
       widthTopOutInCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthTopOutIn.at(c)->GetName()+ ".eps").c_str());
       widthTopOutInCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthTopOutIn.at(c)->GetName()+ ".root").c_str());

       widthBottomInOutCan.at(c)->cd();
       widthBottomInOut.at(c)->Draw("");
       widthBottomInOutCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthBottomInOut.at(c)->GetName()+ ".eps").c_str());
       widthBottomInOutCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthBottomInOut.at(c)->GetName()+ ".root").c_str());

       widthBottomInOutOnlyCan.at(c)->cd();
       widthBottomInOutOnly.at(c)->Draw("");
       widthBottomInOutOnlyCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthBottomInOutOnly.at(c)->GetName()+ ".eps").c_str());
       widthBottomInOutOnlyCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthBottomInOutOnly.at(c)->GetName()+ ".root").c_str());

       widthBottomOutInCan.at(c)->cd();
       widthBottomOutIn.at(c)->Draw("");
       widthBottomOutInCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthBottomOutIn.at(c)->GetName()+ ".eps").c_str());
       widthBottomOutInCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ widthBottomOutIn.at(c)->GetName()+ ".root").c_str());

       chargeTopInOutCan.at(c)->cd();
       chargeTopInOut.at(c)->Draw("");
       chargeTopInOutCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeTopInOut.at(c)->GetName()+ ".eps").c_str());
       chargeTopInOutCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeTopInOut.at(c)->GetName()+ ".root").c_str());

       chargeTopInOutOnlyCan.at(c)->cd();
       chargeTopInOutOnly.at(c)->Draw("");
       chargeTopInOutOnlyCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeTopInOutOnly.at(c)->GetName()+ ".eps").c_str());
       chargeTopInOutOnlyCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeTopInOutOnly.at(c)->GetName()+ ".root").c_str());

       chargeTopOutInCan.at(c)->cd();
       chargeTopOutIn.at(c)->Draw("");
       chargeTopOutInCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeTopOutIn.at(c)->GetName()+ ".eps").c_str());
       chargeTopOutInCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeTopOutIn.at(c)->GetName()+ ".root").c_str());

       chargeBottomInOutCan.at(c)->cd();
       chargeBottomInOut.at(c)->Draw("");
       chargeBottomInOutCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeBottomInOut.at(c)->GetName()+ ".eps").c_str());
       chargeBottomInOutCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeBottomInOut.at(c)->GetName()+ ".root").c_str());

       chargeBottomInOutOnlyCan.at(c)->cd();
       chargeBottomInOutOnly.at(c)->Draw("");
       chargeBottomInOutOnlyCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeBottomInOutOnly.at(c)->GetName()+ ".eps").c_str());
       chargeBottomInOutOnlyCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeBottomInOutOnly.at(c)->GetName()+ ".root").c_str());

       chargeBottomOutInCan.at(c)->cd();
       chargeBottomOutIn.at(c)->Draw("");
       chargeBottomOutInCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeBottomOutIn.at(c)->GetName()+ ".eps").c_str());
       chargeBottomOutInCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ chargeBottomOutIn.at(c)->GetName()+ ".root").c_str());


       timingErrorTopCan.at(c)->cd();
       timingErrorTop.at(c)->DrawNormalized(""); 
       timingErrorTopInOut.at(c)->SetLineColor(kRed); 
       timingErrorTopInOut.at(c)->DrawNormalized("same hist"); 
       timingErrorTopCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ timingErrorTop.at(c)->GetName()+ ".root").c_str());
       timingErrorTopCan.at(c)->SaveAs(("output/"+(string)dir+"/"+ timingErrorTop.at(c)->GetName()+ ".eps").c_str());

       }


      /////////////////////////////////////////////////////



       TCanvas c35("dummyCTplot","dummyCTplot");
       dummyCTplot->Draw("hist"); 
       c35.SaveAs(("output/"+(string)dir+"/"+  "dummyCTplot.root").c_str());
       c35.SaveAs(("output/"+(string)dir+"/"+  "dummyCTplot.eps").c_str());

    //cout << "in here 7"  << endl;

   return 0;
}
