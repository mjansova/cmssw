//./triggerEfficiency test_shallowTrackAndClusterNarrowInfoCRUZETforTriggerEff.root triggerEff 
//./triggerEfficiency test_shallowTrackCRUZET_2017VR.root  triggerEffNewVR 
//./triggerEfficiency test_shallowTrackCRUZET_2017VRNewClusteringTimingTest.root  triggerEffNewVR 
//./triggerEfficiency test_shallowTrackCRUZET_2017VRNewClusteringTimingTestNoQualityCuts.root  triggerEffNewVR 

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


       vector<float>* partition = 0;
       vector<float>* clustercharge = 0;
       vector<float>* clusterwidth = 0;
       vector<float>* clusterlayerwheel = 0;
       vector<float>* clusterstripChargeSubdetid = 0;
       vector<float>* clusterstripCharge = 0;
       vector<float>* clusterstripChargeLayerwheel = 0;
       vector<float>* clusterstripChargeLayerwheel2 = 0;
       vector<float>* clusterstripChargeStripNr = 0;
       vector<float>* clusterstripChargeStripNr2 = 0;
       vector<float>* clusterstripChargeTotWidth = 0;
       vector<float>* clusterstripChargeTotWidth2 = 0;
       vector<float>* clusterstripChargeTotCharge = 0;
       vector<float>* clusterstripChargeTotCharge2 = 0;
       vector<float>* clusterstripChargeLocalTrackPhi = 0;
       vector<float>* clusterstripChargeLocalTrackPhi2 = 0;
       vector<float>* clusterstripChargeGlobalTrackPhi = 0;
       vector<float>* clusterstripChargeGlobalTrackPhi2 = 0;
       vector<float>* clusterstripChargeLocalTrackTheta = 0;
       vector<float>* clusterstripChargeLocalTrackTheta2 = 0;
       vector<float>* clusterstripChargeGlobalTrackTheta = 0;
       vector<float>* clusterstripChargeGlobalTrackTheta2 = 0;
       vector<unsigned>* clusterstripChargeDetid = 0;
       vector<unsigned>* clusterstripChargeDetid2 = 0;
       vector<float>* clusterstripChargeLocalX = 0;
       vector<float>* clusterstripChargeLocalX2 = 0;
       vector<float>* clusterstripChargeLocalY = 0;
       vector<float>* clusterstripChargeLocalY2 = 0;
       vector<float>* tsostrackPt = 0;
       vector<float>* tsostrackPt2 = 0;
       vector<float>* clusterstripChargetrackPt = 0;
       vector<float>* clusterstripChargetrackPt2 = 0;
       vector<float>* tsoslocalphi = 0;
       vector<float>* tsoslocalphi2 = 0;
       vector<float>* tsoslocaltheta = 0;
       vector<float>* tsoslocaltheta2 = 0;
       vector<float>* tsoslocalpitch = 0;
       vector<float>* tsoslocalpitch2 = 0;
       vector<float>* clustersensorThickness = 0;
       vector<float>* clustersensorThickness2 = 0;
       vector<float>* tsosBdotY = 0;
       vector<float>* tsosBdotY2 = 0;
       vector<float>* clusterstripChargelocalpitch = 0;
       vector<float>* clusterstripChargelocalpitch2  = 0;
       vector<float>* clusterstripChargesensorThickness = 0;
       vector<float>* clusterstripChargesensorThickness2  = 0;
       vector<float>* clusterstripChargeBdotY = 0;
       vector<float>* clusterstripChargeBdotY2  = 0;
       unsigned int  nroftracks;
       unsigned int  nrofevents;
       bool passHLTL1SingleMuOpenDTv2;
       bool passHLTL1SingleMuOpenv2;

       vector<float> subpartition;
       vector<float> subpartition2;
       vector<float> subclustercharge;
       vector<float> subclustercharge2;
       vector<float> subclusterwidth;
       vector<float> subclusterwidth2;
       vector<float> subclusterlayerwheel;
       vector<float> subclusterlayerwheel2;
       vector<float> subclusterstripChargeSubdetid;
       vector<float> subclusterstripChargeSubdetid2;
       vector<float> subclusterstripCharge;
       vector<float> subclusterstripCharge2;
       vector<float> subclusterstripChargeLayerwheel;
       vector<float> subclusterstripChargeLayerwheel2;
       vector<float> subclusterstripChargeStripNr;
       vector<float> subclusterstripChargeStripNr2;
       vector<float> subclusterstripChargeTotWidth;
       vector<float> subclusterstripChargeTotWidth2;
       vector<float> subclusterstripChargeTotCharge;
       vector<float> subclusterstripChargeTotCharge2;
       vector<float> subclusterstripChargeLocalTrackPhi;
       vector<float> subclusterstripChargeLocalTrackPhi2;
       vector<float> subclusterstripChargeGlobalTrackPhi;
       vector<float> subclusterstripChargeGlobalTrackPhi2;
       vector<float> subclusterstripChargeLocalTrackTheta;
       vector<float> subclusterstripChargeLocalTrackTheta2;
       vector<float> subclusterstripChargeGlobalTrackTheta;
       vector<float> subclusterstripChargeGlobalTrackTheta2;
       vector<unsigned> subclusterstripChargeDetid;
       vector<unsigned> subclusterstripChargeDetid2;
       vector<float> subclusterstripChargeLocalX;
       vector<float> subclusterstripChargeLocalX2;
       vector<float> subclusterstripChargeLocalY;
       vector<float> subclusterstripChargeLocalY2;
       vector<float> subtsostrackPt;
       vector<float> subtsostrackPt2;
       vector<float> subclusterstripChargetrackPt;
       vector<float> subclusterstripChargetrackPt2;
       vector<float> subtsoslocalphi;
       vector<float> subtsoslocalphi2;
       vector<float> subtsoslocaltheta;
       vector<float> subtsoslocaltheta2;
       vector<float> subtsoslocalpitch;
       vector<float> subtsoslocalpitch2;
       vector<float> subclustersensorThickness;
       vector<float> subclustersensorThickness2;
       vector<float> subtsosBdotY;
       vector<float> subtsosBdotY2;
       vector<float> subclusterstripChargelocalpitch;
       vector<float> subclusterstripChargelocalpitch2;
       vector<float> subclusterstripChargesensorThickness;
       vector<float> subclusterstripChargesensorThickness2;
       vector<float> subclusterstripChargeBdotY;
       vector<float> subclusterstripChargeBdotY2;
       vector<unsigned int>  subnroftracks;
       vector<unsigned int>  subnrofevents;
       vector<bool>  subpassHLTL1SingleMuOpenDTv2;
       vector<bool>  subpassHLTL1SingleMuOpenv2;
       vector<bool>  evsubpassHLTL1SingleMuOpenDTv2;
       vector<bool>  evsubpassHLTL1SingleMuOpenv2;


       t1->SetBranchAddress("clustersubdetid",  &partition );
       t1->SetBranchAddress("clustercharge",  &clustercharge );
       t1->SetBranchAddress("clusterwidth",  &clusterwidth );
       t1->SetBranchAddress("clusterlayerwheel",  &clusterlayerwheel );
       t1->SetBranchAddress("clusterstripChargeSubdetid",  &clusterstripChargeSubdetid );
       t1->SetBranchAddress("clusterstripCharge",  &clusterstripCharge );
       t1->SetBranchAddress("clusterstripChargeLayerwheel",  &clusterstripChargeLayerwheel );
       t1->SetBranchAddress("clusterstripChargeStripNr",  &clusterstripChargeStripNr );
       t1->SetBranchAddress("clusterstripChargeTotWidth",  &clusterstripChargeTotWidth );
       //t1->SetBranchAddress("clusterstripChargeTotCharge",  &clusterstripChargeTotWidth ); //@MJ@ TODO just for now -> switched values
       t1->SetBranchAddress("clusterstripChargeTotCharge",  &clusterstripChargeTotCharge );
       t1->SetBranchAddress("clusterstripChargeLocalTrackPhi",  &clusterstripChargeLocalTrackPhi );
       t1->SetBranchAddress("clusterstripChargeGlobalTrackPhi",  &clusterstripChargeGlobalTrackPhi );
       t1->SetBranchAddress("clusterstripChargeLocalTrackTheta",  &clusterstripChargeLocalTrackTheta );
       t1->SetBranchAddress("clusterstripChargeGlobalTrackTheta",  &clusterstripChargeGlobalTrackTheta );
       t1->SetBranchAddress("clusterstripChargeDetid",  &clusterstripChargeDetid );
       t1->SetBranchAddress("clusterstripChargeLocalX",  &clusterstripChargeLocalX );
       t1->SetBranchAddress("clusterstripChargeLocalY",  &clusterstripChargeLocalY );
       t1->SetBranchAddress("clusterstripChargetrackPt",  &clusterstripChargetrackPt );
       t1->SetBranchAddress("tsostrackPt",  &tsostrackPt );
       t1->SetBranchAddress("tsoslocalphi",  &tsoslocalphi );
       t1->SetBranchAddress("tsoslocaltheta",  &tsoslocaltheta );
       t1->SetBranchAddress("tsoslocalpitch",  &tsoslocalpitch );
       t1->SetBranchAddress("clustersensorThickness",  &clustersensorThickness );
       t1->SetBranchAddress("tsosBdotY",  &tsosBdotY );
       t1->SetBranchAddress("clusterstripChargelocalpitch",  &clusterstripChargelocalpitch );
       t1->SetBranchAddress("clusterstripChargesensorThickness",  &clusterstripChargesensorThickness );
       t1->SetBranchAddress("clusterstripChargeBdotY",  &clusterstripChargeBdotY );
       t1->SetBranchAddress("nroftracks",  &nroftracks );
       t1->SetBranchAddress("nrofevents",  &nrofevents );
       t1->SetBranchAddress("passHLTL1SingleMuOpenDTv2",  &passHLTL1SingleMuOpenDTv2 );
       t1->SetBranchAddress("passHLTL1SingleMuOpenv2",  &passHLTL1SingleMuOpenv2 );

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
           uint32_t up = partition->size();
           subnroftracks.push_back(nroftracks);
           subnrofevents.push_back(nrofevents);
           evsubpassHLTL1SingleMuOpenDTv2.push_back(passHLTL1SingleMuOpenDTv2);
           evsubpassHLTL1SingleMuOpenv2.push_back(passHLTL1SingleMuOpenv2);
           //cout << "passHLTL1SingleMuOpenDTv2 is when filling " << passHLTL1SingleMuOpenDTv2 << endl;
           //cout << "ntracks is when filling " << nroftracks << endl;
           for(uint32_t k=0; k<up;k++)
           {
               if( true )
               {
                   if( true ) //@MJ@ TODO trigger
                   {
                       subpartition.push_back(partition->at(k));
                       subclustercharge.push_back(clustercharge->at(k));
                       subclusterwidth.push_back(clusterwidth->at(k));
                       subclusterlayerwheel.push_back(clusterlayerwheel->at(k));
                       subtsostrackPt.push_back(tsostrackPt->at(k));
                       subtsoslocalphi.push_back(tsoslocalphi->at(k));
                       subtsoslocaltheta.push_back(tsoslocaltheta->at(k));
                       subtsoslocalpitch.push_back(tsoslocalpitch->at(k));
                       subclustersensorThickness.push_back(clustersensorThickness->at(k));
                       subtsosBdotY.push_back(tsosBdotY->at(k));
                       //cout << "filling variable with " << passHLTL1SingleMuOpenDTv2 << endl;
                       subpassHLTL1SingleMuOpenDTv2.push_back(passHLTL1SingleMuOpenDTv2);
                       subpassHLTL1SingleMuOpenv2.push_back(passHLTL1SingleMuOpenv2);
                   }
               }
           }
           //perstrip
           uint32_t upStrip = clusterstripChargeSubdetid->size();
           for(uint32_t k=0; k<upStrip;k++)
           {
               if( true)
               {
                   if(true ) //@MJ@ TODO trigger
                   {
                       subclusterstripChargeSubdetid.push_back(clusterstripChargeSubdetid->at(k));
                       subclusterstripCharge.push_back(clusterstripCharge->at(k));
                       subclusterstripChargeLayerwheel.push_back(clusterstripChargeLayerwheel->at(k));
                       subclusterstripChargeStripNr.push_back(clusterstripChargeStripNr->at(k));
                       subclusterstripChargeTotWidth.push_back(clusterstripChargeTotWidth->at(k));
                       subclusterstripChargeTotCharge.push_back(clusterstripChargeTotCharge->at(k));
                       subclusterstripChargeLocalTrackPhi.push_back(clusterstripChargeLocalTrackPhi->at(k));
                       subclusterstripChargeGlobalTrackPhi.push_back(clusterstripChargeGlobalTrackPhi->at(k));
                       subclusterstripChargeLocalTrackTheta.push_back(clusterstripChargeLocalTrackTheta->at(k));
                       subclusterstripChargeGlobalTrackTheta.push_back(clusterstripChargeGlobalTrackTheta->at(k));
                       subclusterstripChargeDetid.push_back(clusterstripChargeDetid->at(k));
                       subclusterstripChargeLocalX.push_back(clusterstripChargeLocalX->at(k));
                       subclusterstripChargeLocalY.push_back(clusterstripChargeLocalY->at(k));
                       subclusterstripChargetrackPt.push_back(clusterstripChargetrackPt->at(k));
                       subclusterstripChargelocalpitch.push_back(clusterstripChargelocalpitch->at(k));
                       subclusterstripChargesensorThickness.push_back(clusterstripChargesensorThickness->at(k));
                       subclusterstripChargeBdotY.push_back(clusterstripChargeBdotY->at(k));
                   }
               }
           }

   }


   //}maybe




       TH1F* narrowTrackClusterWidthData = new TH1F("narrowTrackClusterWidthData","narrowTrackClusterWidthData", 10, 0, 10);         
       TH1F* usefulCluster = new TH1F("usefulCluster","usefulCluster", 20, 0, 20);         
       TH1F* usefulClustersmuOpenDT = new TH1F("usefulClustersmuOpenDT","usefulClustersmuOpenDT", 20, 0, 20);         
       TH1F* usefulClustersmuOpen = new TH1F("usefulClustersmuOpen","usefulClustersmuOpen", 20, 0, 20);         
       TH1F* allClustersmuOpenDT = new TH1F("usefulClustersmuOpenDT","usefulClustersmuOpenDT", 20, 0, 20);         
       TH1F* allClustersmuOpen = new TH1F("usefulClustersmuOpen","usefulClustersmuOpen", 20, 0, 20);         
       TH1F* allCluster = new TH1F("allCluster","allCluster", 20, 0, 20);         

       uint8_t TIDoffset = 4;
       uint8_t TECoffset = 7;
       bool narrowtrack = false;
    //cout << "in here 2"  << endl;
       uint32_t trackCounter = 0;
       uint32_t trackCounterpassHLTL1SingleMuOpenDTv2=0;
       uint32_t trackCounterpassHLTL1SingleMuOpenv2=0;
       uint32_t evCounterpassHLTL1SingleMuOpenDTv2=0;
       uint32_t evCounterpassHLTL1SingleMuOpenv2=0;

       for(uint32_t t=0; t<subnroftracks.size(); t++)
       {
           trackCounter+=subnroftracks.at(t);
           //cout << " evsubpassHLTL1SingleMuOpenDTv2 " << evsubpassHLTL1SingleMuOpenDTv2.at(t) << endl;
           if(evsubpassHLTL1SingleMuOpenDTv2.at(t))
           {
               trackCounterpassHLTL1SingleMuOpenDTv2+=subnroftracks.at(t);
               evCounterpassHLTL1SingleMuOpenDTv2++;
           }
           if(evsubpassHLTL1SingleMuOpenDTv2.at(t))
           {
               trackCounterpassHLTL1SingleMuOpenv2+=subnroftracks.at(t);
               evCounterpassHLTL1SingleMuOpenv2++;
           }
       }
       cout << "we need " << (float) subnrofevents.size()/trackCounter << "events triggered by cosmics trigger to get one track" << endl;
       cout << "we have " << (float) trackCounter << "tracks" << endl;
       cout << "we have " << (float) subnrofevents.size() << "events" << endl;
       cout << "we need " << (float) evCounterpassHLTL1SingleMuOpenDTv2/trackCounterpassHLTL1SingleMuOpenDTv2 << "events triggered by trackCounterpassHLTL1SingleMuOpenDTv2  tigger to get one track" << endl;
       cout << "we need " << (float) evCounterpassHLTL1SingleMuOpenv2/trackCounterpassHLTL1SingleMuOpenv2 << "events triggered by trackCounterpassHLTL1SingleMuOpenv2 trigger to get one track" << endl;

       bool smuOpenDT = false;
       bool smuOpen = false;

       vector<string> labels {"IB1", "IB2", "OB1", "OB2", "W1a", "W2a", "W3a", "W1b", "W2b", "W3b", "W4b", "W5b", "W6b", "W7b"};
       for(uint32_t m = 0; m<subclustercharge.size(); m++)
       {
           narrowtrack=false;
           if( abs(tan(subtsoslocaltheta.at(m))*cos(subtsoslocalphi.at(m))) < narrowness*(subtsoslocalpitch.at(m)/subclustersensorThickness.at(m)) )
           {
               narrowTrackClusterWidthData->Fill(subclusterwidth.at(m));
               narrowtrack=true;
           }
           smuOpenDT = subpassHLTL1SingleMuOpenDTv2.at(m);
           smuOpen = subpassHLTL1SingleMuOpenv2.at(m);
           //cout << "smuOpenDT is " << subpassHLTL1SingleMuOpenDTv2.at(m) << endl;

               if(subpartition.at(m) == 3 ) //TIB
               {
                   if(subclusterlayerwheel.at(m)>0 && subclusterlayerwheel.at(m)<3)
                   {
                       allCluster->Fill(1);
                       if(smuOpenDT)
                           allClustersmuOpenDT->Fill(1);
                       if(smuOpen)
                           allClustersmuOpen->Fill(1);
                       if(narrowtrack)
                       {
                           usefulCluster->Fill(1);
                           if(smuOpenDT)
                               usefulClustersmuOpenDT->Fill(1);
                           if(smuOpen)
                               usefulClustersmuOpen->Fill(1);
                       }
                   }
                   else
                   {
                       allCluster->Fill(2);
                       if(smuOpenDT)
                           allClustersmuOpenDT->Fill(2);
                       if(smuOpen)
                           allClustersmuOpen->Fill(2);
                       if(narrowtrack)
                       {
                           usefulCluster->Fill(2); 
                           if(smuOpenDT)
                               usefulClustersmuOpenDT->Fill(2);
                           if(smuOpen)
                               usefulClustersmuOpen->Fill(2);
                       }          

                   }

               }
               if(subpartition.at(m) == 5 ) //TOB
               {
                   if(subclusterlayerwheel.at(m)>0 && subclusterlayerwheel.at(m)<5)
                   {
                       allCluster->Fill(3);
                       if(smuOpenDT)
                           allClustersmuOpenDT->Fill(3);
                       if(smuOpen)
                           allClustersmuOpen->Fill(3);
                       if(narrowtrack)
                       {
                           usefulCluster->Fill(3);   
                           if(smuOpenDT)
                               usefulClustersmuOpenDT->Fill(3);
                           if(smuOpen)
                               usefulClustersmuOpen->Fill(3);
                       }        
                   }
                   else
                   {
                       allCluster->Fill(4);
                       if(smuOpenDT)
                           allClustersmuOpenDT->Fill(4);
                       if(smuOpen)
                           allClustersmuOpen->Fill(4);
                       if(narrowtrack)
                       {
                           usefulCluster->Fill(4);
                           if(smuOpenDT)
                               usefulClustersmuOpenDT->Fill(4);
                           if(smuOpen)
                               usefulClustersmuOpen->Fill(4);
                       }       
                   }

               }
               if(subpartition.at(m) == 4 ) //TID
               {
                       cout << "TID filling bin " << TIDoffset + subclusterlayerwheel.at(m) << endl;
                       allCluster->Fill(TIDoffset + subclusterlayerwheel.at(m) );           
                       if(smuOpenDT)
                           allClustersmuOpenDT->Fill( TIDoffset + subclusterlayerwheel.at(m));
                       if(smuOpen)
                           allClustersmuOpen->Fill(TIDoffset + subclusterlayerwheel.at(m));
                       if(narrowtrack)
                       {
                           usefulCluster->Fill(TIDoffset + subclusterlayerwheel.at(m) );           
                           if(smuOpenDT)
                               usefulClustersmuOpenDT->Fill(TIDoffset + subclusterlayerwheel.at(m));
                           if(smuOpen)
                               usefulClustersmuOpen->Fill(TIDoffset + subclusterlayerwheel.at(m));
                       }
               }
               if(subpartition.at(m) == 6 ) //TEC
               {
                       cout << "TEC filling bin " << TECoffset + subclusterlayerwheel.at(m) << " for layer " << subclusterlayerwheel.at(m) << endl;
                       allCluster->Fill(TECoffset + subclusterlayerwheel.at(m) );           
                       if(smuOpenDT)
                           allClustersmuOpenDT->Fill( TECoffset + subclusterlayerwheel.at(m));
                       if(smuOpen)
                           allClustersmuOpen->Fill(TECoffset + subclusterlayerwheel.at(m));
                       if(narrowtrack)
                       {
                           usefulCluster->Fill(TECoffset + subclusterlayerwheel.at(m) );           
                           if(smuOpenDT)
                               usefulClustersmuOpenDT->Fill(TECoffset + subclusterlayerwheel.at(m));
                           if(smuOpen)
                               usefulClustersmuOpen->Fill(TECoffset + subclusterlayerwheel.at(m));
                       }
               }
               
         
       }


       //F
       TH1F* fofCluster = (TH1F*) usefulCluster->Clone();
       TH1F* fofClustersmuOpenDT = (TH1F*) usefulClustersmuOpenDT->Clone();
       TH1F* fofClustersmuOpen = (TH1F*) usefulClustersmuOpen->Clone();
       fofCluster->Divide(allCluster);           
       fofClustersmuOpenDT->Divide(allClustersmuOpenDT);
       fofClustersmuOpen->Divide(allClustersmuOpen);

       //clusters per track
       TH1F*  clustersPerTrackClustersmuOpenDT = (TH1F*) usefulClustersmuOpenDT->Clone();
       clustersPerTrackClustersmuOpenDT->Scale(1.0/trackCounterpassHLTL1SingleMuOpenDTv2);

       vector<float> stripChargeVec;
       TH1F* narrowTrackSharing1Data = new TH1F("narrowTrackSharing1Data", "narrowTrackSharing1Data" , 100, 0, 1 );
       TH1F* narrowTrackSharing2Data = new TH1F("narrowTrackSharing2Data", "narrowTrackSharing2Data" , 100, 0, 1 );

       for(uint32_t m = 0; m<subclusterstripCharge.size(); m++)
       {
           
           if( abs(tan(subclusterstripChargeLocalTrackTheta.at(m))*cos(subclusterstripChargeLocalTrackPhi.at(m))) < narrowness*(subclusterstripChargelocalpitch.at(m)/subclusterstripChargesensorThickness.at(m)) )
           {
               stripChargeVec.push_back(subclusterstripCharge.at(m));
               if(subclusterstripChargeTotWidth.at(m) == stripChargeVec.size())
               {
                   float maxVal = 0;
                   int indexVal = -1;
                   for( uint32_t s =0; s<stripChargeVec.size(); s++)
                   {
                       if(stripChargeVec.at(s) > maxVal)
                       {
                           maxVal = stripChargeVec.at(s);
                           indexVal = s;
                       }

                   }
                   if(indexVal-1>=0)
                   {
                       narrowTrackSharing1Data->Fill((float)stripChargeVec.at(indexVal-1)/maxVal);
                   }
                   if(indexVal+1<stripChargeVec.size())
                       narrowTrackSharing1Data->Fill((float)stripChargeVec.at(indexVal+1)/maxVal);
                   if(indexVal-2>=0)
                       narrowTrackSharing2Data->Fill((float)stripChargeVec.at(indexVal-2)/maxVal);
                   if(indexVal+2<stripChargeVec.size())
                       narrowTrackSharing2Data->Fill((float)stripChargeVec.at(indexVal+2)/maxVal);

                   stripChargeVec.clear();
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
       narrowTrackClusterWidthData->SetMarkerStyle(kFullCircle); 
       narrowTrackClusterWidthData->SetMarkerColor(kBlack); 
       narrowTrackSharing1Data->SetMarkerStyle(kFullCircle);
       narrowTrackSharing1Data->SetMarkerColor(kBlack);
       narrowTrackSharing2Data->SetMarkerStyle(kFullCircle);
       narrowTrackSharing2Data->SetMarkerColor(kBlack);

       narrowTrackClusterWidthData->SetTitle(""); 
       narrowTrackSharing1Data->SetTitle("");
       narrowTrackSharing2Data->SetTitle("");

       narrowTrackSharing1Data->GetXaxis()->SetTitle("#eta #pm 1");
       narrowTrackSharing2Data->GetXaxis()->SetTitle("#eta #pm 2");

       narrowTrackClusterWidthData->SetMaximum(1.5*  narrowTrackClusterWidthData->GetMaximum());
       narrowTrackSharing1Data->SetMaximum(2.0*  narrowTrackSharing1Data->GetMaximum());
       narrowTrackSharing2Data->SetMaximum(1.5*  narrowTrackSharing2Data->GetMaximum());


       for(uint32_t l=0; l<labels.size(); l++)
       {

           usefulCluster->GetXaxis()->SetBinLabel(l+2,labels.at(l).c_str()); 
           fofCluster->GetXaxis()->SetBinLabel(l+2,labels.at(l).c_str()); 
       }

       TCanvas c36("narrowTrackClusterWidth","narrowTrackClusterWidth");
       narrowTrackClusterWidthData->DrawNormalized(""); 
       c36.SaveAs(("output/"+(string)dir+"/"+  "narrowTrackClusterWidth.root").c_str());
       c36.SaveAs(("output/"+(string)dir+"/"+  "narrowTrackClusterWidth.eps").c_str());

       TCanvas c37("narrowTrackSharing1","narrowTrackSharing1");
       narrowTrackSharing1Data->DrawNormalized(""); 
       c37.SaveAs(("output/"+(string)dir+"/"+  "narrowTrackSharing1.root").c_str());
       c37.SaveAs(("output/"+(string)dir+"/"+  "narrowTrackSharing1.eps").c_str());

       TCanvas c38("narrowTrackSharing2","narrowTrackSharing2");
       narrowTrackSharing2Data->DrawNormalized(""); 
       c38.SaveAs(("output/"+(string)dir+"/"+  "narrowTrackSharing2.root").c_str());
       c38.SaveAs(("output/"+(string)dir+"/"+  "narrowTrackSharing2.eps").c_str());


       TCanvas c39("usefulCluster","usefulCluster");
       usefulCluster->SetTitle(""); 
       usefulCluster->Draw("hist"); 
       c39.SaveAs(("output/"+(string)dir+"/"+  "usefulCluster.root").c_str());
       c39.SaveAs(("output/"+(string)dir+"/"+  "usefulCluster.eps").c_str());

       TCanvas c40("allCluster","allCluster");
       allCluster->Draw("hist"); 
       c40.SaveAs(("output/"+(string)dir+"/"+  "allCluster.root").c_str());
       c40.SaveAs(("output/"+(string)dir+"/"+  "allCluster.eps").c_str());

       TCanvas c41("fofCluster","fofCluster");
       fofCluster->SetTitle(""); 
       fofCluster->Draw("hist"); 
       c41.SaveAs(("output/"+(string)dir+"/"+  "fofCluster.root").c_str());
       c41.SaveAs(("output/"+(string)dir+"/"+  "fofCluster.eps").c_str());

       TCanvas c42("fofClustersmuOpenDT","fofClustersmuOpenDT");
       fofClustersmuOpenDT->Draw("hist"); 
       c42.SaveAs(("output/"+(string)dir+"/"+  "fofClustersmuOpenDT.root").c_str());
       c42.SaveAs(("output/"+(string)dir+"/"+  "fofClustersmuOpenDT.eps").c_str());

       TCanvas c43("fofClustersmuOpen","fofClustersmuOpen");
       fofClustersmuOpen->Draw("hist"); 
       c43.SaveAs(("output/"+(string)dir+"/"+  "fofClustersmuOpen.root").c_str());
       c43.SaveAs(("output/"+(string)dir+"/"+  "fofClustersmuOpen.eps").c_str());

       TCanvas c44("clustersPerTrackClustersmuOpenDT","clustersPerTrackClustersmuOpenDT");
       clustersPerTrackClustersmuOpenDT->Draw("hist"); 
       c44.SaveAs(("output/"+(string)dir+"/"+  "clustersPerTrackClustersmuOpenDT.root").c_str());
       c44.SaveAs(("output/"+(string)dir+"/"+  "clustersPerTrackClustersmuOpenDT.eps").c_str());

      /////////////////////////////////////////////////////



       TCanvas c35("dummyCTplot","dummyCTplot");
       dummyCTplot->Draw("hist"); 
       c35.SaveAs(("output/"+(string)dir+"/"+  "dummyCTplot.root").c_str());
       c35.SaveAs(("output/"+(string)dir+"/"+  "dummyCTplot.eps").c_str());

    //cout << "in here 7"  << endl;

   return 0;
}
