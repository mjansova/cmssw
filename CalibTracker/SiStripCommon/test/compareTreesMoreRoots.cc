//./compareTreesMoreRoots 2 onTrackClustercosmicsData0T.root onTrackClustercosmicsMC0T.root TEC yes test > logsMean/log0TdataMcTEC
//./compareTreesMoreRoots 3 onTrackClustercosmicsData0T.root onTrackClustercosmicsMC0T.root onTrackClustercosmicsData38T.root TEC yes test > logsMean/log0TdataMcTEC
//
//TEMPERATURE
//./compareTreesMoreRoots 6 test_shallowTrackAndClusterFullInfoNoPU2016Data.root test_shallowTrackAndClusterFullInfo.rootTemperature_min20_218.4step3.root  test_shallowTrackAndClusterFullInfo.rootTemperature_min10_245.7step3.root  test_shallowTrackAndClusterFullInfo.rootTemperature_def_273step3.root test_shallowTrackAndClusterFullInfo.rootTemperature_pl10_300.3step3.root  test_shallowTrackAndClusterFullInfo.rootTemperature_pl20_327.6step3.root TOB yes lowPU2016Temperature > logsMean/loglowPU2016TemperatureTOB
//
//DEPLETION VOLTAGE
//./compareTreesMoreRoots 6 test_shallowTrackAndClusterFullInfoNoPU2016Data.root test_shallowTrackAndClusterFullInfo.rootDepletionVoltage_min20_136step3.root  test_shallowTrackAndClusterFullInfo.rootDepletionVoltage_min10_153step3.root  test_shallowTrackAndClusterFullInfo.rootDepletionVoltage_def_170step3.root test_shallowTrackAndClusterFullInfo.rootDepletionVoltage_pl10_187step3.root  test_shallowTrackAndClusterFullInfo.rootDepletionVoltage_pl20_204step3.root TOB yes lowPU2016DepletionVoltage > logsMean/loglowPU2016DepletionVoltageTOB
//
//APPLIED VOLTAGE
//./compareTreesMoreRoots 6 test_shallowTrackAndClusterFullInfoNoPU2016Data.root test_shallowTrackAndClusterFullInfo.rootAppliedVoltage_min20_240step3.root  test_shallowTrackAndClusterFullInfo.rootAppliedVoltage_min10_270step3.root  test_shallowTrackAndClusterFullInfo.rootAppliedVoltage_def_300step3.root test_shallowTrackAndClusterFullInfo.rootAppliedVoltage_pl10_330step3.root  test_shallowTrackAndClusterFullInfo.rootAppliedVoltage_pl20_360step3.root TOB yes lowPU2016AppliedVoltage > logsMean/loglowPU2016AppliedVoltageTOB
//
//CHARGE MOBILTY
//./compareTreesMoreRoots 6 test_shallowTrackAndClusterFullInfoNoPU2016Data.root test_shallowTrackAndClusterFullInfo.rootChargeMobility_min20_248step3.root  test_shallowTrackAndClusterFullInfo.rootChargeMobility_min10_279step3.root  test_shallowTrackAndClusterFullInfo.rootChargeMobility_def_310step3.root test_shallowTrackAndClusterFullInfo.rootChargeMobility_pl10_341step3.root  test_shallowTrackAndClusterFullInfo.rootChargeMobility_pl20_372step3.root TOB yes lowPU2016ChargeMobility > logsMean/loglowPU2016ChargeMobilityTOB
//
//CROSS TALK
//./compareTreesMoreRoots 4 test_shallowTrackAndClusterFullInfoNoPU2016DatastripProp.root test_shallowTrackAndClusterFullInfoNoPU2016MCstripProp.root test_shallowTrackAndClusterFullInfo.rootCT_1step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.4v2step3.root   TOB yes lowPU2016CrossTalk > logsMean/loglowPU2016CrossTalkTOB

//CROSTALK tuning l<5
//./compareTreesMoreRoots 12 test_shallowTrackAndClusterFullInfoNoPU2016DatastripProp.root  test_shallowTrackAndClusterFullInfo.rootCT_0.70step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.72step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.74step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.76step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.78step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.80step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.82step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.84step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.86step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.88step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.90step3.root TOB yes lowPU2016CrossTalkTOBL1tuning > logsMean/lowPU2016CrossTalkTOBL1tuningTOB


//CROSTALK tuning l<5
//./compareTreesMoreRoots 6 test_shallowTrackAndClusterFullInfoNoPU2016DatastripProp.root  test_shallowTrackAndClusterFullInfoNoPU2016MCstripProp.root test_shallowTrackAndClusterFullInfo.rootCT_0.80step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.82step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.84step3.root test_shallowTrackAndClusterFullInfo.rootCT_0.86step3.root TOB yes lowPU2016CrossTalkTOBL1tuning > logsMean/lowPU2016CrossTalkTOBL1tuningTOB

//CT other two param TOB
//./compareTreesMoreRoots 6 test_shallowTrackAndClusterFullInfoNoPU2016_newData.root test_shallowTrackAndClusterFullInfoNoPU2016_newMC.root test_shallowTrackAndClusterFullInfo.rootCT3_0.8401step3.root test_shallowTrackAndClusterFullInfo.rootCT3_0.8402step3.root test_shallowTrackAndClusterFullInfo.rootCT3_0.8403step3.root test_shallowTrackAndClusterFullInfo.rootCT3_0.8404step3.root TOB yes lowPU2016CrossTalkTOBL1tuningFiner > logsMean/lowPU2016CrossTalkTOBL1tuningFinerTOB

//one layer test
 //./compareTreesMoreRoots 2 test_shallowTrackAndClusterFullInfoNoPU2016_newData.root test_shallowTrackAndClusterFullInfoNoPU2016_newMC.root TOB yes lowPU2016layer3dataMC > logsMean/lowPU2016layer3dataMCTOB

//TOB tuned
 //./compareTreesMoreRoots 3 test_shallowTrackAndClusterFullInfoNoPU2016_newData.root test_shallowTrackAndClusterFullInfoNoPU2016_newMC.root test_shallowTrackAndClusterFullInfo.rootCT3_0.805step3.root TOB yes lowPU2016dataMCTOBtuned > logsMean/lowPU2016dataMCTOBtunedTOB


//TOB tuned more stats
 //./compareTreesMoreRoots 3 test_shallowTrackAndClusterFullInfoNoPU2016_newData.root test_shallowTrackAndClusterFullInfo.rootCT_defstep3.root test_shallowTrackAndClusterFullInfo.rootCT3_0.805step3.root TOB yes TOBdataDefMCfromScratchTunedMC

//TOBCRUZET tuned
//
 //./compareTreesMoreRoots 3 test_shallowTrackClusterUP17CosmicsData0T.root test_shallowTrackClusterUP17CosmicsMC0T.root test_shallowTrackAndClusterFullInfo.rootCTcruzet_0.805step3.root TOB yes CRUZETdataMCTOBtuned > logsMean/CRUZETdataMCTOBtunedTOB

//TOBCRAFT tuned
//
 //./compareTreesMoreRoots 3 test_shallowTrackCluster2017CRAFTCosmicsRun298647.root test_shallowTrackAndClusterCRAFT_CRAFTMCdefFromDAS.root test_shallowTrackAndClusterFullInfo.rootCTcruzet_0.805step3.root TOB yes CRAFTTdataMCTOBtuned > logsMean/CRAFTdataMCTOBtunedTOB

//gain options/
 //./compareTreesMoreRoots 5 test_shallowTrackAndClusterNarrowInfoNoPU2016_data.root test_shallowTrackAndClusterNarrowInfo.rootCT_defstep3.root test_shallowTrackAndClusterNarrowInfo.rootgainOption1step3.root test_shallowTrackAndClusterNarrowInfo.rootgainOption2step3.root test_shallowTrackAndClusterNarrowInfo.rootgainOption3step3.root  TOB yes noPUdataGainOptions > logsMean/noPUdataGainOptionsTOB

//TOB tuned more stats - DO NOT USE
 //./compareTreesMoreRoots 6 test_shallowTrackAndClusterFullInfoNoPU2016_newData.root test_shallowTrackAndClusterFullInfo.rootCT_defstep3.root test_shallowTrackAndClusterFullInfo.rootCT3_0.805step3.root test_shallowTrackAndClusterNarrowInfo.rootgainOption1Tuned0.805step3.root test_shallowTrackAndClusterNarrowInfo.rootgainOption2Tuned0.805step3.root  test_shallowTrackAndClusterNarrowInfo.rootgainOption3Tuned0.805step3.root TOB yes TOBdataDefMCfromScratchTunedGainOptionsMC

 //./compareTreesMoreRoots 5 test_shallowTrackAndClusterFullInfoNoPU2016_newData.root test_shallowTrackAndClusterFullInfo.rootCT3_0.805step3.root test_shallowTrackAndClusterNarrowInfo.rootgainOption1Tuned0.805step3.root test_shallowTrackAndClusterNarrowInfo.rootgainOption2Tuned0.805step3.root  test_shallowTrackAndClusterNarrowInfo.rootgainOption3Tuned0.805step3.root TOB yes TOBdataDefMCfromScratchTunedGainOptionsMC > logsMean/noPUdataGainOptionsTOBtuned

//check the reco NOISE
 //./compareTreesMoreRoots 3 test_shallowTrackAndClusterNarrowInfoNoPU2016_data.root test_shallowTrackAndClusterNarrowInfo.rootgainOptionNoiseTest2step3.root test_shallowTrackAndClusterNarrowInfo.rootgainOptionNoiseTest23step3.root  TOB yes noPUdataGainOptionsCheckRECO > logsMean/noPUdataGainOptionsTOBcheckRECO
//

// etoADC
 //./compareTreesMoreRoots 6 test_shallowTrackAndClusterNarrowInfoNoPU2016_data.root test_shallowTrackAndClusterNarrowInfo.rootCT_defstep3.root test_shallowTrackAndClusterNarrowInfo.rootelectronPerAdc_min20_198.0step3.root test_shallowTrackAndClusterNarrowInfo.rootelectronPerAdc_min10_222.0step3.root test_shallowTrackAndClusterNarrowInfo.rootelectronPerAdc_pl10_272.0step3.root test_shallowTrackAndClusterNarrowInfo.rootelectronPerAdc_pl20_296step3.root  TOB yes noPUdataElectronPerADC > logsMean/noPUdataElectronPerADCTOB
 
//CRUZET option3
//./compareTreesMoreRoots 3 test_shallowTrackAndClusterNarrowInfoCRUZETdataMoreStatsRunGT.root test_shallowTrackAndClusterNarrowInfoCRUZET.rootCTcruzet_0.805step3.root test_shallowTrackAndClusterNarrowInfoCRUZET.rootgainOption3CTcruzet_0.805step3.root TOB yes distributionsCRUZETDataToTunedCTTOB05Option3
//
//no PU data/MC with noises
 //./compareTreesMoreRoots 2  test_shallowTrackAndClusterSeedInfoNoPU2016_data.root test_shallowTrackAndClusterSeedInfoNoPU2016_MC.root  TOB yes noPUdataNoises > logsMean/noPUdataNoisesTOB
 //./compareTreesMoreRoots 2  test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root test_shallowTrackAndClusterSeedInfo.rootgainOptionNoiseTestFullstep3.root  TOB yes noPUdataNoises > logsMean/noPUdataNoisesTOB


//Options corrected
//
 //./compareTreesMoreRoots 3 test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root test_shallowTrackAndClusterSeedInfo.rootCT_defstep3.root test_shallowTrackAndClusterSeedInfo.rootgainOption3step3.root  TOB yes noPUdataGainOptions > logsMean/noPUdataGainOptionsTOB
//
 //./compareTreesMoreRoots 3 test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root test_shallowTrackAndClusterSeedInfo.rootCT3_0.805step3.root test_shallowTrackAndClusterSeedInfo.rootgainOption3Tuned0.805step3.root  TOB yes TOBdataDefMCfromScratchTunedGainOptionsMC > logsMean/noPUdataGainOptionsTOBtuned


//e to ADC
 //./compareTreesMoreRoots 6 test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root  test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_min20_198.0step3.root test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_min10_222.0step3.root test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_pl10_272.0step3.root test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_pl20_296.0step3.root test_shallowTrackAndClusterSeedInfo.rootelectronPerAdc_def_247.0step3.root TOB yes noPUdataElectronPerADC > logsMean/noPUdataElectronPerADCTOB

//e to ADC tuned
 //./compareTreesMoreRoots 6 test_shallowTrackAndClusterSeedInfoNoPU2016_dataTest.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC198Tuned0.805step3.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC222Tuned0.805step3.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC272Tuned0.805step3.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC296Tuned0.805step3.root test_shallowTrackAndClusterSeedInfo.rootgainOption3ePerADC247Tuned0.805step3.root TOB yes noPUdataElectronPerADCtuned > logsMean/noPUdataElectronPerADCTOBtuned

 //@MJ@ TODO get rid of layer 5 cut!!!!!!
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

    //if(argc != 7)
    //    throw std::runtime_error("Bad number of arguments!");
    
    TString nfiles = argv[1];
    uint32_t nfile = nfiles.Atoi();
    TString file1 = argv[2];

    vector<TString> fileSim;
    for(uint32_t t = 0; t< nfile-1; t++)
    {
        fileSim.push_back(argv[t+3]);
    }
    TString subDet = argv[nfile+2];
    TString norm = argv[nfile+3];
    TString dir = argv[nfile+4];

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

    vector<TFile*> f2;
    vector<TTree*> t2;
    for(uint32_t t =0; t<fileSim.size(); t++)
    {
        TFile* fle = NULL;
        fle =  TFile::Open(fileSim.at(t));
        f2.push_back(fle);
    }
 
    if(f1==NULL)
        throw std::runtime_error("File 1 address not set");
    t1 = dynamic_cast<TTree*>(f1->Get("testTree/tree"));
    if(t1==NULL)
        throw std::runtime_error("Tree 1 address not set");

    for(uint32_t t =0; t<f2.size(); t++)
    {
       if(f2.at(t)==NULL)
           throw std::runtime_error("File 2 address not set t=" +t );
       TTree * tre =NULL;
       tre  = dynamic_cast<TTree*>(f2.at(t)->Get("testTree/tree"));
       if(tre==NULL)
       {
           TString  message = TString("Tree 2 address not set file ") + to_string(t)+ "  " + fileSim.at(t);
           throw std::runtime_error(message);
       }
       t2.push_back(tre);

    }

    TObjArray *branchList = t1->GetListOfBranches();
    TObjArray *branchListT2 = t2.at(0)->GetListOfBranches();

    if( branchList->GetEntries() != branchListT2->GetEntries() ) //@MJ@ TODO could be more sophiticated, liek comparing elements
        throw std::runtime_error("Trees have different content");

    vector<string> observables;
    for(int i = 0; i < branchList->GetEntries(); ++i) 
    { 
       if((string) branchList->At(i)->GetName() == "nrofevents")
          continue;
       if((string) branchList->At(i)->GetName() == "nroftracks")
          continue;

       observables.push_back(branchList->At(i)->GetName());
            //cout << branchList->At(i)->GetName() << endl; 
    }

    for(uint32_t b=0; b<observables.size(); b++)
    for(std::vector<string>::iterator it = observables.begin() ; it != observables.end(); ++it)
    {
       cout << "observables to be read: " << *it << endl;
       if(*it == "nrofevents" || *it == "nroftracks")
           observables.erase(it);
    }

   if(observables.size() == 0)
       throw std::runtime_error("The vector of observables not filled");

 
   vector< vector<float>* > var;
   vector< vector< vector<float>* > > varT2;
   vector< vector<float> > varTot;
   vector< vector< vector<float> > > varTotT2;
   var.resize(observables.size());
   varTot.resize(observables.size());
   
   varT2.resize(t2.size());
   varTotT2.resize(t2.size());

   uint32_t posSdId = 0;
   uint32_t posSdIdLayer = 0;
   uint32_t puId = 0;


   for(uint32_t v=0;v<var.size();v++)
   {
       if(observables.at(v) == "clustersubdetid")
           posSdId = v; 
       if(observables.at(v) == "clusterlayerwheel")
           posSdIdLayer = v; 
       if(observables.at(v) == "PU")
          puId = v; 
       t1->SetBranchAddress(observables.at(v).c_str(), & (var.at(v)) );
       for(uint32_t t=0;t<t2.size();t++)
       {
           varT2.at(t). resize(observables.size());
           varTotT2.at(t).resize(observables.size());
       
           t2.at(t)->SetBranchAddress(observables.at(v).c_str(), & (varT2.at(t).at(v)) );
       }
   }
   cout << "data hist size " << var.size() << " mc hist size " << varT2.at(0).size() << endl;
   //data always first
    

   uint32_t evCount=0;
   
   //cout << "in here a" << endl;
   Int_t nentries = (Int_t)t1->GetEntries();
   Int_t nentriesT2 = (Int_t)t2.at(0)->GetEntries();

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
                   if(var.at(posSdIdLayer)->at(k)>0 && var.at(posSdIdLayer)->at(k)<5) //@MJ@ TODO temporaray solution!!!!
                       varTot.at(j).push_back(var.at(j)->at(k));
               }
   //cout << "in here d" << endl;
           }
       }

   }

   //fill variables from tree 2
   for (Int_t e=0; e<nentriesT2; e++) 
   {
       for(uint32_t t=0;t<t2.size();t++)
       {
	       t2.at(t)->GetEntry(e);

	       for(uint32_t j=0; j<varT2.at(t).size();j++)
	       {
	   //cout << "in here e" << endl;
		   
		   uint32_t up = varT2.at(t).at(posSdId)->size() < varT2.at(t).at(j)->size() ? varT2.at(t).at(posSdId)->size() :varT2.at(t).at(j)->size();
		   for(uint32_t k=0; k<up;k++)
		   {
		       if(observables.at(j) == "clusternumber")
			   varTotT2.at(t).at(j).push_back(varT2.at(t).at(j)->at(k));

		       if(varT2.at(t).at(posSdId)->at(k) == sdId)
		       {
                           if(varT2.at(t).at(posSdIdLayer)->at(k)>0 && varT2.at(t).at(posSdIdLayer)->at(k)<5) //@MJ@ TODO temporaray solution!!!!
			       varTotT2.at(t).at(j).push_back(varT2.at(t).at(j)->at(k));
		       }
	   //cout << "in here f" << endl;
		   }
	       }

        }

   }

   cout << "in here g" << endl;
   vector<TH1F*> hists;
   vector< vector<TH1F*> > histsT2;
   //cout << "in here h" << endl;
   hists.resize(varTot.size());
   histsT2.resize(t2.size());
   //if(hists.size() != histsT2.at(0).size())
   //    throw std::runtime_error("The histogras do not have same size, something went terribly wrong!");

   for(uint32_t h=0; h<varTot.size(); h++) //@MJ@ TODO in case of empty vector do something!!!
   {
       cout << "in here back" << endl;
       if(varTot.at(h).size()==0) //@MJ@ TODO do it better
       {
           cout << "continuing 1 " << endl;
           continue;
       }
       if(varTotT2.at(0).at(h).size()==0) //@MJ@ TODO do it better
       {
           cout << "continuing 2 " << endl;
           continue;
       }

       sort(varTot.at(h).begin(), varTot.at(h).end());
       for(uint32_t t=0;t<t2.size();t++)
       {
           sort(varTotT2.at(t).at(h).begin(), varTotT2.at(t).at(h).end());
       }

    //cout << "in here 1"  << endl;
       float start = varTot.at(h).at(0) < varTotT2.at(0).at(h).at(0) ? varTot.at(h).at(0) : varTotT2.at(0).at(h).at(0);
    //cout << "in here 1.1"  << endl;
       float size1 = varTot.at(h).size() -1;
    //cout << "in here 1.2"  << endl;
       float size2 = varTotT2.at(0).at(h).size() -1;
    //cout << "in here 1.3"  << endl;
       float end = varTot.at(h).at(size1) > varTotT2.at(0).at(h).at(size2) ? varTot.at(h).at(size1) : varTotT2.at(0).at(h).at(size2);
    //cout << "in here 1.4"  << endl;
       hists.at(h) = new TH1F(observables.at(h).c_str(), observables.at(h).c_str(), 1000, start, end );
       //hists.at(h) = new TH1F(observables.at(h).c_str(), observables.at(h).c_str(), 100, 5, 10 );
    //cout << "in here 1.5"  << endl;
       
       for(uint32_t t=0;t<t2.size();t++)
       {
           histsT2.at(t).resize(varTotT2.at(t).size());
           histsT2.at(t).at(h) = new TH1F((observables.at(h)+"T2").c_str(), (observables.at(h)+"T2").c_str(), 1000, start, end );
           //histsT2.at(t).at(h) = new TH1F((observables.at(h)+"T2").c_str(), (observables.at(h)+"T2").c_str(), 100, 5, 10 );
       }

    cout << "in here 2"  << endl;
       for(uint32_t m = 0; m<varTot.at(h).size(); m++)
       {
           hists.at(h)->Fill(varTot.at(h).at(m));
       }
    cout << "in here 3"  << endl;
       for(uint32_t t=0;t<t2.size();t++)
       {
           for(uint32_t m = 0; m<varTotT2.at(t).at(h).size(); m++)
           {
               histsT2.at(t).at(h)->Fill(varTotT2.at(t).at(h).at(m));
           }
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
       hists.at(h)->SetMaximum(3*maxm);
       hists.at(h)->SetTitle("");
       //hists.at(h)->Draw("*H");
       hists.at(h)->Draw("P");

       double n = hists.at(h)->GetEntries();
       //cout << "norm 1 " << hists.at(h)->GetNormFactor() << " netries " << n << end;
       //double n = 80654/5948; 
       for(uint32_t t=0;t<t2.size();t++)
       {
           //histsT2.at(t).at(h)->SetLineColor(2+(2*t));
           histsT2.at(t).at(h)->SetLineColor(2+(t));
           cout << "hist nr " << t << " color " << 2+(2*t) << endl;
           if(norm == "yes")
           {
               histsT2.at(t).at(h)->DrawNormalized("same hist e", n);
           }
           else
               histsT2.at(t).at(h)->Draw("same hist e");
  
           cout << "hist " << observables.at(h) << " mean 1 " << hists.at(h)->GetMean() << "  mean2 " <<  histsT2.at(t).at(h)->GetMean()<< " mean 1 err " << hists.at(h)->GetMeanError() << "  mean2 err " <<  histsT2.at(t).at(h)->GetMeanError() << endl;
 
       }
   //cout << "in here 6"  << endl;
       pad2.cd();
       vector<TH1F*> histD;
       histD.resize(t2.size());
       vector<double> n2;
       n2.resize(t2.size());
       for(uint32_t t=0;t<t2.size();t++)
       {
           histD.at(t)  = (TH1F*) histsT2.at(t).at(h)->Clone() ; //@MJ@ TODO for all root files
           n2.at(t) = histsT2.at(t).at(h)->GetEntries();
           histD.at(t)->Scale(n/n2.at(t));
           histD.at(t)->Divide(hists.at(h));
           //histD.at(t)->SetLineColor(2+(2*t));
           histD.at(t)->SetLineColor(2+(t));
           histD.at(t)->SetMaximum(2);
           histD.at(t)->SetMinimum(0);
           histD.at(t)->SetTitle("");
           if(t==0)
               histD.at(t)->Draw("");
           else
               histD.at(t)->Draw("same");
       }
       cout << "hist " << h << " drawn" << endl; 

       c.Update();
       c.SaveAs(("output/"+(string)dir+"/"+observables.at(h) + (string)subDet + ".root").c_str());
       //c.SaveAs(("output/"+(string)dir+"/"+observables.at(h) + (string)subDet + ".eps").c_str());
       //delete histD; //do i need that actually?!
   }

    //cout << "in here 7"  << endl;

   return 0;
}
