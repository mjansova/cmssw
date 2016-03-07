// -*- C++ -*-
//
// Package:    SiStripBaselineAnalyzer
// Class:      SiStripBaselineAnalyzer
// 
/**\class SiStripBaselineAnalyzer SiStripBaselineAnalyzer.cc Validation/SiStripAnalyzer/src/SiStripBaselineAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Ivan Amos Cali
//         Created:  Mon Jul 28 14:10:52 CEST 2008
//
//
 

// system include files
#include <memory>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "DataFormats/SiStripDigi/interface/SiStripProcessedRawDigi.h"
#include "DataFormats/SiStripDigi/interface/SiStripRawDigi.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"

#include "CondFormats/SiStripObjects/interface/SiStripPedestals.h"
#include "CondFormats/DataRecord/interface/SiStripPedestalsRcd.h"

#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripPedestalsSubtractor.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripCommonModeNoiseSubtractor.h"
#include "RecoLocalTracker/SiStripZeroSuppression/interface/SiStripRawProcessingFactory.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"

#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"
#include "DataFormats/Provenance/interface/Timestamp.h"

//ROOT inclusion
#include "TROOT.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "THStack.h"


//
// class decleration
// class decleration
//

class SiStripBaselineAnalyzer : public edm::EDAnalyzer {
   public:
      explicit SiStripBaselineAnalyzer(const edm::ParameterSet&);
      ~SiStripBaselineAnalyzer();


   private:
      virtual void beginJob() override ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override ;
      
	  std::auto_ptr<SiStripPedestalsSubtractor>   subtractorPed_;
          edm::ESHandle<SiStripPedestals> pedestalsHandle;
          std::vector<int> pedestals;
          uint32_t peds_cache_id;

          bool plotClusters_;
          bool plotBaseline_;
          bool plotBaselinePoints_;
          bool plotRawDigi_;
          bool plotAPVCM_;
          bool plotPedestals_;
	  
          edm::InputTag srcBaseline_;
          edm::InputTag srcBaselinePoints_;
          edm::InputTag srcAPVCM_;
	  edm::InputTag srcProcessedRawDigi_;
      
          edm::Service<TFileService> fs_;
  
	  TH1F* h1BadAPVperEvent_ = NULL;
	  
          TH1F* h1RawDigis_ = NULL;
	  TH1F* h1ProcessedRawDigis_ = NULL;
	  TH1F* h1Baseline_ = NULL;
	  TH1F* h1BaselineDistr_ = NULL;
	  TH1F* h1Clusters_ = NULL;
          TH1F* h1APVCM_;
          TH1F* h1Pedestals_;	  
          TH1F* h1Baselines_;	  
          TH2D* h2Baselines_;	  

          TH1I* h1ClusterMult_; 
          TH1I* h1ClusterCharge_;
          TH1I* h1ClusterWidth_; 
          TH1I* h1ClusterMean_; 
          TH1I* h1ClusterSigma_; 
	  
	  TCanvas* Canvas_;
	  std::vector<TH1F> vProcessedRawDigiHisto_;
	  std::vector<TH1F> vBaselineHisto_;
          std::vector<TH1F> vBaselinePointsHisto_;
	  std::vector<TH1F> vClusterHisto_;
	  
	  uint16_t nModuletoDisplay_;
	  uint16_t actualModule_;

          std::map<uint32_t, uint32_t> idx;
          std::vector<uint32_t> moduleKeys;
          std::vector<uint32_t> moduleValues;
          std::vector<uint32_t> moduleIdx;
          std::vector< std::vector<float> > baselinesMtx;
          //std::vector<TH1F*> tBaseline_; 
  
          TFile *outFile;
          TTree *timeTree;
          Long64_t tDiff;
          ULong64_t timeT;
          ULong64_t timeH;
          ULong64_t timeL;
          UInt_t runNr;
          UInt_t eventNr;
};


SiStripBaselineAnalyzer::SiStripBaselineAnalyzer(const edm::ParameterSet& conf){
   
  srcBaseline_ =  conf.getParameter<edm::InputTag>( "srcBaseline" );
  srcBaselinePoints_ = conf.getParameter<edm::InputTag>( "srcBaselinePoints" );
  srcProcessedRawDigi_ =  conf.getParameter<edm::InputTag>( "srcProcessedRawDigi" );
  srcAPVCM_ =  conf.getParameter<edm::InputTag>( "srcAPVCM" );
  subtractorPed_ = SiStripRawProcessingFactory::create_SubtractorPed(conf.getParameter<edm::ParameterSet>("Algorithms"));
  nModuletoDisplay_ = conf.getParameter<uint32_t>( "nModuletoDisplay" );
  plotClusters_ = conf.getParameter<bool>( "plotClusters" );
  plotBaseline_ = conf.getParameter<bool>( "plotBaseline" );
  plotBaselinePoints_ = conf.getParameter<bool>( "plotBaselinePoints" );
  plotRawDigi_ = conf.getParameter<bool>( "plotRawDigi" );
  plotAPVCM_ = conf.getParameter<bool>( "plotAPVCM" );
  plotPedestals_ = conf.getParameter<bool>( "plotPedestals" );

  h1BadAPVperEvent_ = fs_->make<TH1F>("BadAPV/Event","BadAPV/Event", 2001, -0.5, 2000.5);
  h1BadAPVperEvent_->SetXTitle("# Modules with Bad APVs");
  h1BadAPVperEvent_->SetYTitle("Entries");
  h1BadAPVperEvent_->SetLineWidth(2);
  h1BadAPVperEvent_->SetLineStyle(2);

  h1APVCM_ = fs_->make<TH1F>("APVCM","APVCM", 2048, -1023.5, 1023.5);
  h1APVCM_->SetXTitle("APV CM [adc]");
  h1APVCM_->SetYTitle("Entries");
  h1APVCM_->SetLineWidth(2);
  h1APVCM_->SetLineStyle(2);

  h1Pedestals_ = fs_->make<TH1F>("Pedestals","Pedestals", 2048, -1023.5, 1023.5);
  h1Pedestals_->SetXTitle("Pedestals [adc]");
  h1Pedestals_->SetYTitle("Entries");
  h1Pedestals_->SetLineWidth(2);
  h1Pedestals_->SetLineStyle(2);
  
  h1Baselines_ = fs_->make<TH1F>("Baselines","Baselines", 2048, -1024, 1024);
  h1Baselines_->SetXTitle("Baselines [adc]");
  h1Baselines_->SetYTitle("Entries");
  h1Baselines_->SetLineWidth(2);
  h1Baselines_->SetLineStyle(2);
  
  h2Baselines_ = fs_->make<TH2D>("BaselinesVsModules","BaselinesVsModules", 1324, -300, 1024, 15000, 0, 15000);
  h2Baselines_->SetXTitle("Baselines [adc]");
  h2Baselines_->SetYTitle("Modules");
  h2Baselines_->SetLineWidth(2);
  h2Baselines_->SetLineStyle(2);
  
  h1ClusterMult_ = fs_->make<TH1I>("ClusterMult","Cluster Multiplicity;nClusters;nEvents", 100, 0, 500000);
  h1ClusterCharge_ = fs_->make<TH1I>("ClusterCharge","Cluster Charge;Cluster Charge;nCluster", 100, 0, 5000);
  h1ClusterWidth_ = fs_->make<TH1I>("ClusterWidth","Cluster Width;Cluster Width;nCluster", 128, 0, 128);
  h1ClusterMean_ = fs_->make<TH1I>("ClusterMean","Cluster Mean;Cluster Mean;nCluster", 128, 0, 128);
  h1ClusterSigma_ = fs_->make<TH1I>("ClusterSigma","Cluster Sigma;Cluster Sigma;nCluster", 60, 0, 50);

  outFile = new TFile("timeInfos.root","recreate");
  timeTree = new TTree("timeTree","timeTree");
  timeTree->Branch("moduleKeys", "std::vector<uint32_t>", &moduleKeys);
  timeTree->Branch("moduleValues", "std::vector<uint32_t>", &moduleValues);
  timeTree->Branch("moduleIdx", "std::vector<uint32_t>", &moduleIdx);
  timeTree->Branch("baselinesMtx", "std::vector< std::vector<float> >", &baselinesMtx);
  timeTree->Branch("tDiff", &tDiff, "tDiff/L");
  timeTree->Branch("timeT", &timeT, "timeT/l");
  timeTree->Branch("timeH", &timeH, "timeH/l");
  timeTree->Branch("timeL", &timeL, "timeL/l");
  timeTree->Branch("runNr", &runNr, "runNr/i");
  timeTree->Branch("eventNr", &eventNr, "eventNr/i");
}


SiStripBaselineAnalyzer::~SiStripBaselineAnalyzer()
{
 
   

}

void
SiStripBaselineAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
   //std::cout << "NEW EVENT!" << std::endl; 
   moduleKeys.clear();
   moduleValues.clear();
   moduleIdx.clear();
   baselinesMtx.clear();
   
   bool ClusterDists = false;
   using namespace edm;
   if(plotPedestals_&&actualModule_ ==0){
      uint32_t p_cache_id = es.get<SiStripPedestalsRcd>().cacheIdentifier();
      if(p_cache_id != peds_cache_id) {
	es.get<SiStripPedestalsRcd>().get(pedestalsHandle);
	peds_cache_id = p_cache_id;
      }
      
      
      std::vector<uint32_t> detIdV;
      pedestalsHandle->getDetIds(detIdV);
      for(uint32_t i=0; i < detIdV.size(); ++i){
        pedestals.clear();
        SiStripPedestals::Range pedestalsRange = pedestalsHandle->getRange(detIdV[i]);
        pedestals.resize((pedestalsRange.second- pedestalsRange.first));
	pedestalsHandle->allPeds(pedestals, pedestalsRange);
	for(uint32_t it=0; it < pedestals.size(); ++it) h1Pedestals_->Fill(pedestals[it]);
      }
   }

   if(plotAPVCM_){
     edm::Handle<edm::DetSetVector<SiStripProcessedRawDigi> > moduleCM;
     edm::InputTag CMLabel("siStripZeroSuppression:APVCM");
     e.getByLabel(srcAPVCM_,moduleCM);

     edm::DetSetVector<SiStripProcessedRawDigi>::const_iterator itCMDetSetV =moduleCM->begin();
     for (; itCMDetSetV != moduleCM->end(); ++itCMDetSetV){  
       edm::DetSet<SiStripProcessedRawDigi>::const_iterator  itCM= itCMDetSetV->begin();
       for(;itCM != itCMDetSetV->end(); ++itCM) h1APVCM_->Fill(itCM->adc()); //@MJ@ TODO filling of baselines
     }
   }
   if(!plotRawDigi_) return;
   subtractorPed_->init(es); 
 
   edm::Handle< edm::DetSetVector<SiStripRawDigi> > moduleRawDigi;
   e.getByLabel(srcProcessedRawDigi_,moduleRawDigi);
 
   edm::Handle<edm::DetSetVector<SiStripProcessedRawDigi> > moduleBaseline;
   if(plotBaseline_) e.getByLabel(srcBaseline_, moduleBaseline); 

//here (2 lines)
   edm::Handle<edm::DetSetVector<SiStripProcessedRawDigi> > moduleBaselinePoints;
   if(plotBaselinePoints_) e.getByLabel(srcBaselinePoints_, moduleBaselinePoints); 
   
   edm::Handle<edmNew::DetSetVector<SiStripCluster> > clusters;
   if(plotClusters_){
   	edm::InputTag clusLabel("siStripClusters");
   	e.getByLabel(clusLabel, clusters);
   }

   char detIds[20];
   char evs[20];
   char runs[20];    
   char times[20];    
   

   TFileDirectory sdProcessedRawDigis_= fs_->mkdir("ProcessedRawDigis");
   TFileDirectory sdRawDigis_= fs_->mkdir("RawDigis");
   TFileDirectory sdBaseline_= fs_->mkdir("Baseline");
   TFileDirectory sdBaselineDistr_= fs_->mkdir("BaselineDistr");
   TFileDirectory sdtBaseline_= fs_->mkdir("tBaseline");
   TFileDirectory sdBaselinePoints_= fs_->mkdir("BaselinePoints");
   TFileDirectory sdClusters_= fs_->mkdir("Clusters");
   

//here
//


   edm::RunNumber_t const run = e.id().run();
   edm::EventNumber_t const event = e.id().event();
   runNr = static_cast<UInt_t>(run);
   eventNr = static_cast<UInt_t>(event);
   
   edm::Timestamp const timestamp = e.time();
   edm::TimeValue_t timeVal = timestamp.value();
   timeT = static_cast<ULong64_t>(timeVal);
   static const ULong64_t t = timeT;
   tDiff = t - timeT;
   timeH = timestamp.unixTime();
   timeL = timestamp.microsecondOffset();

   //ULong64_t time2= timeH;
   //time2 = time2 << 32;
   //time2 += timeL;
   //std::cout << "tL: " << timeL << " t " << time <<  " t2 " << time2 << std::endl;

   edm::DetSetVector<SiStripProcessedRawDigi>::const_iterator itDSBaseline;
   if(plotBaseline_) itDSBaseline = moduleBaseline->begin();
   edm::DetSetVector<SiStripRawDigi>::const_iterator itRawDigis = moduleRawDigi->begin();
   
   uint32_t NBabAPVs = moduleRawDigi->size();     
   h1BadAPVperEvent_->Fill(NBabAPVs);
   
   for (; itRawDigis != moduleRawDigi->end(); ++itRawDigis) {
      //std::cout << "modules to display: " << nModuletoDisplay_ << std::endl;
      if(actualModule_ > nModuletoDisplay_) {return;}
      uint32_t detId = itRawDigis->id;
      moduleIdx.push_back(detId);
      
      //std::cout << "time : " << time << "t diff" << tDiff << std::endl;
      sprintf(detIds,"%u", detId);
      sprintf(evs,"%llu", event);
      sprintf(runs,"%u", run);
      sprintf(times,"%llu", timeT);
      
      uint32_t mapSize = idx.size();
      //std::cout << "map size is: " << mapSize << std::endl;
      /*std::pair<std::map<uint32_t,uint32_t>::iterator,bool>  newElem =*/
      idx.insert(std::pair<uint32_t, uint32_t>(detId,mapSize+1));
      /* if(newElem.second)
      {
        TH1F* h1tmp = sdtBaseline_.make<TH1F>(detIds, detIds, 30000, -30000, 30000);
        h1tmp->SetXTitle("time");
	h1tmp->SetYTitle("baseline");
	h1tmp->SetMaximum(20);
	h1tmp->SetMinimum(-30);
	h1tmp->SetLineWidth(2);
	h1tmp->SetLineStyle(2);
	h1tmp->SetLineColor(kPink+7);
        tBaseline_.push_back(h1tmp);
      }
      */
      //std::cout << "detId: " << detId << std::endl;
	  
      if(plotBaseline_){
//	std::cout << "bas id: " << itDSBaseline->id << " raw id: " << detId << std::endl;
	if(itDSBaseline->id != detId){
		itDSBaseline = moduleBaseline->find(detId);
                if(itDSBaseline->id != detId){ if(plotBaseline_)itDSBaseline++; continue; }
//                else std::cout << "Resynched..." << std::endl;
	}	  
      }
      
    
      actualModule_++;
      //std::cout << "event nr: " << e.id().event() << "run nr: " <<  e.id().run()<< std::endl;
      //std::cout << "processing module N: " << actualModule_<< " detId: " << detId << " event: "<< event << std::endl; 

	  
      edm::DetSet<SiStripRawDigi>::const_iterator itRaw = itRawDigis->begin(); 
      bool restAPV[6] = {0,0,0,0,0,0};
      int strip =0, totADC=0;
      int minAPVRes = 7, maxAPVRes = -1;
      for(;itRaw != itRawDigis->end(); ++itRaw, ++strip){
	    float adc = itRaw->adc();
	    totADC+= adc;
	    if(strip%128 ==127){
      		//std::cout << "totADC " << totADC << std::endl;
	      int APV = strip/128;
	      if(totADC!= 0){
      	    	restAPV[APV] = true;
      			totADC =0;
      			if(APV>maxAPVRes) maxAPVRes = APV;
      			if(APV<minAPVRes) minAPVRes = APV;
      	      }
	    }
      }

      uint16_t bins =768;
      float minx = -0.5, maxx=767.5;
      if(minAPVRes !=7){
      	minx = minAPVRes * 128 -0.5;
      	maxx = maxAPVRes * 128 + 127.5;
      	bins = maxx-minx;
      }
      

      char* dHistoName = Form("Id%s_run%s_ev%s_t%s",detIds, runs, evs, times);
      //char* dHistoName2 = Form("Id%s_run%s"s);
      h1ProcessedRawDigis_ = sdProcessedRawDigis_.make<TH1F>(dHistoName,dHistoName, bins, minx, maxx); 
      h1RawDigis_ = sdRawDigis_.make<TH1F>(dHistoName,dHistoName, bins, minx, maxx); 

      edm::DetSet<SiStripRawDigi>::const_iterator itRaw2 = itRawDigis->begin(); 
      int strip2=0;
      for(; itRaw2 != itRawDigis->end(); ++itRaw2, ++strip2){
            h1RawDigis_->Fill(strip2,itRaw2->adc());
      }
      
      if(plotBaseline_){
	h1Baseline_ = sdBaseline_.make<TH1F>(dHistoName,dHistoName, bins, minx, maxx); 
        h1Baseline_->SetXTitle("strip#");
	h1Baseline_->SetYTitle("ADC");
	h1Baseline_->SetMaximum(1024.);
	h1Baseline_->SetMinimum(-300.);
	h1Baseline_->SetLineWidth(2);
	h1Baseline_->SetLineStyle(2);
	h1Baseline_->SetLineColor(2);
      }
      
      if(plotBaseline_){
	h1BaselineDistr_ = sdBaselineDistr_.make<TH1F>(dHistoName,dHistoName, 1324, -300., 1024); 
        h1BaselineDistr_->SetXTitle("baseline");
	h1BaselineDistr_->SetYTitle("count");
	h1BaselineDistr_->SetMaximum(100);
	h1BaselineDistr_->SetMinimum(0);
	h1BaselineDistr_->SetLineWidth(2);
	h1BaselineDistr_->SetLineStyle(2);
	h1BaselineDistr_->SetLineColor(9);
      }

      if(plotClusters_){
        h1Clusters_ = sdClusters_.make<TH1F>(dHistoName,dHistoName, bins, minx, maxx);
	  
        h1Clusters_->SetXTitle("strip#");
        h1Clusters_->SetYTitle("ADC");
        h1Clusters_->SetMaximum(1024.);
        h1Clusters_->SetMinimum(-300.);
        h1Clusters_->SetLineWidth(2);
	h1Clusters_->SetLineStyle(2);
	h1Clusters_->SetLineColor(3);
      }

      h1ProcessedRawDigis_->SetXTitle("strip#");  
      h1ProcessedRawDigis_->SetYTitle("ADC");
      h1ProcessedRawDigis_->SetMaximum(1024.);
      h1ProcessedRawDigis_->SetMinimum(-300.);
      h1ProcessedRawDigis_->SetLineWidth(2);
      h1RawDigis_->SetXTitle("strip#");  
      h1RawDigis_->SetYTitle("ADC");
      h1RawDigis_->SetMaximum(1024.);
      h1RawDigis_->SetMinimum(-300.);
      h1RawDigis_->SetLineWidth(2);

       
      std::vector<int16_t> ProcessedRawDigis(itRawDigis->size());
      subtractorPed_->subtract( *itRawDigis, ProcessedRawDigis);
      //std::cout << "raw digi size: " << itRawDigis->size() << std::endl;

      edm::DetSet<SiStripProcessedRawDigi>::const_iterator  itBaseline;
      if(plotBaseline_) itBaseline = itDSBaseline->begin(); 
      
      std::vector<int16_t>::const_iterator itProcessedRawDigis;
     
      strip =0;
      int32_t num = 0;
      int32_t denom = 0;  
      std::vector<float> blsInModule;  
      for(itProcessedRawDigis = ProcessedRawDigis.begin();itProcessedRawDigis != ProcessedRawDigis.end(); itProcessedRawDigis++){ 
       	if(restAPV[strip/128]){
          //std::cout << "bla" << std::endl;
	  float adc = *itProcessedRawDigis;     
	  h1ProcessedRawDigis_->Fill(strip, adc);
	  if(plotBaseline_){
            //std::cout << "filling baseline " << std::endl;
	    h1Baseline_->Fill(strip, itBaseline->adc());
            //std::cout << "strip: " << strip << " moduel: " << idx.find(detId)->second << " baselien: " << itBaseline->adc() << std::endl;
            if(strip%128 == 0)
            {
                //std::cout << "strip: " << strip << std::endl;
                //std::cout << "filling baseline 2" << std::endl;
                //@MJ@ TODO fill histogram
                //std::cout << "filling hist h1baseline distr, bin:  " << 300 +  itBaseline->adc() << std::endl;
                //std::cout << "element in map: " << idx.find(detId)->second << std::endl;
	        //h1BaselineDistr_->Fill(itBaseline->adc());
	        h1Baselines_->Fill(itBaseline->adc());
	        h2Baselines_->Fill(itBaseline->adc(), idx.find(detId)->second);
                num = num + itBaseline->adc();
                denom++;
                blsInModule.push_back(itBaseline->adc());
            }
            //std::cout << "det Id: " << detId << "strip " << strip << " baseline: " << itBaseline->adc() << std::endl;
	    ++itBaseline;
	  }
	 }
	++strip;
      }
      //tBaseline_.at((idx.find(detId)->second)-1)->Fill(tDiff, num/denom);
      baselinesMtx.push_back(blsInModule);	  
      if(plotBaseline_) ++itDSBaseline; 
      if(plotClusters_){
          int nclust = 0;
	  edmNew::DetSetVector<SiStripCluster>::const_iterator itClusters = clusters->begin();
	  for ( ; itClusters != clusters->end(); ++itClusters ){
		for ( edmNew::DetSet<SiStripCluster>::const_iterator clus =	itClusters->begin(); clus != itClusters->end(); ++clus){
		  if(itClusters->id() == detId){
		    int firststrip = clus->firstStrip();
	            //std::cout << "Found cluster in detId " << detId << " " << firststrip << " " << clus->amplitudes().size() << " -----------------------------------------------" << std::endl;		
     		    strip=0;
		    for( auto itAmpl = clus->amplitudes().begin(); itAmpl != clus->amplitudes().end(); ++itAmpl){
		      h1Clusters_->Fill(firststrip+strip, *itAmpl);
		      ++strip;
		    }
		  }
                  
                  //cluster plots from here on
                   
                  if(ClusterDists == false){
		    nclust++;

     		    int strip2=0;
                    double charge = 0;
                    double mean = 0;
                    double sigma = 0;
		    for( auto itAmpl = clus->amplitudes().begin(); itAmpl != clus->amplitudes().end(); ++itAmpl){
                      charge += *itAmpl;
		      ++strip2;
                      mean += strip2*(*itAmpl);
		      sigma += strip2*strip2*(*itAmpl);
		    }
                    h1ClusterCharge_->Fill(charge);
                    h1ClusterWidth_->Fill(strip2);
                    mean = mean/charge;
                    h1ClusterMean_->Fill(mean);
                    sigma = TMath::Power((sigma/charge-mean*mean),0.5);
                    h1ClusterSigma_->Fill(sigma);
                  }              
		}
	  }
        if(ClusterDists==false)
        {
          h1ClusterMult_->Fill(nclust); 
          ClusterDists = true;
        }
      }

    if (plotBaseline_)
    {
        delete h1Baseline_;
        delete  h1BaselineDistr_;
        h1Baseline_ = NULL;
        h1BaselineDistr_ = NULL;
    }
    if (plotClusters_)
    {
        delete h1Clusters_;
        h1Clusters_ = NULL;
    }
    delete h1RawDigis_;
    delete h1ProcessedRawDigis_;
    h1RawDigis_ = NULL;
    h1ProcessedRawDigis_ = NULL;
 
    }		

actualModule_ =0;  
//tDiff =3;
//std::cout << "t DIff: " << tDiff << std::endl;

//std::cout << "tL: " << timeL << " t " << timeT << std::endl;
std::transform( idx.begin(), idx.end(), std::back_inserter( moduleKeys ), [](std::pair<uint32_t, uint32_t> const & p) { return p.first; } );
std::transform( idx.begin(), idx.end(), std::back_inserter( moduleValues ), [](std::pair<uint32_t, uint32_t> const & r) { return r.second; } );
timeTree->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void SiStripBaselineAnalyzer::beginJob()
{
  
  
actualModule_ =0;  
   
 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
SiStripBaselineAnalyzer::endJob() {

    timeTree->Print();
    outFile->Write();     
}

//define this as a plug-in
DEFINE_FWK_MODULE(SiStripBaselineAnalyzer);

