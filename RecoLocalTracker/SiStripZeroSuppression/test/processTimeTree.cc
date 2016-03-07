#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <sstream> 

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

int main(int argc, char *argv[])
{
    TFile *f = new TFile("timeInfos.root");
    TTree *t1 = (TTree*)f->Get("timeTree");

    vector<unsigned int> *moduleKeys = NULL;
    vector<unsigned int> *moduleValues = NULL;
    vector<unsigned int> *moduleIdx = NULL;
    vector< vector<float> > *baselinesMtx = NULL;
    Long64_t tDiff;
    ULong64_t timeT; 
    ULong64_t timeH; 
    ULong64_t timeL;
    uint32_t runNr; 
    uint32_t eventNr; 

    t1->SetBranchAddress("tDiff",&tDiff);
    t1->SetBranchAddress("moduleKeys",&moduleKeys);
    t1->SetBranchAddress("moduleValues",&moduleValues);
    t1->SetBranchAddress("baselinesMtx",&baselinesMtx);
    t1->SetBranchAddress("moduleIdx",&moduleIdx);
    t1->SetBranchAddress("timeT",&timeT);
    t1->SetBranchAddress("timeH",&timeH);
    t1->SetBranchAddress("timeL",&timeL);
    t1->SetBranchAddress("runNr",&runNr);
    t1->SetBranchAddress("eventNr",&eventNr);
  
    vector<ULong64_t> timeSorted;
    vector<Long64_t> timeSortedDiff;
    std::map<ULong64_t,uint32_t> timeMap;
    std::map<ULong64_t,uint32_t> timeMapDiff;
    std::map< uint32_t, uint32_t > detIdMap;
    std::map<uint32_t, TGraph*> graphMap;

    Int_t nentries = (Int_t)t1->GetEntries();
   
    TH1D* timeDiff = NULL;
    timeDiff = new TH1D("timeDiff", "timeDiff", 8590150, 0, 17180300289); //@MJ@ TODO this is nasty!
    //read all entries and fill the histograms
    for (Int_t i=0; i<nentries; i++) 
    {
         t1->GetEntry(i);
         timeSorted.push_back(timeT);
         timeMap.insert( std::pair<ULong64_t,uint32_t>(timeT,i) );
         
         std::transform( moduleValues->begin(), moduleValues->end(), moduleKeys->begin(), 
                    std::inserter(detIdMap, detIdMap.end() ), std::make_pair<uint32_t const&,uint32_t const&> );
   /*     for(uint32_t det = 10; det < 50; det++ ) //@MJ@ TODO user defined
        {
            if(i == 0)
            {
                //TGraph* graph = new TGraph(nentries*(baselinesMtx->at(det - 1).size()));
                TGraph* graph = new TGraph(50);//(nentries);
                graph->SetName((TString) static_cast<ostringstream*>( &(ostringstream() << det) )->str()); //@MJ@ TODO real det ID!
                graphMap.insert( std::pair<uint32_t, TGraph*>(det,graph) );
            
            }
            std::vector<float> baselines =  baselinesMtx->at(det - 1);
            //uint32_t point = i*baselines.size(); //@MJ TODO@
            for(uint32_t b = 0; b<baselines.size(); b++)
            {
                if(b==3 && i < 50)
                (graphMap.find(det)->second)->SetPoint(i, timeT, baselines.at(b));//+b*10);
            }
        }
*/
        //std::cout << "event: " << eventNr << " time: " << timeT << std::endl;
    }

    ULong64_t prev = 0;
    uint32_t j = 0;
    for (std::map<ULong64_t,uint32_t>::iterator it=timeMap.begin(); it!=timeMap.end(); ++it)
    {
        std::cout << it->first << " => " << it->second << '\n';
        t1->GetEntry(it->second);
        
        for(uint32_t det = 10; det < 50; det++ ) //@MJ@ TODO user defined
        {
            if(prev == 0)
            {
                //Graph* graph = new TGraph(nentries*(baselinesMtx->at(det - 1).size()));
                TGraph* graph = new TGraph(nentries);
                graph->SetName((TString) static_cast<ostringstream*>( &(ostringstream() << det) )->str()); //@MJ@ TODO real det ID!
                graphMap.insert( std::pair<uint32_t, TGraph*>(det,graph) );
            
            }
            std::vector<float> baselines =  baselinesMtx->at(det - 1);
            //uint32_t point = i*baselines.size(); //@MJ TODO@
            for(uint32_t b = 0; b<baselines.size(); b++)
            {
                if(b==2)
                {
                (graphMap.find(det)->second)->SetPoint(j, it->first , baselines.at(b));//+b*10);
                //cout << "baseline: " << baselines.at(b) << "time: " << it->first  << endl;
                }
            }
        }

        j++;


        if(prev != 0)
            timeMapDiff.insert( std::pair<ULong64_t,uint32_t> (it->first - prev,it->second));
        prev = it->first;
    }
    for (std::map<ULong64_t,uint32_t>::iterator it2=timeMapDiff.begin(); it2!=timeMapDiff.end(); ++it2)
    {
        std::cout << it2->first << " => " << it2->second << '\n';
        Double_t binContent = static_cast<Double_t>(it2->first);
        cout << "bin content: " << binContent << std::endl;
        timeDiff->Fill(binContent);
    }

   
    //cout << "size unsorted time " << timeSorted.size() << endl;
    std::sort (timeSorted.begin(), timeSorted.end());
    //cout << "size sorted time " << timeSorted.size() << endl;
    
    
    for(uint32_t t = 1; t < timeSorted.size(); t++)
    {
        Long64_t diff  = timeSorted.at(t) - timeSorted.at(t-1);
        timeSortedDiff.push_back(diff);
        //cout << "diff " << diff << endl;
    }
  
    std::sort (timeSortedDiff.begin(), timeSortedDiff.end());
    //cout << "size sorted time diff" << timeSortedDiff.size() << endl;
    
    for (std::vector<Long64_t>::iterator it=timeSortedDiff.begin(); it!=timeSortedDiff.end(); ++it)
        std::cout << "time difference between events: " << *it << std::endl;


    TString rootOut = "output_plots.root";
    TFile fi(rootOut,"RECREATE");
    timeDiff->Write();
    for (std::map<uint32_t, TGraph*>::iterator it3=graphMap.begin(); it3!=graphMap.end(); ++it3)
    {
        it3->second->Write();
    }
}
