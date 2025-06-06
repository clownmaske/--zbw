#include <vector>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TROOT.h"
//#include "TCanvas.h"
#include "JPSimOutput.hh" 
//#include "JPWaveformPreprocess.h"

void IntegralWave(TString inputfile, TString baselinefile, TString outputfile)
{
    //Open inputfile
    TFile* f = new TFile(inputfile);
    TFile* fb = new TFile(baselinefile);

    //Create a file to save the waveform
    TFile* f1 = new TFile(outputfile, "recreate");

    //Get the SimTriggerInfo and Readout Tree
    TTree* t = (TTree*)f->Get("Readout");
    //TTree* t1 = (TTree*)f->Get("SimTriggerInfo");
    TTree* tb = (TTree*)fb->Get("BaselineTree");

    //Claim vector and struct
    std::vector<UInt_t>* ChannelId = new std::vector<UInt_t>;
    std::vector<UInt_t>* Waveform = new std::vector<UInt_t>;
    //vector<JPSimTruthTree_t>* truthList = new vector<JPSimTruthTree_t>;
    //vector<JPSimPE_t>* PEList = new vector<JPSimPE_t>;
    std::vector<Double_t>* Baseline = new std::vector<Double_t>;

    //SetBranchAddress
    t->SetBranchAddress("ChannelId", &ChannelId);
    t->SetBranchAddress("Waveform", &Waveform);
    //t1->SetBranchAddress("PEList", &PEList);
    //t1->SetBranchAddress("truthList", &truthList);
    tb->SetBranchAddress("Baseline", &Baseline);

    //Get windowSize
    t->GetEntry(1);
    //t1->GetEntry(0);
    Int_t windowSize = Waveform->size() / ChannelId->size();
    cout << "WindowSize is " << windowSize << endl;

    TH1F* h1 = new TH1F("IntegralWave", "IntegralWave;Integral(mV*ns);Count", 100, 0, 2000);
    
    //For All Trigger
    for (Int_t i = 0; i < 1920; i++)
    {
        t->GetEntry(i);
        tb->GetEntry(i);

        for (size_t j = 0; j < ChannelId->size(); j++) 
        {
            UInt_t PMTId = ChannelId->at(j);
            Double_t baseline = Baseline->at(j);
            
            //Integral of Waveform
            double Integral_wave = 0.0;
            const Int_t offset = windowSize * j;
            for (Int_t k = 0; k < windowSize; ++k)
            {
                Integral_wave += (baseline - Waveform->at(offset + k));
            }

            // Fill the Histogram
            h1->Fill(Integral_wave);
            cout << "Trigger" << i << " ChannelId" << PMTId << " Integral is " << Integral_wave << std::endl;
        }

    }
    h1->Write();
    f1->Write();
    f1->Close();
    f->Close();
    fb->Close();
}