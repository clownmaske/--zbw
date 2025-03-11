Example1
void ReadTrack()
{
	// Add File
	TChain* to = new TChain("SimTriggerInfo","");
	to->Add("Your output file of the Simulation");
	
	vector<JPSimTruthTree_t>* truthList = new vector<JPSimTruthTree_t>;
	to->SetBranchAddress("truthList", &truthList);
	for(Int_t i = 0 ; i < to->GetEntries(); i++)
	{
		to->GetEntry(i);
		for(Int_t j = 0 ; j < truthList->size(); j ++)
		{	
			// Get Track
			std::vector<JPSimTrack_t> trackList = (*truthList)[j].trackList;
			for(Int_t k = 0 ; k < trackList.size() ; k++)
			{
				// Get Step Point
				std::vector<JPSimStepPoint_t> StepPoints = trackList.at(k).StepPoints;
				// Get the number of the detected scintillation PE
				if(trackList.at(k).nPdgId == 0&&StepPoints.at(0).nProcessSubType == 22&&trackList.at(k).bDetectedPhoton == 1)	
          			cout<<"StepPoint number is "<<StepPoints.size()<<endl;
				// Loop the StepPoints and get the position of each step point
				for(Int_t m = 0 ; m < StepPoints.size(); m++)
          			{
                			x = StepPoints.at(m).fX;                    
                			y = StepPoints.at(m).fY;               
                			z = StepPoints.at(m).fZ; 
				}
			}
		}
	}
}


Example2
#include <vector>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TROOT.h"
//#include "JPSimOutput.hh"
using namespace std;

void ReadWaveform(TString inputfile="JPSim_output.root", TString outputfile="JPSim_waveform.root")
{  
	// Load the library
	//gSystem->Load("../../utils/libJPSIMOUTPUT.so");

	// Open the output file of JPSim
	TFile* f = new TFile(inputfile);

	//Create a file to save the waveform
	TFile* f1 = new TFile(outputfile,"recreate");

	// Get the SimTrigger Tree
	TTree* t = (TTree*)f->Get("Readout");

	std::vector<UInt_t>* ChannelId = new std::vector<UInt_t>;
	std::vector<UInt_t>* Waveform = new std::vector<UInt_t>;

	// Create the two branches to be read

	t->SetBranchAddress("ChannelId", &ChannelId);
	t->SetBranchAddress("Waveform", &Waveform);

	// Get windowSize
	t->GetEntry(0);
	Int_t windowSize = Waveform->size()/ChannelId->size();

	// Read the TriggerReadout tree
	for(Int_t i=0; i<t->GetEntries(); i++)
	{
		t->GetEntry(i);

		// Create the superposition histogram
		char name[100];
		char name1[100];
		sprintf(name,"Trigger_%d", i);
		sprintf(name1,"Superposition_%d", i);
		//cout<<name1<<endl;
		TH1F* h0 = new TH1F(name1, name1, windowSize, 0, windowSize);
		f1->mkdir(name);
		f1->cd(name);

		for(size_t j=0; j<ChannelId->size(); j++)
		{
			UInt_t PMTId = ChannelId->at(j);

			char pmtname[1000];
			sprintf(pmtname,"PMT_%d",PMTId);
			TH1F* h1 = new TH1F(pmtname,pmtname, windowSize, 0, windowSize);
			
			for(Int_t k=1; k<=windowSize; k++)
			{
				h1->SetBinContent(k, Waveform->at(k-1+windowSize*j));
			}

			h0->Add(h1);
		}

		f1->cd("../");
	
	}


	f1->Write();
	f1->Close();
	f->Close();

}