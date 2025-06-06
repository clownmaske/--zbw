#include <vector>
#include <map>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h" 
#include "TSystem.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "JPSimOutput.hh" 
//#include "JPWaveformPreprocess.h"

void setJPStyle() {
    Int_t jpFont = 42;
    Double_t jpWidth = 2.00;
    Double_t jpTSize = 0.06;

    gROOT->SetStyle("Plain");
    TStyle* jpStyle = new TStyle("jpStyle", "Jinping plots style");

    jpStyle->SetFrameFillColor(0);
    jpStyle->SetFrameBorderMode(0);
    jpStyle->SetPadBorderMode(0);
    jpStyle->SetPadColor(0);
    jpStyle->SetCanvasBorderMode(0);
    jpStyle->SetCanvasColor(0);
    jpStyle->SetStatColor(0);
    jpStyle->SetLegendBorderSize(0);

    jpStyle->SetPadTopMargin(0.135);
    jpStyle->SetPadRightMargin(0.08);
    jpStyle->SetPadBottomMargin(0.16);
    jpStyle->SetPadLeftMargin(0.15);

    jpStyle->SetTextFont(jpFont);
    jpStyle->SetTextSize(jpTSize);
    jpStyle->SetLabelFont(jpFont, "x");
    jpStyle->SetLabelFont(jpFont, "y");
    jpStyle->SetLabelFont(jpFont, "z");
    jpStyle->SetLabelSize(jpTSize, "x");
    jpStyle->SetLabelSize(jpTSize, "y");
    jpStyle->SetLabelSize(jpTSize, "z");
    jpStyle->SetTitleFont(jpFont);
    jpStyle->SetTitleFont(jpFont, "x");
    jpStyle->SetTitleFont(jpFont, "y");
    jpStyle->SetTitleFont(jpFont, "z");
    jpStyle->SetTitleSize(1.2 * jpTSize, "x");
    jpStyle->SetTitleSize(1.2 * jpTSize, "y");
    jpStyle->SetTitleSize(1.2 * jpTSize, "z");
    jpStyle->SetTitleSize(1.5 * jpTSize, "");

    jpStyle->SetNdivisions(505, "x");
    jpStyle->SetNdivisions(505, "y");
    jpStyle->SetNdivisions(505, "z");

    jpStyle->SetLabelOffset(0.010, "X");
    jpStyle->SetLabelOffset(0.010, "Y");
    jpStyle->SetLabelOffset(0.010, "Z");

    TGaxis::SetMaxDigits(3);

    jpStyle->SetOptStat(0);
    jpStyle->SetOptTitle(1);
    jpStyle->SetOptFit(1);

    jpStyle->SetTitleOffset(0.95, "X");
    jpStyle->SetTitleOffset(0.95, "Y");
    jpStyle->SetTitleOffset(0.85, "Z");

    jpStyle->SetTitleFillColor(0);
    jpStyle->SetTitleStyle(0);
    jpStyle->SetTitleBorderSize(0);
    jpStyle->SetTitleFont(jpFont, "title");
    jpStyle->SetTitleX(0.1);
    jpStyle->SetTitleY(0.985);

    gROOT->SetStyle("jpStyle");
    gROOT->ForceStyle();
}

void ReadWaveform2(TString outputfile)
{  

    setJPStyle();

    // ====================== 多文件处理设置 ======================
    const TString input_dir = "/home/zbw/Desktop/Data/";
    const int first_file_num = 0;
    const int last_file_num = 25;

    // 创建TChain处理多文件
    TChain* readoutChain = new TChain("Readout");
    TChain* simTriggerChain = new TChain("SimTriggerInfo");

    // 添加所有输入文件
    for (int i = first_file_num; i <= last_file_num; ++i) {
        TString filename;
        if (i == 0) {
            filename = input_dir + "SlowLS.root";
        }
        else {
            filename = input_dir + Form("SlowLS_%d.root", i);
        }

        readoutChain->Add(filename);
        simTriggerChain->Add(filename);
        std::cout << "Added file: " << filename << std::endl;
    }

    // ====================== 设置分支地址 ======================
    std::vector<UInt_t>* ChannelId = nullptr;
    std::vector<UInt_t>* Waveform = nullptr;
    vector<JPSimPE_t>* PEList = nullptr;
    

    readoutChain->SetBranchAddress("ChannelId", &ChannelId);
    readoutChain->SetBranchAddress("Waveform", &Waveform);
    simTriggerChain->SetBranchAddress("PEList", &PEList);

    // ====================== 创建输出文件 ======================
    TFile* f1 = new TFile(outputfile, "recreate");

    //Get windowSize
    readoutChain->GetEntry(1);
    //t1->GetEntry(0);
    Int_t windowSize = Waveform->size() / ChannelId->size();
    cout << "WindowSize is " << windowSize << endl;


    //Get Double-charge PMTId
    std::vector<std::vector<int>> DoublechargePMTs;
    //for (Int_t i = 0; i < t->GetEntries(); i++)
    const Long64_t totalEntries = simTriggerChain->GetEntries();
    for (Long64_t i = 0; i < totalEntries; ++i)
    {
        simTriggerChain->GetEntry(i);

        //Count PMTs
        std::unordered_map<int, int> pmtCounter;
        for (const auto& pe : *PEList)
        {
            ++pmtCounter[pe.PMTId];
        }

        // Select double-charge PMTs
        std::vector<int> currentUniquePMTs;
        //for (const auto& entry : pmtCounter)
        for (const auto& [pmtId, count] : pmtCounter)
        {
            if (count == 2)
            {
                currentUniquePMTs.push_back(pmtId);
            }
        }
        DoublechargePMTs.push_back(std::move(currentUniquePMTs));
    }
    std::cout << "Double-charge PMTs Ready!" << endl;

	// Read the TriggerReadout tree
    const Long64_t nEntries = readoutChain->GetEntries();
    for (Long64_t i = 0; i < 10; ++i)
    {
        if (i == 3497)continue;
        if (i == 5797)continue;
        if (i == 8246)continue;
        if (i == 9910)continue;

        readoutChain->GetEntry(i);
        simTriggerChain->GetEntry(i);

        std::unordered_map<int, int> pmtCounterCurrent;
        for (const auto& pe : *PEList)
        {
            ++pmtCounterCurrent[pe.PMTId];
        }


        const auto& validPMTs1 = DoublechargePMTs[i];
        std::unordered_map<UInt_t, std::pair<Double_t, Double_t>> pmtTruthTimes1;
        std::unordered_map<UInt_t, int> pmtCounter1;

        // Get PE Number
        for (const auto& pe : *PEList)
        {
            ++pmtCounter1[pe.PMTId];
        }

        // Get Double PEs of Certain PMT
        for (const auto& pe : *PEList)
        {
            const UInt_t pmtId = pe.PMTId;

            // Only Get Double-charge PMTs
            if (pmtCounter1[pmtId] != 2) continue;
            if (!std::count(validPMTs1.begin(), validPMTs1.end(), pmtId)) continue;

            // Get First PEtime
            if (pmtTruthTimes1[pmtId].first == 0)
            {
                pmtTruthTimes1[pmtId].first = pe.HitPosInWindow;
            }
            else {
                // Get Double PETime
                if (pe.HitPosInWindow < pmtTruthTimes1[pmtId].first)
                {
                    pmtTruthTimes1[pmtId].second = pmtTruthTimes1[pmtId].first;
                    pmtTruthTimes1[pmtId].first = pe.HitPosInWindow;
                }
                else {
                    pmtTruthTimes1[pmtId].second = pe.HitPosInWindow;
                }
            }
        }

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
            UInt_t pmtId = ChannelId->at(j);

            //Select Double-charge PMTs
            std::unordered_set<int> validPMTsSet(validPMTs1.begin(), validPMTs1.end());
            if (!validPMTsSet.count(pmtId)) continue;
            if (pmtTruthTimes1[pmtId].second == 0) continue;

			char pmtname[1000];
            char pmtname2[1000];
			sprintf(pmtname,"PMT_%d Waveform",PMTId);
            sprintf(pmtname2, "PMT_%d Waveform;Time(ns);Voltage(mV)", PMTId);
			TH1F* h1 = new TH1F(pmtname,pmtname2, windowSize, 0, windowSize);
			
			for(Int_t k=1; k<=windowSize; k++)
			{
				h1->SetBinContent(k, Waveform->at(k-1+windowSize*j));
			}

			h0->Add(h1);
            std::cout << "Trigger" << i << "PMT" << j << std::endl;
		}

		f1->cd("../");
	
	}


	f1->Write();
	f1->Close();
    delete readoutChain;
    delete simTriggerChain;

}