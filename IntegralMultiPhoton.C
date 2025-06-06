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
#include "TCanvas.h"
#include "TLegend.h"
#include "JPSimOutput.hh"

void IntegralMultiPhoton(TString baselinefile, TString outputfile) 
{  
    // ====================== 多文件处理设置 ======================
    const TString input_dir = "/home/zbw/Desktop/Data/";  // 修改为实际路径
    const int first_file_num = 0;
    const int last_file_num = 25;

    // 关键修改3：创建TChain处理多文件
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
    }

    // ====================== 初始化基线文件 ======================
    TFile* fb = new TFile(baselinefile);
    TTree* tb = (TTree*)fb->Get("BaselineTree");
    TFile* f1 = new TFile(outputfile, "recreate");

    // ====================== 设置分支地址 ======================
    std::vector<UInt_t>* ChannelId = nullptr;
    std::vector<UInt_t>* Waveform = nullptr;
    vector<JPSimPE_t>* PEList = nullptr;
    std::vector<Double_t>* Baseline = nullptr;

    readoutChain->SetBranchAddress("ChannelId", &ChannelId);
    readoutChain->SetBranchAddress("Waveform", &Waveform);
    simTriggerChain->SetBranchAddress("PEList", &PEList);
    tb->SetBranchAddress("Baseline", &Baseline);

    // ====================== 获取动态窗口大小 ======================
    readoutChain->GetEntry(0);
    const Int_t windowSize = Waveform->size() / ChannelId->size();
    std::cout << "WindowSize: " << windowSize << std::endl;

    // Create histograms
    std::map<int, TH1F*> rawHists;
    std::map<int, TH1F*> normHists; 
    const int bins = 100;
    const float xmin = 0, xmax = 2000;
    for (int n = 1; n <= 5; ++n) 
    {
        TString rawName = Form("h%dphoton_raw", n);
        TString rawTitle = Form("Raw Integral (%d photons);Integral(mV*ns);Proportion", n);
        rawHists[n] = new TH1F(rawName, rawTitle, bins, xmin, xmax);

        TString normName = Form("h%dphoton_norm", n);
        TString normTitle = Form("Normalized (%d photons)", n);
        normHists[n] = new TH1F(normName, normTitle, bins, xmin, xmax);
    }

    // PMT Count

    const Long64_t totalEntries = simTriggerChain->GetEntries();
    std::vector<std::map<int, std::vector<int>>> pmtGroupsList(totalEntries);

    for (Long64_t i = 0; i < totalEntries; ++i) 
    { 
        simTriggerChain->GetEntry(i);

        std::unordered_map<int, int> pmtCounter;
        for (const auto& pe : *PEList) {
            ++pmtCounter[pe.PMTId];
        }

        std::map<int, std::vector<int>> currentGroups;
        for (const auto& [pmtId, count] : pmtCounter) {
            if (count >= 1 && count <= 5) {
                currentGroups[count].push_back(pmtId);
            }
        }
        pmtGroupsList[i] = std::move(currentGroups);
    }
    std::cout << "PMTs Ready! Total events: " << totalEntries << std::endl;

    //For All Trigger
    const Long64_t nEntries = readoutChain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) 
    { 
        if (i == 3497)continue;
        if (i == 5797)continue;
        if (i == 8246)continue;
        if (i == 9910)continue;

        readoutChain->GetEntry(i);
        tb->GetEntry(i);
        simTriggerChain->GetEntry(i);

        const auto& pmtGroups = pmtGroupsList[i];

        for (int n = 1; n <= 5; ++n) 
        {
            auto it = pmtGroups.find(n);
            if (it == pmtGroups.end()) continue;

            std::unordered_set<UInt_t> targetPMTs;
            for (auto id : it->second) 
            {
                targetPMTs.insert(static_cast<UInt_t>(id));
            }

            for (size_t j = 0; j < ChannelId->size(); ++j) 
            {
                UInt_t pmtId = ChannelId->at(j);
                if (!targetPMTs.count(pmtId)) continue;

                // Integral caculate
                Double_t baseline = Baseline->at(j);
                double integral = 0.0;
                const Int_t offset = windowSize * j;
                for (Int_t k = 0; k < windowSize; ++k) {
                    integral += (baseline - Waveform->at(offset + k));
                }

                rawHists[n]->Fill(integral);
            }
            std::cout << "PMTs of " << i << "_" << n << " Finish!" << endl;
        }
    }

    int colors[] = { kBlack, kRed, kBlue, kGreen + 2, kMagenta };
    for (int n = 1; n <= 5; ++n) 
    {
        rawHists[n]->SetLineColor(colors[n - 1]);
        rawHists[n]->SetLineWidth(2);

        normHists[n] = (TH1F*)rawHists[n]->Clone(normHists[n]->GetName());
        if (normHists[n]->Integral() > 0) 
        {
            normHists[n]->Scale(1.0 / normHists[n]->Integral());
        }
        normHists[n]->SetLineColor(colors[n - 1]);
        normHists[n]->SetLineWidth(2);
    }

    // Draw Combine
    TCanvas* cRaw = new TCanvas("cRaw", "Raw Distributions", 1200, 800);
    TCanvas* cNorm = new TCanvas("cNorm", "Normalized Distributions", 1200, 800);

    cRaw->cd();
    TLegend* legRaw = new TLegend(0.7, 0.7, 0.9, 0.9);
    bool firstRaw = true;
    for (int n = 1; n <= 5; ++n) 
    {
        if (firstRaw)
        {
            rawHists[n]->Draw("HIST");
            firstRaw = false;
        }
        else 
        {
            rawHists[n]->Draw("HIST SAME");
        }
        legRaw->AddEntry(rawHists[n], Form("%d photons", n), "l");
    }
    legRaw->Draw();

    cNorm->cd();
    TLegend* legNorm = new TLegend(0.7, 0.7, 0.9, 0.9);
    bool firstNorm = true;
    for (int n = 1; n <= 5; ++n) {
        if (firstNorm) {
            normHists[n]->Draw("HIST");
            firstNorm = false;
        }
        else {
            normHists[n]->Draw("HIST SAME");
        }
        legNorm->AddEntry(normHists[n], Form("%d photons", n), "l");
    }
    legNorm->Draw();

    f1->cd();
    for (int n = 1; n <= 5; ++n) 
    {
        rawHists[n]->Write();
        normHists[n]->Write();
    }
    cRaw->Write();
    cNorm->Write();

    f1->Close();
    delete readoutChain;
    delete simTriggerChain;
    fb->Close();

}