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






void PMTNumberCount(TString baselinefile, TString outputfile)
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

    // ====================== 初始化基线文件 ======================
    TFile* fb = new TFile(baselinefile);
    TTree* tb = (TTree*)fb->Get("BaselineTree");

    // ====================== 设置分支地址 ======================
    std::vector<UInt_t>* ChannelId = nullptr;
    std::vector<UInt_t>* Waveform = nullptr;
    vector<JPSimPE_t>* PEList = nullptr;
    std::vector<Double_t>* Baseline = nullptr;

    readoutChain->SetBranchAddress("ChannelId", &ChannelId);
    readoutChain->SetBranchAddress("Waveform", &Waveform);
    simTriggerChain->SetBranchAddress("PEList", &PEList);
    tb->SetBranchAddress("Baseline", &Baseline);

    // ====================== 创建输出文件 ======================
    TFile* f1 = new TFile(outputfile, "recreate");


    std::unordered_map<int, int> photonCountStats; 

    const int maxPhotonCount = 20;
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

        std::unordered_map<int, int> currentPMTCounts;

        for (const auto& pe : *PEList)
        {
            currentPMTCounts[pe.PMTId]++;
        }

        for (const auto& [pmtId, count] : currentPMTCounts)
        {
            if (count > 0 && count <= maxPhotonCount) 
            {
                photonCountStats[count]++;
            }
        }
    }

    // 确定直方图范围
    int minCount = 1;
    int maxCount = maxPhotonCount;
    for (const auto& [count, freq] : photonCountStats) {
        if (count < minCount) minCount = count;
        if (count > maxCount) maxCount = count;
    }
    maxCount = std::min(maxCount, maxPhotonCount);

    // 创建直方图
    TH1F* hPhotonCounts = new TH1F("hPhotonCounts",
        "PMT Photon Count Distribution;Photon Count;Number of PMTs",
        maxCount - minCount + 1, minCount - 0.5, maxCount + 0.5);

    // 填充直方图
    for (int count = minCount; count <= maxCount; ++count) {
        if (photonCountStats.find(count) != photonCountStats.end()) {
            hPhotonCounts->SetBinContent(
                count - minCount + 1,
                photonCountStats[count]
            );
        }
    }

    // 设置显示选项
    hPhotonCounts->SetFillColor(kBlue);
    //hPhotonCounts->SetStats(true);

    TCanvas* c1 = new TCanvas("c1", "Canvas", 800, 600);

    // 保存结果
    f1->cd();
    hPhotonCounts->Write();
    hPhotonCounts->Draw();
    c1->SaveAs("PMTCount.pdf");

    // 打印统计摘要
    std::cout << "\n======== Photon Count Statistics ========\n";
    for (int count = minCount; count <= maxCount; ++count) {
        if (photonCountStats.count(count)) {
            std::cout << count << " photons: "
                << photonCountStats[count] << " PMTs\n";
        }
    }
    std::cout << "======================================\n";

    // 清理资源
    f1->Close();
    delete readoutChain;
    delete simTriggerChain;
}