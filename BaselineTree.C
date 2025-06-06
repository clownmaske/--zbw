#include <vector>
#include <map>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h" 
#include "TSystem.h"
#include "TH1F.h"
#include "TROOT.h"
//#include "TCanvas.h"
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

void BaselineTree(TString outputfile)
{
    setJPStyle();

    // 打开输入文件
    TChain* chain = new TChain("Readout");
    for (int i = 0; i <= 25; ++i)
    {
        TString filename;
        if (i == 0)
        {
            filename = "/home/zbw/Desktop/Data/SlowLS.root";     // 第一个文件无序号
        }
        else
        {
            filename = Form("/home/zbw/Desktop/Data/SlowLS_%d.root", i);  // SlowLS_1.root 到 SlowLS_25.root
        }
        chain->Add(filename);
    }

    // 创建输出文件
    TFile f1(outputfile, "recreate");

    // 设置分支地址
    std::vector<UInt_t> ChannelId, * pChannelId = &ChannelId;
    std::vector<UInt_t> Waveform, * pWaveform = &Waveform;
    chain->SetBranchAddress("ChannelId", &pChannelId);
    chain->SetBranchAddress("Waveform", &pWaveform);

    // 获取窗口大小
    chain->GetEntry(1);
    const Int_t windowSize = pWaveform->size() / pChannelId->size();
    std::cout << "WindowSize is " << windowSize << std::endl;

    // 创建输出树
    TTree baselineTree("BaselineTree", "BaselineTree");
    Int_t Trigger;
    std::vector<UInt_t> ChannelIdVec;
    std::vector<Double_t> BaselineVec;
    baselineTree.Branch("Trigger", &Trigger);
    baselineTree.Branch("ChannelId", &ChannelIdVec);
    baselineTree.Branch("Baseline", &BaselineVec);

    // 重用内存的临时变量
    std::vector<bool> is_pulse(windowSize, false);
    std::vector<UInt_t> noise;
    noise.reserve(windowSize);

    TH1F* hNoise = new TH1F("hNoise",
        "PMT Noise Distribution;Voltage of Noise(mV);Count",
        16,-8,8);

    // 处理事件循环
    const Long64_t nEntries = chain->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i)
    {
        chain->GetEntry(i);
        Trigger = i;

        const size_t nChannels = pChannelId->size();
        ChannelIdVec.clear();
        BaselineVec.clear();
        ChannelIdVec.reserve(nChannels);
        BaselineVec.reserve(nChannels);

        for (size_t j = 0; j < nChannels; ++j)
        {
            const UInt_t PMTId = (*pChannelId)[j];
            const UInt_t* waveform = &(*pWaveform)[j * windowSize];

            // 重置标记数组
            std::fill(is_pulse.begin(), is_pulse.end(), false);

            // 标记脉冲区域
            for (Int_t m = 0; m < windowSize; ++m)
            {
                if (waveform[m] < 916)
                {
                    const Int_t start = std::max(m - 20, 0);
                    const Int_t end = std::min(m + 20, windowSize - 1);
                    for (Int_t mm = start; mm <= end; ++mm)
                    {
                        is_pulse[mm] = true;
                    }
                }
            }

            // 收集噪声数据
            noise.clear();
            for (Int_t m = 0; m < windowSize; ++m)
            {
                if (!is_pulse[m])
                {
                    noise.push_back(waveform[m]);
                }
            }

            // 计算基线
            Double_t baseline = 921.1;  // 默认值
            if (!noise.empty())
            {
                baseline = 0.0;
                for (const auto& val : noise)
                {
                    baseline += val;
                }
                baseline /= noise.size();
                for (const auto& vall : noise)
                {
                    hNoise->Fill(vall - baseline);
                }
            }
            else
            {
                std::cerr << "\n\nError：Trigger " << i
                    << " ChannelId " << PMTId
                    << " No noise data！\n"
                    << " Push Enter to continue ..."
                    << std::endl;

                // 清除输入缓冲区
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
                // 等待用户输入
                std::cin.get();
            }

            ChannelIdVec.push_back(PMTId);
            BaselineVec.push_back(baseline);
            std::cout << "Trigger" << i << "ChannelId" << j << "_" << baseline << std::endl;
        }

        baselineTree.Fill();
    }

    // 写入输出并清理
    hNoise->Write();
    f1.Write();
    f1.Close();
    delete chain;
}