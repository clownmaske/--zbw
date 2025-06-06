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
#include "TGraph.h"
#include "JPSimOutput.hh" 
#include "TFitResult.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TCanvas.h" 
//#include "TNamed.h"
//#include "TCanvas.h"
//#include "JPWaveformPreprocess.h"
//#include "DictOutput_rdict.h"  // Include dictionary


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

void SingleAnalysis(TString baselinefile, TString outputfile)
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

    // ====================== 获取窗口大小 ======================
    readoutChain->GetEntry(0);
    const Int_t windowSize = Waveform->size() / ChannelId->size();
    cout << "WindowSize is " << windowSize << endl;

    // ====================== 创建输出文件 ======================
    TFile* f1 = new TFile(outputfile, "recreate");

    // ====================== 初始化数据结构 ======================
    std::vector<std::vector<int>> SinglechargePMTs;
    const Long64_t totalEntries = simTriggerChain->GetEntries();  // 动态获取事件总数

    // ====================== PMT预处理 ======================
    for (Long64_t i = 0; i < totalEntries; ++i) {
        simTriggerChain->GetEntry(i);

        std::unordered_map<int, int> pmtCounter;
        for (const auto& pe : *PEList) {
            ++pmtCounter[pe.PMTId];
        }

        std::vector<int> currentUniquePMTs;
        for (const auto& [pmtId, count] : pmtCounter) {
            if (count == 1) {
                currentUniquePMTs.push_back(pmtId);
            }
        }
        SinglechargePMTs.push_back(std::move(currentUniquePMTs));
    }
    std::cout << "Single-charge PMTs Ready! Total events: " << totalEntries << endl;
    

    const int minCrossCount = 3;
    std::vector<double> riseTimes, pulseWidths;
    TH1F* hPeak = new TH1F("hPeak", " Single Peak Distribution;Voltage of Peak(mV);Count", 50, 0, 50);
    TH1F* hRiseTime = new TH1F("hRiseTime", " RiseTime Distribution;deltat(ns);Count", 100, 0, 10);
    TH1F* hPulseWidth = new TH1F("hPulseWidth", "PulseWidth Distribution;tw(ns);Count", 200, 0, 20);
    TH1F* hPre1 = new TH1F("hPre1", "Prediction Error;t_1-deltat-tTruth(ns);Count", 100, -5, 5);
    TH1F* hPre2 = new TH1F("hPre2", "t_2-tw-deltat-tTruth(ns);t_2-tw-deltat-tTruth(ns);Count", 100, -10, 10);

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

        const auto& singlePMTs = SinglechargePMTs[i];
        std::unordered_set<UInt_t> validPMTs;
        for (auto id : singlePMTs)
        {
            validPMTs.insert(static_cast<UInt_t>(id));
        }

        // Only get Single-charge PMTs
        std::unordered_map<UInt_t, Double_t> pmtTruthTimes;
        for (const auto& pe : *PEList)
        {
            if (std::count(validPMTs.begin(), validPMTs.end(), pe.PMTId) > 0)
            {
                pmtTruthTimes[pe.PMTId] = pe.HitPosInWindow;
            }
        }

        for (size_t j = 0; j < ChannelId->size(); j++)
        {
            UInt_t PMTId = ChannelId->at(j);
            Double_t baseline = Baseline->at(j);
            Double_t t1_threshold = baseline - 5;
            Double_t t2_threshold = baseline - 5;
            Double_t t_1;
            Double_t t_2;

            //Select Single-charge PMTs
            if (!validPMTs.count(PMTId)) continue;

            //Get tTruth
            const Double_t tTruth = pmtTruthTimes[PMTId];

            //Peak of Waveform and Get t1 t2 and Draw
            Int_t Minima_wave = 1000;
            std::vector<UInt_t> waveformData;
            const Int_t offset = windowSize * j;
            double integral = 0.0;

            for (Int_t k = 0; k < windowSize; ++k)
            {
                const UInt_t val = Waveform->at(offset + k);
                if (val < Minima_wave) Minima_wave = val;
                waveformData.push_back(Waveform->at(offset + k));
                integral += (baseline - Waveform->at(offset + k));
            }
            Double_t Peak_wave = baseline - Minima_wave;
            hPeak->Fill(Peak_wave);

            //findCrossing Function
            auto findCrossing = [minCrossCount](const std::vector<UInt_t>& waveformData, Double_t threshold, bool forward)
            {
                int crossCount = 0;
                int foundIdx = -1;
                int step = forward ? 1 : -1;
                int start = forward ? 0 : waveformData.size() - 1;
                int end = forward ? waveformData.size() : -1;

                for (int t = start; t != end; t += step) {
                    if (waveformData[t] < threshold) {
                        if (++crossCount >= minCrossCount) {
                            foundIdx = forward ?
                                (t - minCrossCount + 1) :
                                (t + minCrossCount - 1);
                            break;
                        }
                    }
                    else {
                        crossCount = 0;
                    }
                }
                return foundIdx;
            };

            int raw_t1 = findCrossing(waveformData, t1_threshold, true);
            int raw_t2 = findCrossing(waveformData, t2_threshold, false);


            auto interpolateTime = [](const std::vector<UInt_t>& data,
                int raw_t, Double_t threshold, bool isRising)
            {
                if (raw_t <= 0 || raw_t >= data.size() - 1)
                    return static_cast<Double_t>(raw_t);

                const Double_t y0 = data[raw_t - (isRising ? 0 : 1)];
                const Double_t y1 = data[raw_t + (isRising ? 1 : 0)];
                return raw_t + (threshold - y0) / (y1 - y0);
            };

            t_1 = interpolateTime(waveformData, raw_t1, t1_threshold, false);
            t_2 = interpolateTime(waveformData, raw_t2, t2_threshold, true);
            if (raw_t1 == -1 || raw_t2 == -1) continue;
            if (!(t_1 < t_2)) continue;

            //const Double_t rt = - 3.25595e-08 * Peak_wave * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                //+ 6.51398e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                //- 0.000505115 * Peak_wave * Peak_wave * Peak_wave
                //+ 0.019321 * Peak_wave * Peak_wave 
                //- 0.392735 * Peak_wave 
                //+ 7.67088;
            //const Double_t pw = -2.99672e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                //+ 0.00047288 * Peak_wave * Peak_wave * Peak_wave
                //- 0.0289518 * Peak_wave * Peak_wave
                //+ 0.935336 * Peak_wave
                //- 2.10849;

            const Double_t rt = + 8.1201e-11 * integral * integral * integral * integral
                - 1.46696e-07 * integral * integral * integral
                + 8.91729e-05 * integral * integral
                - 0.0246859 * integral
                + 6.50041;
            const Double_t pw = + 5.41956e-10 * integral * integral * integral * integral
                - 2.83408e-07 * integral * integral * integral
                - 5.80337e-05 * integral * integral
                + 0.0653895 * integral
                + 0.937227;

            const Double_t p1 = t_1 - 4.874 - tTruth;
            //const Double_t p1 = t_1 - rt - tTruth;
            const Double_t p2 = t_2 - rt - pw - tTruth;


            hRiseTime->Fill(t_1 - tTruth);
            hPulseWidth->Fill(t_2 - t_1);
            hPre1->Fill(p1);
            hPre2->Fill(p2);

            cout << "Trigger" << i << "|ChannelId" << PMTId << "|Risetime" << rt << "|PulseWidth" << pw << "|" << t_1 << "|" << tTruth << std::endl;

        }

    }
    hPeak->Write();
    hRiseTime->Write();
    hPulseWidth->Write();
    hPre1->Write();
    hPre2->Write();
    f1->Write();
    f1->Close();
    delete readoutChain;
    delete simTriggerChain;
    fb->Close();
}