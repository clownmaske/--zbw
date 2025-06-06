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

struct PeakInfo
{
    int position;
    double height;
};

//生成函数
double generateRandomTime(TRandom3& rng) {
    double u = rng.Uniform();
    double t = 0.0;
    double tolerance = 1e-8;
    int maxIter = 1000;
    for (int iter = 0; iter < maxIter; ++iter) {
        double term1 = 26.7 * (1 - exp(-t / 26.7));
        double term2 = (26.7 * 1.67) / (26.7 + 1.67) * (1 - exp(-t / 26.7 - t / 1.67));
        double cdf = (26.7 + 1.67) / (26.7 * 26.7) * (term1 - term2);
        double f = cdf - u;
        if (std::abs(f) < tolerance) {
            break;
        }
        // 计算导数即n(t)
        double n_t = (26.7 + 1.67) / (26.7 * 26.7) * (exp(-t / 26.7) - exp(-t * (1.0 / 1.67 + 1.0 / 26.7)));
        if (n_t < 1e-10)
        {
            // 防止除零
            t = (u > cdf) ? t + 1.0 : 0.0;
            continue;
        }
        t -= f / n_t;
        if (t < 0) { t = 0; break; }
    }
    return t;
}


void TimeCollection3(TString baselinefile, TString outputfile)
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

    //Create New Tree
    TTree* WaveformAnalysis = new TTree("WaveformAnalysis", "WaveformAnalysis");
    std::vector<UInt_t> TriggerVec1;
    std::vector<UInt_t> PMTIdVec1;
    std::vector<Double_t> TimepreVec1;
    //std::vector<UInt_t> PETypeVec1;
    WaveformAnalysis->Branch("Trigger", &TriggerVec1);
    WaveformAnalysis->Branch("PMTId", &PMTIdVec1);
    WaveformAnalysis->Branch("PETime", &TimepreVec1);
    //WaveformAnalysis->Branch("PEType", &PETypeVec1);

    //Get Single-charge PMTId
    std::vector<std::vector<int>> SinglechargePMTs;
    const Long64_t totalEntries = simTriggerChain->GetEntries();
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
    std::cout << "Single-charge PMTs Ready!" << endl;

    //Get Double-charge PMTId
    std::vector<std::vector<int>> DoublechargePMTs;
    //for (Int_t i = 0; i < t->GetEntries(); i++)
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

    const int minCrossCount = 3;
    const double amplitudeThreshold = 2.0;
    
    TH1F* hPre1 = new TH1F("hPre1", "Prediction Error;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre2 = new TH1F("hPre2", "Double Prediction Error1;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre3 = new TH1F("hPre3", "Double Prediction Error2;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre4 = new TH1F("hPre4", "Double Prediction Error3;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre41 = new TH1F("hPre41", "Double Prediction Error3;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre5 = new TH1F("hPre5", "Double Prediction Error4;PreTimeError(ns);Counts", 200, -10, 10);
    TH1F* hPre51 = new TH1F("hPre51", "Double Prediction Error4;PreTimeError(ns);Counts", 200, -10, 10);
    TH1F* hPre6 = new TH1F("hPre6", "First PETime;PreTimeError(ns);Counts", 500, 0, 500);
    TH1F* hPre7 = new TH1F("hPre7", "t_Sampling;t_Sampling(ns);Counts", 500, 0, 500);

    TH1F* hPre10 = new TH1F("hPre10", "Prediction Error;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre20 = new TH1F("hPre20", "Double Prediction Error1;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre30 = new TH1F("hPre30", "Double Prediction Error2;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre40 = new TH1F("hPre40", "Double Prediction Error3;PreTimeError(ns);Counts", 100, -5, 5);
    TH1F* hPre50 = new TH1F("hPre50", "Double Prediction Error4;PreTimeError(ns);Counts", 200, -10, 10);
    

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

        TriggerVec1.clear();
        PMTIdVec1.clear();
        TimepreVec1.clear();
        //PETypeVec1.clear();



        std::unordered_map<int, int> pmtCounterCurrent;
        for (const auto& pe : *PEList)
        {
            ++pmtCounterCurrent[pe.PMTId];
        }

        // Only get Single-charge PMTs
        const auto& singlePMTs = SinglechargePMTs[i];
        std::unordered_set<UInt_t> validPMTs;
        for (auto id : singlePMTs)
        {
            validPMTs.insert(static_cast<UInt_t>(id));
        }

        std::unordered_map<UInt_t, Double_t> pmtTruthTimes;
        for (const auto& pe : *PEList)
        {
            if (std::count(validPMTs.begin(), validPMTs.end(), pe.PMTId) > 0)
            {
                pmtTruthTimes[pe.PMTId] = pe.HitPosInWindow;
            }
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

            //const Double_t rt = +8.1201e-11 * integral * integral * integral * integral
                //- 1.46696e-07 * integral * integral * integral
                //+ 8.91729e-05 * integral * integral
                //- 0.0246859 * integral
                //+ 6.50041;
            const Double_t rt = 1.63502e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                - 0.000252396 * Peak_wave * Peak_wave * Peak_wave
                + 0.0142514 * Peak_wave * Peak_wave
                - 0.370663 * Peak_wave
                + 8.3774;
            if (rt < 0)continue;
            const Double_t p1 = t_1 - rt;
            const Double_t p10 = t_1 - 4.874;
            //const Double_t p2 = t_2 - rt - pw - tTruth;
            if (p1 < 200)continue;

            // Fill the Vector
            TriggerVec1.push_back(i);
            PMTIdVec1.push_back(PMTId);
            TimepreVec1.push_back(p1);
            
            hPre1->Fill(p1 - tTruth);
            hPre10->Fill(p10 - tTruth);

            cout << "Trigger" << i << "|ChannelId" << PMTId << "|Pretime|" << p1 << std::endl;

        }

        for (size_t j = 0; j < ChannelId->size(); ++j)
        {
            const UInt_t pmtId = ChannelId->at(j);
            Double_t baseline = Baseline->at(j);
            Double_t t1_threshold = baseline - 5;
            Double_t t2_threshold = baseline - 5;
            Double_t t_1;
            Double_t t_2;

            //Select Double-charge PMTs
            std::unordered_set<int> validPMTsSet(validPMTs1.begin(), validPMTs1.end());
            if (!validPMTsSet.count(pmtId)) continue;
            if (pmtTruthTimes1[pmtId].second == 0) continue;

            //Get tTruth
            const auto& [t1Truth, t2Truth] = pmtTruthTimes1[pmtId];

            //Get Waveform
            Int_t Minima_wave = 1000;
            std::vector<UInt_t> waveformData;
            const Int_t offset = windowSize * j;

            for (Int_t k = 0; k < windowSize; ++k)
            {
                const UInt_t val = Waveform->at(offset + k);
                if (val < Minima_wave) Minima_wave = val;
                waveformData.push_back(Waveform->at(offset + k));
            }
            Double_t Peak_wave = baseline - Minima_wave;

            //Find Risetime
            auto findRisingEdges = [](const std::vector<UInt_t>& data,
                Double_t threshold,
                int numEdges,
                int requiredCrossCount)
            {
                std::vector<int> edges;
                int crossCount = 0;
                bool inEdge = false;

                for (size_t t = 0; t < data.size(); ++t) {
                    if (data[t] < threshold) {
                        crossCount++;
                        if (!inEdge && crossCount >= requiredCrossCount) {
                            edges.push_back(t - requiredCrossCount + 1);
                            inEdge = true;
                            crossCount = 0;
                            if (edges.size() >= numEdges) break;
                        }
                    }
                    else {
                        inEdge = false;
                        crossCount = 0;
                    }
                }
                return edges;
            };

            auto risingEdges = findRisingEdges(waveformData, baseline - 5, 2, 3);

            if (risingEdges.size() == 2)
            {
                int raw_t1 = -1, raw_t2 = -1;
                raw_t1 = risingEdges[0];
                raw_t2 = risingEdges[1];

                if (raw_t1 == -1 || raw_t2 == -1) continue;
                if (raw_t2 <= raw_t1) continue;

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
                t_2 = interpolateTime(waveformData, raw_t2, t2_threshold, false);

                // 峰值检测验证有效性
                std::vector<PeakInfo> peaks;
                bool inPeak = false;
                int peakStart = 0;
                double currentMax = 0.0;
                int maxPos = 0;

                for (int t = 0; t < windowSize; ++t) {
                    double amplitude = baseline - waveformData[t];

                    if (amplitude > amplitudeThreshold) {
                        if (!inPeak) {
                            inPeak = true;
                            peakStart = t;
                            currentMax = amplitude;
                            maxPos = t;
                        }
                        else {
                            if (amplitude > currentMax) {
                                currentMax = amplitude;
                                maxPos = t;
                            }
                        }
                    }
                    else {
                        if (inPeak) {
                            inPeak = false;
                            if ((t - peakStart) >= minCrossCount) {
                                peaks.push_back({ maxPos, currentMax });
                            }
                        }
                    }
                }

                // 处理最后一个峰
                if (inPeak && (windowSize - peakStart) >= minCrossCount) {
                    peaks.push_back({ maxPos, currentMax });
                }

                if (peaks.size() == 1) continue;
                if (peaks.size() == 2)
                {

                    PeakInfo SecPeak = peaks[0];
                    const Double_t secpeak = SecPeak.height;
                    const Double_t deltat2 = 1.63502e-06 * secpeak * secpeak * secpeak * secpeak
                        - 0.000252396 * secpeak * secpeak * secpeak
                        + 0.0142514 * secpeak * secpeak
                        - 0.370663 * secpeak
                        + 8.3774;
                    if (deltat2 < 0)continue;
                    PeakInfo SecPeak1 = peaks[1];
                    const Double_t secpeak1 = SecPeak1.height;
                    const Double_t deltat3 = 1.63502e-06 * secpeak1 * secpeak1 * secpeak1 * secpeak1
                        - 0.000252396 * secpeak1 * secpeak1 * secpeak1
                        + 0.0142514 * secpeak1 * secpeak1
                        - 0.370663 * secpeak1
                        + 8.3774;
                    if (deltat3 < 0)continue;
                    if (SecPeak.position < t_1 || SecPeak.position > t_1 + 5)continue;
                    if (SecPeak1.position < t_2 || SecPeak1.position > t_2 + 5)continue;

                    const Double_t p2 = t_1 - deltat2;
                    const Double_t p3 = t_2 - deltat3;
                    const Double_t p20 = t_1 - 4.874;
                    const Double_t p30 = t_2 - 4.874;


                    if (p2 < 200)continue;
                    if (p3 < 200)continue;



                    // Fill the Vector
                    TriggerVec1.push_back(i);
                    TriggerVec1.push_back(i);
                    PMTIdVec1.push_back(pmtId);
                    PMTIdVec1.push_back(pmtId);
                    TimepreVec1.push_back(p2);
                    TimepreVec1.push_back(p3);

                    hPre2->Fill(p2 - t1Truth);
                    hPre3->Fill(p3 - t2Truth);
                    hPre20->Fill(p20 - t1Truth);
                    hPre30->Fill(p30 - t2Truth);

                    cout << "Trigger" << i << "|ChannelId" << pmtId << "|p2" << p2 << "|p3" << p3 << std::endl;
                };
            };
            if (risingEdges.size() == 1)
            {
                //findCrossing Function
                auto findCrossing = [minCrossCount](const std::vector<UInt_t>& data,
                    Double_t threshold, bool forward)
                {
                    int crossCount = 0;
                    int foundIdx = -1;
                    int step = forward ? 1 : -1;
                    int start = forward ? 0 : data.size() - 1;
                    int end = forward ? data.size() : -1;

                    for (int t = start; t != end; t += step)
                    {
                        if (data[t] < threshold)
                        {
                            if (++crossCount >= minCrossCount)
                            {
                                foundIdx = forward ?
                                    (t - minCrossCount + 1) :
                                    (t + minCrossCount - 1);
                                break;
                            }
                        }
                        else
                        {
                            crossCount = 0;
                        }
                    }
                    return foundIdx;
                };
                int raw_t1 = findCrossing(waveformData, t1_threshold, true);
                int raw_t2 = findCrossing(waveformData, t2_threshold, false);

                if (raw_t1 == -1 || raw_t2 == -1 || raw_t2 <= raw_t1) continue;

                if (raw_t1 < 180)continue;
                if (raw_t2 < 180)continue;

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

                //峰值寻找
                std::vector<PeakInfo> peaks;
                bool inPeak = false;
                int peakStart = 0;
                double currentMax = 0.0;
                int maxPos = 0;

                for (int t = 0; t < windowSize; ++t) {
                    double amplitude = baseline - waveformData[t];

                    if (amplitude > amplitudeThreshold) {
                        if (!inPeak) {
                            inPeak = true;
                            peakStart = t;
                            currentMax = amplitude;
                            maxPos = t;
                        }
                        else {
                            if (amplitude > currentMax) {
                                currentMax = amplitude;
                                maxPos = t;
                            }
                        }
                    }
                    else {
                        if (inPeak) {
                            inPeak = false;
                            if ((t - peakStart) >= minCrossCount) {
                                peaks.push_back({ maxPos, currentMax });
                            }
                        }
                    }
                }

                // 处理最后一个峰
                if (inPeak && (windowSize - peakStart) >= minCrossCount) {
                    peaks.push_back({ maxPos, currentMax });
                }



                Double_t p4;
                Double_t p5;
                if (peaks.size() == 1)
                {
                    const Double_t pw = -3.60278e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                        + 0.000573038 * Peak_wave * Peak_wave * Peak_wave
                        - 0.034562 * Peak_wave * Peak_wave
                        + 1.06207 * Peak_wave
                        - 3.15224;
                    const Double_t rt = 1.63502e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                        - 0.000252396 * Peak_wave * Peak_wave * Peak_wave
                        + 0.0142514 * Peak_wave * Peak_wave
                        - 0.370663 * Peak_wave
                        + 8.3774;
                    if (rt < 0)continue;
                    p4 = t_1 - rt;
                    p5 = t_2 - pw - rt;
                    hPre4->Fill(p4 - t1Truth);
                    hPre5->Fill(p5 - t2Truth);
                };
                if (peaks.size() == 2)
                {
                    PeakInfo SecPeak = peaks[0];
                    const Double_t secpeak = SecPeak.height;
                    const Double_t rt1 = 1.63502e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                        - 0.000252396 * Peak_wave * Peak_wave * Peak_wave
                        + 0.0142514 * Peak_wave * Peak_wave
                        - 0.370663 * Peak_wave
                        + 8.3774;
                    if (rt1 < 0)continue;
                    if (SecPeak.position < t_1 || SecPeak.position > t_1 + 5)continue;

                    PeakInfo SecPeak1 = peaks[1];
                    const Double_t secpeak1 = SecPeak1.height;
                    const Double_t pw = -3.60278e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                        + 0.000573038 * Peak_wave * Peak_wave * Peak_wave
                        - 0.034562 * Peak_wave * Peak_wave
                        + 1.06207 * Peak_wave
                        - 3.15224;
                    const Double_t rt2 = 1.63502e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                        - 0.000252396 * Peak_wave * Peak_wave * Peak_wave
                        + 0.0142514 * Peak_wave * Peak_wave
                        - 0.370663 * Peak_wave
                        + 8.3774;
                    if (rt2 < 0)continue;
                    if (SecPeak1.position < t_2 || SecPeak1.position > t_2 + 5)continue;
                    
                    p4 = t_1 - rt1;
                    p5 = t_2 - pw - rt2;
                    hPre41->Fill(p4 - t1Truth);
                    hPre51->Fill(p5 - t2Truth);
                };

                
                const Double_t p40 = t_1 - 4.874;
                const Double_t p50 = t_2 - 8.716 - 4.874;
                // Fill the Vector
                TriggerVec1.push_back(i);
                TriggerVec1.push_back(i);
                PMTIdVec1.push_back(pmtId);
                PMTIdVec1.push_back(pmtId);
                TimepreVec1.push_back(p4);
                TimepreVec1.push_back(p5);

                
                hPre40->Fill(p40 - t1Truth);
                hPre50->Fill(p50 - t2Truth);

      

                cout << "Trigger" << i << "|ChannelId" << pmtId << "|p4" << p4 << "|p5" << p5 << std::endl;



            };
        }

        // ====================== 处理三个及以上光电子 ======================
        for (size_t j = 0; j < ChannelId->size(); ++j)
        {
            const UInt_t pmtId = ChannelId->at(j);
            Double_t baseline = Baseline->at(j);

            // 获取当前PMT的光电子数目
            auto it = pmtCounterCurrent.find(pmtId);
            if (it == pmtCounterCurrent.end() || it->second < 3) continue;
            int numPE = it->second;

            // 获取波形数据
            std::vector<UInt_t> waveformData;
            const Int_t offset = windowSize * j;

            for (Int_t k = 0; k < windowSize; ++k)
            {
                waveformData.push_back(Waveform->at(offset + k));
            }

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

            // 找到第一个交点t_1
            Double_t t1_threshold = baseline - 5;
            int raw_t1 = findCrossing(waveformData, t1_threshold, true);
            if (raw_t1 == -1) continue;
            if (raw_t1 < 200) continue;

            auto interpolateTime = [](const std::vector<UInt_t>& data,
                int raw_t, Double_t threshold, bool isRising)
            {
                if (raw_t <= 0 || raw_t >= data.size() - 1)
                    return static_cast<Double_t>(raw_t);

                const Double_t y0 = data[raw_t - (isRising ? 0 : 1)];
                const Double_t y1 = data[raw_t + (isRising ? 1 : 0)];
                return raw_t + (threshold - y0) / (y1 - y0);
            };

            Double_t t_1 = interpolateTime(waveformData, raw_t1, t1_threshold, false);

            const Double_t rt = 4.354;
            const Double_t timepre_start = t_1 - rt;
            // 填充第一个预测时间
            TriggerVec1.push_back(i);
            PMTIdVec1.push_back(pmtId);
            TimepreVec1.push_back(timepre_start);

            hPre6->Fill(timepre_start);

        
            cout << "Trigger" << i << "|ChannelId" << pmtId << "|Pretime_start" << timepre_start << std::endl;


            // 生成n-1个随机时间
            TRandom3 rng(std::chrono::system_clock::now().time_since_epoch().count());
            for (int m = 0; m < numPE - 1; ++m) {
                double t_random = generateRandomTime(rng);
                Double_t timepre = timepre_start + t_random;
                if (timepre > 500)continue;

                //     填充到Tree中
                TriggerVec1.push_back(i);
                PMTIdVec1.push_back(pmtId);
                TimepreVec1.push_back(timepre);

                hPre7->Fill(timepre);

            }


        }
        // Fill the Tree
        WaveformAnalysis->Fill();
        
    }
    
    hPre1->Write();
    hPre2->Write();
    hPre3->Write();
    hPre4->Write();
    hPre5->Write();
    hPre41->Write();
    hPre51->Write();
    hPre6->Write();
    hPre7->Write();
    hPre10->Write();
    hPre20->Write();
    hPre30->Write();
    hPre40->Write();
    hPre50->Write();
    
    f1->Write();
    f1->Close();
    delete readoutChain;
    delete simTriggerChain;
    fb->Close();
}

