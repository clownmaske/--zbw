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

struct PeakInfo
{
    int position;
    double height;
};

void DoubleAnalysis2(TString baselinefile, TString outputfile)
{
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

    //Get Double-charge PMTId
    std::vector<std::vector<int>> DoublechargePMTs;
    const Long64_t totalEntries = simTriggerChain->GetEntries();
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
    TH1F* hPre1 = new TH1F("hPre1", "t_1-deltat-tTruth(ns);t_1-deltat-tTruth(ns);Counts", 100, -10, 10);
    TH1F* hPre2 = new TH1F("hPre2", "t_2-deltat-tTruth(ns);t_2-deltat-tTruth(ns);Counts", 100, -10, 10);
    TH1F* hPre3 = new TH1F("hPre3", "t_1-deltat-tTruth(ns);t_1-deltat-tTruth(ns);Counts", 100, -10, 10);
    TH1F* hPre4 = new TH1F("hPre4", "t_2-tw-deltat-tTruth(ns);t_2-tw-deltat-tTruth(ns);Counts", 100, -10, 10);

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

        for (size_t j = 0; j < ChannelId->size(); ++j)
        {
            const UInt_t pmtId = ChannelId->at(j);
            Double_t baseline = Baseline->at(j);
            Double_t t1_threshold = baseline - 5;
            Double_t t2_threshold = baseline - 5;
            Int_t t_1;
            Int_t t_2;

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
                    const Double_t deltat2 = -3.25595e-08 * secpeak * secpeak * secpeak * secpeak * secpeak
                        + 6.51398e-06 * secpeak * secpeak * secpeak * secpeak
                        - 0.000505115 * secpeak * secpeak * secpeak
                        + 0.019321 * secpeak * secpeak
                        - 0.392735 * secpeak
                        + 7.67088;
                    PeakInfo SecPeak1 = peaks[1];
                    const Double_t secpeak1 = SecPeak1.height;
                    const Double_t deltat3 = -3.25595e-08 * secpeak1 * secpeak1 * secpeak1 * secpeak1 * secpeak1
                        + 6.51398e-06 * secpeak1 * secpeak1 * secpeak1 * secpeak1
                        - 0.000505115 * secpeak1 * secpeak1 * secpeak1
                        + 0.019321 * secpeak1 * secpeak1
                        - 0.392735 * secpeak1
                        + 7.67088;

                    const Double_t p2 = raw_t1 - deltat2 - t1Truth;
                    const Double_t p3 = raw_t2 - deltat3 - t2Truth;

                    hPre1->Fill(p2);
                    hPre2->Fill(p3);

                    cout << "Trigger" << i << "|ChannelId" << pmtId << "|p2|" << p2 << "|p3|" << p3 << std::endl;
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

                t_1 = raw_t1;
                t_2 = raw_t2;

                const Double_t pw = -2.99672e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                    + 0.00047288 * Peak_wave * Peak_wave * Peak_wave
                    - 0.0289518 * Peak_wave * Peak_wave
                    + 0.935336 * Peak_wave
                    - 2.10849;
                const Double_t rt = -3.25595e-08 * Peak_wave * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                    + 6.51398e-06 * Peak_wave * Peak_wave * Peak_wave * Peak_wave
                    - 0.000505115 * Peak_wave * Peak_wave * Peak_wave
                    + 0.019321 * Peak_wave * Peak_wave 
                    - 0.392735 * Peak_wave 
                    + 7.67088;

                const Double_t p4 = t_1 - rt - t1Truth;
                const Double_t p5 = t_2 - pw - rt - t2Truth;

                hPre3->Fill(p4);
                hPre4->Fill(p5);



                cout << "Trigger" << i << "|ChannelId" << pmtId << "|p4|" << p4 << "|p5|" << p5 << std::endl;

            };
            
        }
    }
    hPre1->Write();
    hPre2->Write();
    hPre3->Write();
    hPre4->Write();
    f1->Write();
    f1->Close();
    delete readoutChain;
    delete simTriggerChain;
    fb->Close();
}