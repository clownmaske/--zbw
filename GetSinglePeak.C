#include <vector>
#include <map>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TGraph.h"
#include "JPSimOutput.hh" 
//#include "TNamed.h"
//#include "TCanvas.h"
//#include "JPWaveformPreprocess.h"
//#include "DictOutput_rdict.h"  // Include dictionary

void GetSinglePeak(TString inputfile, TString baselinefile, TString outputfile)
{
    //Open inputfile
    TFile* f = new TFile(inputfile);
    TFile* fb = new TFile(baselinefile);

    //Create a file to save the waveform
    TFile* f1 = new TFile(outputfile, "recreate");

    //Get the SimTriggerInfo and Readout Tree
    TTree* t = (TTree*)f->Get("Readout");
    TTree* t1 = (TTree*)f->Get("SimTriggerInfo");
    TTree* tb = (TTree*)fb->Get("BaselineTree");

    //Claim vector and struct
    std::vector<UInt_t>* ChannelId = new std::vector<UInt_t>;
    std::vector<UInt_t>* Waveform = new std::vector<UInt_t>;
    //vector<JPSimTruthTree_t>* truthList = new vector<JPSimTruthTree_t>;
    Int_t PEList_ = 0;
    vector<JPSimPE_t>* PEList = new vector<JPSimPE_t>;
    std::vector<Double_t>* Baseline = new std::vector<Double_t>;

    //SetBranchAddress
    t->SetBranchAddress("ChannelId", &ChannelId);
    t->SetBranchAddress("Waveform", &Waveform);
    t1->SetBranchAddress("PEList", &PEList);
    //t1->SetBranchAddress("truthList", &truthList);
    tb->SetBranchAddress("Baseline", &Baseline);

    //Get windowSize
    t->GetEntry(1);
    //t1->GetEntry(0);
    Int_t windowSize = Waveform->size() / ChannelId->size();
    cout << "WindowSize is " << windowSize << endl;

    TH1F* h1 = new TH1F("SinglePeak", "SinglePeak", 100, 0, 50);

    //Get PMTId
    std::vector<std::vector<int>> SinglechargePMTs;
    for (Int_t i = 0; i < 400; i++)
    {
        t1->GetEntry(i);
        PEList_ = PEList->size();

        //Count PMTs
        std::unordered_map<int, int> pmtCounter;
        for (const auto& pe : *PEList)
        {
            ++pmtCounter[pe.PMTId];
        }

        // Select single-charge PMTs
        std::vector<int> currentUniquePMTs;
        //for (const auto& entry : pmtCounter)
        for (const auto& [pmtId, count] : pmtCounter)
        {
            if (count == 1)
            {
                currentUniquePMTs.push_back(pmtId);
            }
        }
        SinglechargePMTs.push_back(std::move(currentUniquePMTs));

        //std::cout << "Entry " << i
            //<< " | Total PE: " << PEList_
            //<< " | Single-charge PMTs: " << SinglechargePMTs.back().size()
            //<< '\n' << std::endl;

    }

    std::cout << "Single-charge PMTs Ready!" << endl;

    //Create Histograms
    const int minCrossCount = 3;
    std::vector<Double_t> scatterX, scatterY;

    //For All Trigger
    for (Int_t i = 0; i < 400; i++)
    {
        t->GetEntry(i);
        tb->GetEntry(i);
        t1->GetEntry(i);

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
            Int_t t_1;
            Int_t t_2;

            //Select Single-charge PMTs
            if (!validPMTs.count(PMTId)) continue;

            //Get tTruth
            const Double_t tTruth = pmtTruthTimes[PMTId];

            //Peak of Waveform and Get t1 t2 and Draw
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

            const Double_t pw = t_2 - t_1;
            
            const Double_t xValue = Peak_wave;
            const Double_t yValue = t_2-t_1;
            scatterX.push_back(xValue);
            scatterY.push_back(yValue);

            // Fill the Histogram
            h1->Fill(Peak_wave);
            cout << "Trigger" << i << " ChannelId" << PMTId << " Peak is " << Peak_wave << std::endl;
        }

    }
    TGraph* scatterPlot = new TGraph(scatterX.size(), scatterX.data(), scatterY.data());
    //scatterPlot->GetYaxis()->SetRangeUser(-10.0, 10.0);
    scatterPlot->SetTitle("Peak-deltat Function;Peak(mV);t_2-t_1(ns)");
    scatterPlot->SetMarkerStyle(kFullCircle);
    scatterPlot->SetMarkerColor(kBlue);
    scatterPlot->SetMarkerSize(0.5);


    // 定义线性拟合函数
    TF1* f11 = new TF1("f1", "[0]*x + [1]");
    f11->SetLineColor(kRed);

    // 设置合理初始参数（根据数据范围估算）
    if (!scatterX.empty()) {
        f11->SetParameters(
            (scatterY.back() - scatterY[0]) / (scatterX.back() - scatterX[0]), // 初始斜率
            scatterY[0] // 初始截距
        );
    }

    // 执行拟合
    if (scatterX.size() > 1) {
        scatterPlot->Fit(f11, "Q"); // Q: 静默模式
        double a = f11->GetParameter(0);
        double b = f11->GetParameter(1);
        std::cout << "Fit result: y = (" << a << ")x + (" << b << ")" << endl;
    }


    f1->cd();
    scatterPlot->Write("scatterPlot");
    delete scatterPlot;
    h1->Write();
    f1->Write();
    f1->Close();
    f->Close();
    fb->Close();
}