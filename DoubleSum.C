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

void DoubleSum(TString inputfile, TString baselinefile, TString outputfile) {
    TFile* f = new TFile(inputfile);
    TFile* fb = new TFile(baselinefile);
    TFile* f1 = new TFile(outputfile, "recreate");

    TTree* t = (TTree*)f->Get("Readout");
    TTree* t1 = (TTree*)f->Get("SimTriggerInfo");
    TTree* tb = (TTree*)fb->Get("BaselineTree");

    std::vector<UInt_t>* ChannelId = new std::vector<UInt_t>;
    std::vector<UInt_t>* Waveform = new std::vector<UInt_t>;
    vector<JPSimPE_t>* PEList = new vector<JPSimPE_t>;
    std::vector<Double_t>* Baseline = new std::vector<Double_t>;

    t->SetBranchAddress("ChannelId", &ChannelId);
    t->SetBranchAddress("Waveform", &Waveform);
    t1->SetBranchAddress("PEList", &PEList);
    tb->SetBranchAddress("Baseline", &Baseline);

    t->GetEntry(1);
    Int_t windowSize = Waveform->size() / ChannelId->size();

    std::vector<std::vector<int>> DoublechargePMTs;
    for (Int_t i = 0; i < 400; i++) {
        t1->GetEntry(i);
        std::unordered_map<int, int> pmtCounter;
        for (const auto& pe : *PEList) ++pmtCounter[pe.PMTId];

        std::vector<int> currentUniquePMTs;
        for (const auto& [pmtId, count] : pmtCounter) {
            if (count == 2) currentUniquePMTs.push_back(pmtId);
        }
        DoublechargePMTs.push_back(std::move(currentUniquePMTs));
    }

    // 扩展的波形存储结构
    std::map<int, std::vector<Double_t>> sumWaveforms;
    std::map<int, int> waveformCounts;
    const int preSamples = 50;
    const int postSamples = 150;

    for (Int_t i = 0; i < 400; i++) {
        t->GetEntry(i);
        tb->GetEntry(i);
        t1->GetEntry(i);

        const auto& validPMTs = DoublechargePMTs[i];
        std::unordered_map<UInt_t, std::pair<Double_t, Double_t>> pmtTruthTimes;
        std::unordered_map<UInt_t, int> pmtCounter;

        for (const auto& pe : *PEList) ++pmtCounter[pe.PMTId];
        for (const auto& pe : *PEList) {
            if (pmtCounter[pe.PMTId] != 2) continue;
            if (!std::count(validPMTs.begin(), validPMTs.end(), pe.PMTId)) continue;

            auto& times = pmtTruthTimes[pe.PMTId];
            if (times.first == 0) times.first = pe.HitPosInWindow;
            else {
                if (pe.HitPosInWindow < times.first) {
                    times.second = times.first;
                    times.first = pe.HitPosInWindow;
                }
                else times.second = pe.HitPosInWindow;
            }
        }

        for (size_t j = 0; j < ChannelId->size(); ++j) {
            const UInt_t pmtId = ChannelId->at(j);
            if (!std::count(validPMTs.begin(), validPMTs.end(), pmtId)) continue;
            if (pmtTruthTimes[pmtId].second == 0) continue;

            const auto& [t1Truth, t2Truth] = pmtTruthTimes[pmtId];
            Double_t deltaT = t2Truth - t1Truth;

            // 寻找匹配的时间差区间
            int targetK = -1;
            for (int k = 0; k <= 50; ++k) {
                if (deltaT >= (k - 0.2) && deltaT < (k + 0.2)) {
                    targetK = k;
                    break;
                }
            }
            if (targetK == -1) continue;

            // 波形对齐处理
            int t1_idx = static_cast<int>(t1Truth);
            int start = t1_idx - preSamples;
            int end = t1_idx + postSamples;
            if (start < 0 || end >= windowSize) continue;

            // 提取并叠加波形
            std::vector<Double_t> alignedWave;
            const Int_t offset = windowSize * j;
            for (int k = start; k <= end; ++k) {
                alignedWave.push_back(Waveform->at(offset + k) - Baseline->at(j));
            }

            if (sumWaveforms.find(targetK) == sumWaveforms.end()) {
                sumWaveforms[targetK] = std::vector<Double_t>(alignedWave.size(), 0.0);
            }
            for (size_t n = 0; n < alignedWave.size(); ++n) {
                sumWaveforms[targetK][n] += alignedWave[n];
            }
            waveformCounts[targetK]++;
            std::cout << "Trigger" << i << "ChannlId" << j << " over!" << endl;
        }
    }

    // 保存平均波形
    f1->cd();
    for (const auto& [k, waveSum] : sumWaveforms) {
        if (waveformCounts[k] == 0) continue;

        TH1F* hWave = new TH1F(
            Form("hWaveDeltaT%d", k),
            Form("Average Waveform #DeltaT=%.1f-%.1f;t(ns);voltage(mV)", k - 0.2, k + 0.2),
            waveSum.size(), 0, waveSum.size()
        );

        for (size_t i = 0; i < waveSum.size(); ++i) {
            hWave->SetBinContent(i + 1, waveSum[i] / waveformCounts[k]);
        }
        hWave->Write();
    }

    f1->Close();
    f->Close();
    fb->Close();
}