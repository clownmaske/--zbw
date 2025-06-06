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

void FunctionSingle3(TString baselinefile, TString outputfile)
{
    setJPStyle();
    //Open baselineputfile
    TFile* fb = new TFile(baselinefile);
    TTree* tb = (TTree*)fb->Get("BaselineTree");

    //Create a file to save the waveform
    TFile* f1 = new TFile(outputfile, "recreate");

    const TString input_dir = "/home/zbw/Desktop/Data/";                     
    const int first_file_num = 0;                       
    const int last_file_num = 25;                       
    TChain* readoutChain = new TChain("Readout");       
    TChain* simTriggerChain = new TChain("SimTriggerInfo"); 

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

    //Claim vector and struct
    std::vector<UInt_t>* ChannelId = new std::vector<UInt_t>;
    std::vector<UInt_t>* Waveform = new std::vector<UInt_t>;
    //vector<JPSimTruthTree_t>* truthList = new vector<JPSimTruthTree_t>;
    Int_t PEList_ = 0;
    vector<JPSimPE_t>* PEList = new vector<JPSimPE_t>;
    std::vector<Double_t>* Baseline = new std::vector<Double_t>;

    //SetBranchAddress
    readoutChain->SetBranchAddress("ChannelId", &ChannelId);
    readoutChain->SetBranchAddress("Waveform", &Waveform);
    simTriggerChain->SetBranchAddress("PEList", &PEList);
    //t1->SetBranchAddress("truthList", &truthList);
    tb->SetBranchAddress("Baseline", &Baseline);

    //Get windowSize
    readoutChain->GetEntry(1);
    //t1->GetEntry(0);
    Int_t windowSize = Waveform->size() / ChannelId->size();
    cout << "WindowSize is " << windowSize << endl;

    TH1F* h1 = new TH1F("SinglePeak", "SinglePeak;Voltage(mV);Counts", 100, 0, 100);

    const Long64_t nEntries = readoutChain->GetEntries();
    const Long64_t nSimTriggerEntries = simTriggerChain->GetEntries();
    const Long64_t nBaselineEntries = tb->GetEntries();

    if (nEntries != nSimTriggerEntries || nEntries != nBaselineEntries) {
        std::cerr << "Error!"
            << "\nReadout : " << nEntries
            << "\nSimTrigger : " << nSimTriggerEntries
            << "\nBaseline : " << nBaselineEntries
            << std::endl;
        return;
    }
    //Get PMTId
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
            if (count == 1) currentUniquePMTs.push_back(pmtId);
        }
        SinglechargePMTs.push_back(std::move(currentUniquePMTs));
    }
    std::cout << "Single-charge PMTs Ready!" << endl;

    //Create Histograms
    const int minCrossCount = 3;
    std::vector<Double_t> scatterX, scatterY;

    //For All Trigger
    for(Long64_t i = 0; i < nEntries; ++i)
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
            if (j >= Baseline->size()) {
                std::cerr << "Error " << i
                    << " nBaseline (" << Baseline->size()
                    << ")> ChannelId (" << ChannelId->size() << ")"
                    << std::endl;
                continue;
            }

            UInt_t PMTId = ChannelId->at(j);
            Double_t baseline = Baseline->at(j);
            Double_t t1_threshold = baseline - 5;
            Double_t t2_threshold = baseline - 5;
            Double_t t_1;
            Double_t t_2;

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

            // Fill the Histogram
            h1->Fill(Peak_wave);

            //Select Single-charge PMTs
            if (!validPMTs.count(PMTId)) continue;

            //Get tTruth
            const Double_t tTruth = pmtTruthTimes[PMTId];


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
            //const Double_t xValue = integral;
            const Double_t yValue = t_2 - t_1;
            //const Double_t yValue = t_1 - tTruth;
            scatterX.push_back(xValue);
            scatterY.push_back(yValue);

           
            cout << "Trigger" << i << " ChannelId" << PMTId << " Peak is " << Peak_wave << "|" << t_1 << "|" << tTruth << std::endl;
        }

    }

    // ==================== 新增分箱平均处理部分 ====================
    // 创建Profile直方图进行分箱平均
    const double xmin = *std::min_element(scatterX.begin(), scatterX.end());
    const double xmax = *std::max_element(scatterX.begin(), scatterX.end());
    const int nBins = 100; // 可根据数据量调整分箱数量

    TProfile* hprof = new TProfile("hprof", "Profile of tw vs Peak;Peak(mV);tw(ns)",
        nBins, xmin, xmax);
    for (size_t i = 0; i < scatterX.size(); ++i) {
        hprof->Fill(scatterX[i], scatterY[i]);
    }

    // 提取分箱平均数据
    std::vector<Double_t> avgX, avgY, avgYErr;
    for (int bin = 1; bin <= hprof->GetNbinsX(); ++bin) {
        if (hprof->GetBinEntries(bin) > 0) {
            avgX.push_back(hprof->GetBinCenter(bin));
            avgY.push_back(hprof->GetBinContent(bin));
            avgYErr.push_back(hprof->GetBinError(bin));
        }
    }

    // 创建带误差的TGraph
    TGraph* avgGraph = new TGraph(avgX.size(), avgX.data(), avgY.data());
    avgGraph->SetTitle("Binned Average Data");
    avgGraph->SetMarkerStyle(kFullCircle);
    avgGraph->SetMarkerColor(kRed);
    avgGraph->SetMarkerSize(1.2);
    avgGraph->SetLineWidth(2);

    // ==================== 多项式拟合部分 ====================
    // 创建并配置多项式函数（3次到7次）
    const int maxPolyOrder = 4;
    std::vector<TF1*> polyFunctions;
    std::vector<int> colors = { kBlue, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2, kYellow + 2, kRed };

    for (int order = 4; order <= maxPolyOrder; ++order) {
        TString formula = TString::Format("pol%d", order);
        TF1* f = new TF1(TString::Format("poly%d", order), formula, xmin, xmax);
        f->SetLineColor(colors[order - 3]);
        f->SetLineWidth(2);
        polyFunctions.push_back(f);
    }

    // 执行多阶多项式拟合
    std::vector<TFitResultPtr> fitResults;
    TString fitOption = "SQ";
    for (auto& f : polyFunctions) {
        fitResults.push_back(avgGraph->Fit(f, fitOption));
        fitOption += "+"; // 后续拟合叠加显示
    }

    // ==================== 增强可视化 ====================
    TCanvas* cFit = new TCanvas("cFit", "Fitting Results", 1920, 1080);
    cFit->SetGrid();


    // 子图1: 原始数据分布
    cFit->cd(1);
    TGraph* scatterPlot = new TGraph(scatterX.size(), scatterX.data(), scatterY.data());
    scatterPlot->SetTitle("Function of tw vs Peak;Peak(mV);tw(ns)");
    scatterPlot->SetMarkerStyle(kFullCircle);
    scatterPlot->SetMarkerColor(kBlue);
    scatterPlot->SetMarkerSize(0.3);
    scatterPlot->Draw("AP");
    avgGraph->Draw("P same");
    // 获取拟合参数
    double p0 = polyFunctions[0]->GetParameter(0);
    double p1 = polyFunctions[0]->GetParameter(1);
    double p2 = polyFunctions[0]->GetParameter(2);
    double p3 = polyFunctions[0]->GetParameter(3);
    double p4 = polyFunctions[0]->GetParameter(4);

    // 创建一个文字框显示参数
    TPaveText* paramBox = new TPaveText(0.15, 0.6, 0.5, 0.85, "NDC");  // NDC表示归一化坐标
    paramBox->SetFillColor(0);
    paramBox->SetBorderSize(1);
    paramBox->AddText(Form("p0 = %.6e", p0));
    paramBox->AddText(Form("p1 = %.6e", p1));
    paramBox->AddText(Form("p2 = %.6e", p2));
    paramBox->AddText(Form("p3 = %.6e", p3));
    paramBox->AddText(Form("p4 = %.6e", p4));
    paramBox->Draw();




    // ==================== 控制台输出增强 ====================
    std::cout << "\n========= Fit Summary =========\n";
    for (size_t i = 0; i < polyFunctions.size(); ++i) {
        std::cout << "\nPolynomial Order " << i + 3 << ":\n";
        std::cout << "Chi2/NDF: " << polyFunctions[i]->GetChisquare()
            << "/" << polyFunctions[i]->GetNDF()
            << " = " << polyFunctions[i]->GetChisquare() / polyFunctions[i]->GetNDF() << "\n";
        std::cout << "Parameters:\n";
        for (int p = 0; p <= i + 3; ++p) {
            std::cout << "  p" << p << " = " << polyFunctions[i]->GetParameter(p)
                << " - " << polyFunctions[i]->GetParError(p) << "\n";
        }
    }

    // ==================== 保存结果 ====================
    cFit->Write();
    avgGraph->Write("BinnedAverageData");
    scatterPlot->Write("scatterPlot");
    delete scatterPlot;
    h1->Write();
    f1->Write();
    f1->Close();
    delete readoutChain;
    delete simTriggerChain;
    fb->Close();
}