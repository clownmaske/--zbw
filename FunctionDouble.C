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
#include "JPSimOutput.hh" 
#include "TFitResult.h"
#include "TLegend.h"
#include "TSpline.h"
#include "TCanvas.h"
//#include "JPWaveformPreprocess.h"
//#include "DictOutput_rdict.h"  // Include dictionary

void FunctionDouble(TString inputfile, TString baselinefile, TString outputfile)
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

	//Get PMTId
	std::vector<std::vector<int>> DoublechargePMTs;

	//for (Int_t i = 0; i < t->GetEntries(); i++)
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

		//std::cout << "Entry " << i
			//<< " | Total PE: " << PEList_
			//<< " | Double-charge PMTs: " << DoublechargePMTs.back().size()
			//<< '\n' << std::endl;
	}

	std::cout << "Double-charge PMTs Ready!" << endl;

	//Create Histograms
	const int minCrossCount = 3;
	TH1F* hT2_T2ruthDiff = new TH1F("hT2_T2Truth", "t222 - t2Truth(ns)", 100, -20, 20);
	std::vector<Double_t> scatterX, scatterY;

	//Get Waveform
	for (Int_t i = 0; i < 400; i++)
	{
		t->GetEntry(i);
		tb->GetEntry(i);
		t1->GetEntry(i);

		const auto& validPMTs = DoublechargePMTs[i];
		std::unordered_map<UInt_t, std::pair<Double_t, Double_t>> pmtTruthTimes;
		std::unordered_map<UInt_t, int> pmtCounter;

		// Get PE Number
		for (const auto& pe : *PEList)
		{
			++pmtCounter[pe.PMTId];
		}

		// Get Double PEs of Certain PMT
		for (const auto& pe : *PEList)
		{
			const UInt_t pmtId = pe.PMTId;

			// Only Get Double-charge PMTs
			if (pmtCounter[pmtId] != 2) continue;
			if (!std::count(validPMTs.begin(), validPMTs.end(), pmtId)) continue;

			// Get First PEtime
			if (pmtTruthTimes[pmtId].first == 0)
			{
				pmtTruthTimes[pmtId].first = pe.HitPosInWindow;
			}
			else {
				// Get Double PETime
				if (pe.HitPosInWindow < pmtTruthTimes[pmtId].first)
				{
					pmtTruthTimes[pmtId].second = pmtTruthTimes[pmtId].first;
					pmtTruthTimes[pmtId].first = pe.HitPosInWindow;
				}
				else {
					pmtTruthTimes[pmtId].second = pe.HitPosInWindow;
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
			std::unordered_set<int> validPMTsSet(validPMTs.begin(), validPMTs.end());
			if (!validPMTsSet.count(pmtId)) continue;
			if (pmtTruthTimes[pmtId].second == 0) continue;

			//Get tTruth
			const auto& [t1Truth, t2Truth] = pmtTruthTimes[pmtId];
			//if (t2Truth - t1Truth > 16|| t2Truth - t1Truth < 10)continue;

			//Get t1 t2 and Draw
			Int_t Minima_wave = 1000;
			double Integral_wave = 0.0;
			std::vector<UInt_t> waveformData;
			const Int_t offset = windowSize * j;
			for (Int_t k = 0; k < windowSize; ++k)
			{
				const UInt_t val = Waveform->at(offset + k);
				if (val < Minima_wave) Minima_wave = val;
				waveformData.push_back(Waveform->at(offset + k));
				Integral_wave += (baseline - Waveform->at(offset + k));
			}
			Double_t Peak_wave = baseline - Minima_wave;

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

			if ( t_2 - t_1 > 30 || t_2 - t_1 < 15 )continue;

			const Double_t pw = t_2 - t_1;
			const Double_t deltaT1 = 
				- 5.56892e-18 * Integral_wave * Integral_wave * Integral_wave * Integral_wave * Integral_wave * Integral_wave * Integral_wave
				+ 6.84601e-15 * Integral_wave * Integral_wave * Integral_wave * Integral_wave * Integral_wave * Integral_wave
				- 1.06188e-12 * Integral_wave * Integral_wave * Integral_wave * Integral_wave * Integral_wave
				- 2.06007e-9 * Integral_wave * Integral_wave * Integral_wave * Integral_wave
				+ 1.32727e-6 * Integral_wave * Integral_wave * Integral_wave
				- 0.000370489 * Integral_wave * Integral_wave
				+ 0.0733712 * Integral_wave
				+ 4.39355;
			//const Double_t deltaT1 = 
				//+ 5.38622e-9 * pw * pw * pw * pw
				//- 3.46382e-6 * pw * pw * pw
				//+ 0.000718674 * pw * pw
				//- 0.0553169 * pw
				//+ 14.5957;
			const Double_t t_22 = t_2 - deltaT1;
			const Double_t deltaT2 =
				- 1.1097e-5 * pw * pw * pw * pw * pw
				+ 0.00126809 * pw * pw * pw * pw
				- 0.0583122 * pw * pw * pw
				+ 1.3538 * pw * pw
				- 15.9434 * pw
				+ 76.3843;

			const Double_t xValue = Integral_wave;
			const Double_t yValue = t2Truth - t_1;
			scatterX.push_back(xValue);
			scatterY.push_back(yValue);

			hT2_T2ruthDiff->Fill(t_22 - deltaT2 - t2Truth);
			cout << "Trigger" << i << "|ChannelId" << pmtId << "|pw|" << pw << "|t_1|" << deltaT1 << "|t_2|" << t_22 - t2Truth << std::endl;
		}
	}

	// ==================== 新增分箱平均处理部分 ====================
	// 创建Profile直方图进行分箱平均
	const double xmin = *std::min_element(scatterX.begin(), scatterX.end());
	const double xmax = *std::max_element(scatterX.begin(), scatterX.end());
	const int nBins = 50; // 可根据数据量调整分箱数量

	TProfile* hprof = new TProfile("hprof", "Profile of t2Truth - t_1 vs Integral;Integral(ns*mV);t2Truth - t_1 (ns)",
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
	TGraphErrors* avgGraph = new TGraphErrors(avgX.size(), avgX.data(), avgY.data(),
		nullptr, avgYErr.data());
	avgGraph->SetTitle("Binned Average Data");
	avgGraph->SetMarkerStyle(kFullCircle);
	avgGraph->SetMarkerColor(kRed);
	avgGraph->SetMarkerSize(1.2);
	avgGraph->SetLineWidth(2);

	// ==================== 多项式拟合部分 ====================
	// 创建并配置多项式函数（3次到7次）
	const int maxPolyOrder = 7;
	std::vector<TF1*> polyFunctions;
	std::vector<int> colors = { kBlue, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2, kYellow + 2, kRed };

	for (int order = 3; order <= maxPolyOrder; ++order) {
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
	cFit->Divide(2, 2);

	// 子图1: 原始数据分布
	cFit->cd(1);
	TGraph* scatterPlot = new TGraph(scatterX.size(), scatterX.data(), scatterY.data());
	scatterPlot->SetTitle("Data Distribution;Integral(ns*mV);t2Truth - t_1 (ns)");
	scatterPlot->SetMarkerStyle(kFullCircle);
	scatterPlot->SetMarkerColor(kBlue);
	scatterPlot->SetMarkerSize(0.3);
	scatterPlot->Draw("AP");
	avgGraph->Draw("P same");

	// 子图2: 所有多项式拟合比较
	cFit->cd(2);
	avgGraph->GetHistogram()->SetTitle("Polynomial Fitting Comparison");
	avgGraph->Draw("AP");

	TLegend* leg = new TLegend(0.65, 0.6, 0.9, 0.9);
	leg->AddEntry(avgGraph, "Binned Average", "p");

	for (size_t i = 0; i < polyFunctions.size(); ++i) {
		polyFunctions[i]->Draw("same");
		leg->AddEntry(polyFunctions[i],
			TString::Format("Order %d", i + 3), "l");
	}
	leg->Draw();

	// 子图3: 残差分析
	cFit->cd(3);
	TH1F* hResiduals = new TH1F("hResiduals", "Fit Residuals;Residual (ns);Entries",
		100, -5, 5);
	for (size_t i = 0; i < avgX.size(); ++i) {
		double res = avgY[i] - polyFunctions.back()->Eval(avgX[i]);
		hResiduals->Fill(res);
	}
	hResiduals->SetFillColor(kAzure - 3);
	hResiduals->Draw();

	// 子图4: 拟合参数趋势
	cFit->cd(4);
	TGraph* paramTrend = new TGraph();
	for (size_t i = 0; i < polyFunctions.size(); ++i) {
		double chi2 = polyFunctions[i]->GetChisquare();
		paramTrend->SetPoint(i, i + 3, chi2);
	}
	paramTrend->SetTitle("Fit Quality vs Polynomial Order;Order;#chi^{2}");
	paramTrend->SetMarkerStyle(kFullCircle);
	paramTrend->SetMarkerColor(kRed);
	paramTrend->Draw("APL");

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
	for (auto& f : polyFunctions) f->Write();
	hResiduals->Write();
	paramTrend->Write("Chi2Trend");
	scatterPlot->Write("scatterPlot");
	hT2_T2ruthDiff->Write();
	delete scatterPlot;
	f1->Write();
	f1->Close();
	f->Close();
	fb->Close();
}