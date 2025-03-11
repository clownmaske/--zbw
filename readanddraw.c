#include <vector>
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TApplication.h"

using namespace std;

void ReadWaveform(TString inputfile = "JPSim_output.root", 
                 TString outputfile = "JPSim_waveform.root") 
{  
    // 初始化显示参数
    const int MAX_DISPLAY = 5;  // 最多显示5个波形
    int displayCount = 0;        // 已显示波形计数器
    bool showPlots = true;       // 控制是否显示图形
    
    TApplication* app = nullptr;
    if(showPlots) {
        app = new TApplication("app", 0, nullptr);
    }

    TFile* f = TFile::Open(inputfile);
    TTree* t = (TTree*)f->Get("Readout");

    vector<UInt_t>* ChannelId = nullptr;
    vector<UInt_t>* Waveform = nullptr;
    t->SetBranchAddress("ChannelId", &ChannelId);
    t->SetBranchAddress("Waveform", &Waveform);

    t->GetEntry(0);
    Int_t windowSize = Waveform->size()/ChannelId->size();

    // 创建主画布
    TCanvas* mainCanvas = nullptr;
    if(showPlots) {
        mainCanvas = new TCanvas("c_main", "波形显示", 1200, 800);
        mainCanvas->Divide(2,2);  // 2x2网格显示4个波形
    }

    // 遍历触发事件
    bool earlyExit = false;
    for(Int_t i=0; i<t->GetEntries() && !earlyExit; i++) {
        t->GetEntry(i);

        // 遍历PMT通道
        for(size_t j=0; j<ChannelId->size() && !earlyExit; j++) {
            // 只显示前5个波形
            if(displayCount >= MAX_DISPLAY) {
                earlyExit = true;
                break;
            }

            UInt_t PMTId = ChannelId->at(j);
            
            // 创建单个PMT直方图
            TH1F* h1 = new TH1F(
                Form("PMT_%d", PMTId),
                Form("PMT %d 波形;采样点;幅度", PMTId),
                windowSize, 0, windowSize
            );

            // 填充数据
            for(Int_t k=1; k<=windowSize; k++) {
                h1->SetBinContent(k, Waveform->at(k-1 + windowSize*j));
            }

            // 绘制前4个波形到主画布
            if(displayCount < 4) {
                mainCanvas->cd(displayCount+1);
                h1->SetLineColor(displayCount+1);
                h1->Draw("hist");
                displayCount++;
            } 
            // 第5个波形单独弹出窗口
            else if(displayCount == 4) {
                TCanvas* c_final = new TCanvas("c_final", "第5个波形", 600, 400);
                h1->Draw("hist");
                displayCount++;
                c_final->Update();
            }
        }
    }

    // 自动关闭逻辑
    if(showPlots) {
        mainCanvas->Update();
        cout << "已显示5个波形，5秒后自动关闭..." << endl;
        gSystem->Sleep(5000);  // 显示5秒
        delete mainCanvas;
        delete app;
    }

    f->Close();
}