// 读取波形数据
TFile* file = new TFile("waveform.root");
TTree* tree = (TTree*)file->Get("waveformTree");
Readout_t* readout = nullptr;
tree->SetBranchAddress("Readout", &readout);

// 获取第一个事件的波形
tree->GetEntry(0);

// 创建直方图
TH1D* hWaveform = new TH1D("hWaveform", "Waveform;Time (ns);ADC Counts", 
                           readout->Waveform.size(), 0, windowSize);

// 填充数据
for (Int_t i = 0; i < readout->Waveform.size(); ++i) {
    hWaveform->SetBinContent(i + 1, readout->Waveform[i]);
}

// 绘制波形
hWaveform->Draw("HIST");