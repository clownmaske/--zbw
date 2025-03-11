void plot_root_file() {
    // 1. 打开文件
    TFile *file = TFile::Open("data.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // 2. 查看内容
    file->ls();

    // 3. 读取直方图
    TH1F *hist = (TH1F*)file->Get("h_energy");
    if (hist) {
        TCanvas *c1 = new TCanvas("c1", "Histogram", 800, 600);
        hist->Draw();
        c1->SaveAs("histogram.png");
    }

    // 4. 读取 TTree
    TTree *tree = (TTree*)file->Get("events");
    if (tree) {
        TCanvas *c2 = new TCanvas("c2", "Tree Plot", 800, 600);
        tree->Draw("energy");  // 直接绘制变量
        c2->SaveAs("energy_distribution.png");
    }

    // 5. 关闭文件
    file->Close();
}