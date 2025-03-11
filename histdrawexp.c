// 保存为 "example.C" 后，在 ROOT 中执行：
// .x example.C

void example() {
    // 1. 初始化随机数生成器（ROOT 默认使用 TRandom3）
    gRandom = new TRandom3(); // 使用当前时间作为种子

    // 2. 创建直方图
    // 参数：名称, 标题, bins数量, 范围下限, 范围上限
    TH1D *hist = new TH1D(
        "myHist", 
        "Example Histogram;Value;Entries",  // 标题格式：标题;X轴标签;Y轴标签
        20,   // 20 个 bin
        0,     // X轴最小值
        100    // X轴最大值
    );

    // 3. 填充 100 个随机数据（均匀分布）
    for (int i = 0; i < 100; ++i) {
        double value = gRandom->Uniform(0, 100); // 生成 [0, 100) 的随机数
        hist->Fill(value);
    }

    // 4. 绘制直方图
    TCanvas *canvas = new TCanvas("canvas", "My Canvas", 800, 600);
    hist->Draw();      // 默认绘制样式为条形图
    canvas->Update();  // 确保图像刷新

    // 5. 保存为图片（可选）
    canvas->SaveAs("histogram.png");
    std::cout << "直方图已保存为 histogram.png" << std::endl;
}