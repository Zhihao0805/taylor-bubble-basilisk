清除上一次运行的结果，包括可执行文件
rm -rf intermediate
rm -f run
rm -f restart
rm -f log.txt log2.txt *.log
mkdir -p intermediate
rm -f run


dimensionless 编译
qcc -O2 -Wall -disable-dimensions taylor_clean_benchmark_axi_steady.c -lm -o run
//./run MAXLEVEL Re Ca Lz_over_R TMAX 2> log.txt
./run 9 1 0.016 40 2 | tee intermediate/run.log

echo $BASILISK
/home/zhihao/basilisk/src


============================
新增内容：clean benchmark 改动与使用说明（不影响上面原有命令）
============================

一、已做的代码修改（概览）
1) taylor_clean_benchmark_axi_steady.c
   - 新增“论文一致”的薄膜厚度度量：
     * 在 r=0 轴上寻找气泡前后端位置（f<0.5 的轴上区间）
     * 取中点 x_mid = 0.5*(x_front + x_rear)
     * 在 x_mid 截面上找气泡最大半径 r_gas_max
     * h_inf/R = (R - r_gas_max)/R
   - 新增单次运行的汇总输出：
     * 写入 intermediate/summary.csv
     * 字段：Ca,h_over_R,h_taylor_over_R,rel_error
   - 适配 include：
     * 优先使用 adapt_wavelet.h
     * 若不存在则回退 adapt.h

2) 新增脚本（scripts/）
   - scripts/run_clean_sweep.py
     * Python 批量扫 Ca：
       Ca = [0.0015, 0.004, 0.008, 0.016, 0.044, 0.097]
     * 每次运行会生成 intermediate/summary.csv
     * 最终合并为 intermediate/clean_sweep.csv
   - scripts/plot_fig5.py
     * 读取 intermediate/clean_sweep.csv
     * 生成 log-log 图：clean 点 + Taylor 实线 + 4^(2/3)*Taylor 点划线
     * 输出 intermediate/clean_benchmark_fig5.png


二、文件详细使用方法
1) taylor_clean_benchmark_axi_steady.c
   - 编译（与原来一致）：
     qcc -O2 -Wall -disable-dimensions taylor_clean_benchmark_axi_steady.c -lm -o run
   - 单次运行（保持 CLI 不变）：
     ./run MAXLEVEL Re Ca Lz_over_R TMAX | tee intermediate/run.log
   - 运行结束后自动生成：
     intermediate/summary.csv
     * 内容示例：
       Ca,h_over_R,h_taylor_over_R,rel_error
       0.016,0.0xxx,0.0yyy,0.0zzz

2) scripts/run_clean_sweep.py
   - 用途：批量扫 Ca 并生成汇总 CSV
   - 默认运行（使用 ./run 可执行文件）：
     python3 scripts/run_clean_sweep.py
   - 可选参数：
     --run ./run
     --maxlevel 9
     --re 1.0
     --lz-over-r 40
     --tmax 2.0
     --intermediate intermediate
     --csv intermediate/clean_sweep.csv
   - 输出：
     intermediate/clean_sweep.csv

3) scripts/plot_fig5.py
   - 用途：绘制 Fig.5 风格的 log-log 图
   - 运行：
     python3 scripts/plot_fig5.py
   - 可选参数：
     --csv intermediate/clean_sweep.csv
     --out intermediate/clean_benchmark_fig5.png
   - 输出图像：
     intermediate/clean_benchmark_fig5.png


三、验证指标与曲线（供对照）
1) Taylor (Aussillous & Quéré, 2000)
   h_inf/R = 1.34*Ca^(2/3) / (1 + 1.34*2.5*Ca^(2/3))
2) Ratulowski & Chang 上限：
   (h_inf/R)_max = 4^(2/3) * Taylor
3) 误差准则（论文 clean case）：
   最大相对误差 < 8%
