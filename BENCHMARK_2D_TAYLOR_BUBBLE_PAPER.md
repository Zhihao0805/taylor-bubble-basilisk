# 2D Taylor Bubble Paper Benchmark (Section 4)

## 1) 论文基准提取（来源）
- 参考文献文件：`/home/zhihao/projects/Taylor_bubble/2D_taylor_bubble.pdf`
- 目标章节：Section 4, *2-dimensional Taylor bubble*

按 PDF Section 4 提取到的关键设定：
- 工况：2D 直通道内已形成气泡/液滴，低 Re 与低 Ca。
- 边界条件：入口恒定速度、出口固定压力、壁面 no-slip。
- 无量纲参数：`Re = 0.1` 固定，`Ca` 扫描约 `[0.0005, 0.1]`，并取 `rho_d/rho_c = 1`、`mu_d/mu_c = 1`。
- AMR 指标（论文）：`u_x, u_y, cs, f`，阈值示例 `eps_u=0.007`, `eps_cs=1e-3`, `eps_f=1e-1`。
- 收敛性研究（论文）：
  - 最大网格层级 `LEV=8..12`（对应 `H/Delta` 从约 25 到 410）
  - AMR 阈值敏感性
  - Poisson 容差敏感性（给出 `TOL=1e-5` 后收敛稳定）

论文 Section 4 方程（用于本仓复现）：
- Eq. (7) Bretherton 薄膜关系：
  - `2t/H = 0.643*(3Ca)^(2/3)`
- Eq. (8) Aussillous 扩展：
  - `2t/H = [0.643*(3Ca)^(2/3)] / [1 + 0.643*2.50*(3Ca)^(2/3)]`
- Eq. (9) 速度比关系：
  - `Ud/Uf = 1/(1 - 2t/H)`

## 2) 仓库映射与最小改动路径

### 新增求解器（paper 专用）
- `taylor_benchmark_2Dpaper.c`

实现点：
- 2D 通道（入口速度、出口压力、上下壁 no-slip）
- `Re/Ca` 入参控制
- 计算并输出 `Ud, Uf, Ud/Uf, 2t/H(由Ud/Uf反推), Eq(7)/Eq(8)/Eq(9)理论值, 相对误差`
- 输出网格统计 `LEV, Ncells` 与 Poisson `TOLERANCE`
- 输出时间序列与末段统计（最后 30% 时间窗均值）
- 生成“气泡居中”MP4（每帧按实时 `xcm` 平移重采样 `f` 场）

### 自动化脚本
- `scripts/run_2dpaper_campaigns.py`
  - Campaign 1: 网格收敛（`grid_convergence.csv`）
  - Campaign 2: Poisson 容差敏感性（`tol_sensitivity.csv` + `tol_timeseries.csv`）
  - Campaign 3: Ca 扫描（`ca_sweep.csv`）
  - 代表算例生成居中 MP4
- `scripts/plot_2dpaper_results.py`
  - `grid_convergence.png`
  - `tol_sensitivity.png`
  - `ca_sweep_vs_theory.png`

## 3) 构建/运行命令

仓库已有 markdown 指令中明确命令（原文流程）：
```bash
qcc -O2 -Wall -disable-dimensions taylor_clean_benchmark_axi_steady.c -lm -o run
./run MAXLEVEL Re Ca Lz_over_R TMAX
./run 9 1 0.016 40 2 | tee intermediate/run.log
python3 scripts/run_clean_sweep.py
python3 scripts/plot_fig5.py
```

本次 Section-4 复现实验新增命令：
```bash
qcc -O2 -Wall -disable-dimensions taylor_benchmark_2Dpaper.c -lm -o run2d
python3 scripts/run_2dpaper_campaigns.py \
  --run ./run2d \
  --grid-ca 0.01 --grid-levs 8,9,10 --grid-tmax 0.08 \
  --tol-ca 0.01 --tol-lev 9 --tol-tmax 0.1 \
  --ca-values 0.002,0.005,0.01,0.02,0.05,0.1 \
  --ca-lev 8 --ca-lev-low 9 --ca-low-threshold 0.003 --ca-tmax 0.08 \
  --movie-ca 0.01 --movie-lev 9 --movie-tmax 0.15
python3 scripts/plot_2dpaper_results.py --intermediate intermediate
```

## 4) Campaign 参数与结果

### Campaign 1: Grid convergence
输出：`intermediate/grid_convergence.csv`

| LEV | H/Delta | Ud/Uf(mean,last30%) |
|---:|---:|---:|
| 8 | 25.6 | 0.8311 |
| 9 | 51.2 | 1.0697 |
| 10 | 102.4 | 1.1593 |

图：`intermediate/grid_convergence.png`

### Campaign 2: Poisson tolerance sensitivity
输出：`intermediate/tol_sensitivity.csv`, `intermediate/tol_timeseries.csv`

| TOL | Ud/Uf(mean,last30%) |
|---:|---:|
| 1e-3 | 1.0007 |
| 1e-5 | 0.9915 |
| 1e-7 | 0.9920 |

`Ud/Uf(t)` 波动统计（`tol_timeseries.csv`）：
- `1e-3`: std ≈ 0.1907
- `1e-5`: std ≈ 0.1862
- `1e-7`: std ≈ 0.1867

图：`intermediate/tol_sensitivity.png`

### Campaign 3: Ca sweep
输出：`intermediate/ca_sweep.csv`

| Ca | LEV | Ud/Uf(sim) | Ud/Uf(Auss Eq8-9) | rel_err_auss |
|---:|---:|---:|---:|---:|
| 0.002 | 9 | 0.9465 | 1.0206 | 0.0725 |
| 0.005 | 8 | 0.7724 | 1.0369 | 0.2551 |
| 0.01 | 8 | 0.8311 | 1.0568 | 0.2136 |
| 0.02 | 8 | 0.9080 | 1.0859 | 0.1638 |
| 0.05 | 8 | 0.9203 | 1.1427 | 0.1946 |
| 0.1 | 8 | 0.9278 | 1.2012 | 0.2276 |

误差汇总：
- `max_rel_err_auss = 0.2551`
- `mean_rel_err_auss = 0.1879`

图：`intermediate/ca_sweep_vs_theory.png`

## 5) 居中气泡 MP4
- 文件：`intermediate/movie_lev9_re0p1_ca0p01_tol1p0em05_centered.mp4`
- `ffprobe`:
  - width: 800
  - height: 800
  - duration: 0.32 s
  - frames: 8

说明：
- 在 `movie` 事件中每帧计算气泡质心 `xcm`，将 `f` 场按 `xcm` 平移重采样到 `f_centered`，再输出视频，因此气泡在视频中保持居中。

## 6) 与论文的偏差（明确记录）
1. 网格收敛上限：
- 论文做到 `LEV=12`；本次在当前计算资源下 `LEV=11` 过慢，最终采用 `LEV=8..10`。

2. Ca 范围：
- 论文范围到 `Ca=0.0005`；本次为保证整套流水线可完成，采用 `Ca=0.002..0.1`。

3. AMR 变量：
- 论文包含 `cs` 阈值；本实现采用掩膜壁面（非 embed `cs` 路线），AMR 用 `u_x, u_y, f`，并保留论文同量级阈值 (`eps_u=0.007`, `eps_f=0.1`)。

4. 定标加速：
- 为避免低 Ca 下过小时间步，采用 `Uin=0.1` 并保持无量纲关系 `mu=Uin/Re`, `sigma=mu*Uin/Ca`，保证 `Re/Ca` 定义一致，同时提升运行效率。

5. 稳态判据：
- 论文正文主要给出时间演化与终值对比；本实现采用“最后 30% 时间窗均值”作为稳态统计值。

## 7) 验证清单（PASS/FAIL）
- Build（新增 2D paper 驱动）成功：PASS
- 生成三个 CSV：PASS
  - `intermediate/grid_convergence.csv`
  - `intermediate/tol_sensitivity.csv`
  - `intermediate/ca_sweep.csv`
- 生成三张图：PASS
  - `intermediate/grid_convergence.png`
  - `intermediate/tol_sensitivity.png`
  - `intermediate/ca_sweep_vs_theory.png`
- 生成居中 MP4：PASS
  - `intermediate/movie_lev9_re0p1_ca0p01_tol1p0em05_centered.mp4`
- runtime 大文件未入 git：PASS（`intermediate/`、`*.png`、`*.mp4` 已忽略）

