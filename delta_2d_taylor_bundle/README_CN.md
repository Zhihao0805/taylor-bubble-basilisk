# Delta 运行打包说明（2D Taylor bubble）

本目录用于把“前处理 + 计算 + 后处理 + 可视化 + 自动 git”集中到一个地方，方便你后续在 Delta 上 clone 后直接运行。

## 目录内容
- `params_delta.env`
- `01_preprocess_build.sh`
- `02_run_campaigns.sh`
- `03_postprocess_plot.sh`
- `04_check_outputs.sh`
- `05_autogit.sh`
- `run_delta.sbatch`

## 一次性前置条件
1. 你在 Delta 节点上已经能用 `qcc`（Basilisk 环境可用）。
2. 你有 `python3`、`matplotlib`、`ffmpeg`（用于 MP4 检查）。
3. 仓库根目录执行，且保留默认忽略规则（`intermediate/`、`*.png`、`*.mp4` 不提交）。

## 推荐执行顺序（交互式）
```bash
cd ~/projects/Taylor_bubble
bash delta_2d_taylor_bundle/01_preprocess_build.sh
bash delta_2d_taylor_bundle/02_run_campaigns.sh
bash delta_2d_taylor_bundle/03_postprocess_plot.sh
bash delta_2d_taylor_bundle/04_check_outputs.sh
bash delta_2d_taylor_bundle/05_autogit.sh "Delta run: 2D Taylor bubble campaigns"
```

## 推荐执行顺序（批处理）
```bash
cd ~/projects/Taylor_bubble
sbatch delta_2d_taylor_bundle/run_delta.sbatch
```

## 当前默认参数（轻量可跑版）
见 `params_delta.env`，核心是：
- `Re=0.1`
- Grid convergence: `LEV=8,9,10`, `Ca=0.01`
- Tolerance sweep: `TOL=1e-3,1e-5,1e-7`
- Ca sweep: `Ca=0.002,0.005,0.01,0.02,0.05,0.1`
- 代表视频：`Ca=0.01`, 且气泡居中输出 MP4

这些参数是为了保证在算力有限时能稳定跑完完整流程。若你在 Delta 算力更充足，可在 `params_delta.env` 中改为更“接近论文”的更高分辨率/更低 Ca 方案。

## 关键输出文件
- CSV：
  - `intermediate/grid_convergence.csv`
  - `intermediate/tol_sensitivity.csv`
  - `intermediate/ca_sweep.csv`
- 图：
  - `intermediate/grid_convergence.png`
  - `intermediate/tol_sensitivity.png`
  - `intermediate/ca_sweep_vs_theory.png`
- 视频：
  - `intermediate/movie_..._centered.mp4`（气泡居中）

## 自动 git 说明
- `05_autogit.sh` 内部调用仓库已有 `scripts/autopush.sh`。
- 只会提交源码/脚本/markdown；不会提交 `intermediate/` 及大文件产物（受 `.gitignore` 保护）。
