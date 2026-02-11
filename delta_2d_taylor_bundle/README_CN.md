# Delta 运行打包说明（2D Taylor bubble）

本目录用于把“前处理 + 计算 + 后处理 + 可视化 + 自动 git”集中到一个地方，方便你后续在 Delta 上 clone 后直接运行。
目标是：你只要按步骤执行脚本，就能稳定产出三组 CSV、三张图和居中 MP4。

## 目录内容
- `params_delta.env`
- `taylor_benchmark_2Dpaper_bundle.c`
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

## 关于 C 文件
- 现在打包目录内已包含可编译入口：`delta_2d_taylor_bundle/taylor_benchmark_2Dpaper_bundle.c`。
- `01_preprocess_build.sh` 会直接用该文件编译 `run2d`。

## 从 0 开始（Delta）
```bash
git clone <你的仓库地址> Taylor_bubble
cd Taylor_bubble
git checkout <你要跑的分支>
```

如果你用 module 环境（示例）：
```bash
module purge
# module load gcc python
# source 你的 Basilisk 环境脚本
which qcc
```
确保 `qcc` 能找到后再执行下面流程。

## 推荐执行顺序（交互式）
```bash
cd ~/projects/Taylor_bubble
bash delta_2d_taylor_bundle/01_preprocess_build.sh
bash delta_2d_taylor_bundle/02_run_campaigns.sh
bash delta_2d_taylor_bundle/03_postprocess_plot.sh
bash delta_2d_taylor_bundle/04_check_outputs.sh
bash delta_2d_taylor_bundle/05_autogit.sh "Delta run: 2D Taylor bubble campaigns"
```

每一步作用：
- `01`：备份旧 `intermediate/` 到 `backups/时间戳/`，并编译 `run2d`
- `02`：执行三组 benchmark campaign（grid / tolerance / Ca）
- `03`：绘图
- `04`：检查产物是否齐全 + MP4 元数据 + git 状态
- `05`：自动提交推送（走仓库的 `scripts/autopush.sh`）

## 推荐执行顺序（批处理）
```bash
cd ~/projects/Taylor_bubble
sbatch delta_2d_taylor_bundle/run_delta.sbatch
```

批处理完成后，建议再手动执行一次：
```bash
bash delta_2d_taylor_bundle/04_check_outputs.sh
```
确认无误再执行自动 git：
```bash
bash delta_2d_taylor_bundle/05_autogit.sh "Delta run: 2D Taylor bubble campaigns"
```

## 当前默认参数（轻量可跑版）
见 `params_delta.env`，核心是：
- `Re=0.1`
- Grid convergence: `LEV=8,9,10`, `Ca=0.01`
- Tolerance sweep: `TOL=1e-3,1e-5,1e-7`
- Ca sweep: `Ca=0.002,0.005,0.01,0.02,0.05,0.1`
- 代表视频：`Ca=0.01`, 且气泡居中输出 MP4

这些参数是为了保证在算力有限时能稳定跑完完整流程。若你在 Delta 算力更充足，可在 `params_delta.env` 中改为更“接近论文”的更高分辨率/更低 Ca 方案。

## 如何改参数（非常关键）
统一改 `delta_2d_taylor_bundle/params_delta.env`：
- 网格收敛：`GRID_LEVS`, `GRID_CA`, `GRID_TMAX`
- Poisson 容差：`TOL_VALUES`, `TOL_LEV`, `TOL_TMAX`
- Ca 扫描：`CA_VALUES`, `CA_LEV`, `CA_LEV_LOW`, `CA_LOW_THRESHOLD`, `CA_TMAX`
- 视频：`MOVIE_CA`, `MOVIE_LEV`, `MOVIE_TMAX`

修改后直接重新跑：
```bash
bash delta_2d_taylor_bundle/01_preprocess_build.sh
bash delta_2d_taylor_bundle/02_run_campaigns.sh
bash delta_2d_taylor_bundle/03_postprocess_plot.sh
bash delta_2d_taylor_bundle/04_check_outputs.sh
```

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

## 常见问题
1. `qcc not found`
- 说明 Basilisk 环境没加载。先 `which qcc` 确认。

2. 运行太慢
- 降低 `LEV`、缩短 `*_TMAX`、把 `CA_VALUES` 里极小 Ca 去掉。

3. 想更接近论文
- 提升 `GRID_LEVS` 上限，降低 `CA_VALUES` 最小值（如到 `5e-4`），并增加 walltime。

4. 自动 git 失败
- 检查 SSH/权限；你也可以先 `git status` 手动提交，再推送。
