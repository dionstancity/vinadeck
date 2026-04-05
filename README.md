# VinaDock Studio

一个参考 `SwissDock` 工作流的本地可视化 `AutoDock Vina` / 官方 `Vina-GPU` 界面，适合在自己的电脑上完成结构下载、受体清洗、`Meeko` 自动准备、搜索盒设置、单任务或批量 docking，以及结果和历史任务查看。

## 当前功能

- 上传一个或多个 `ligand.pdbqt`，自动切换单任务 / 批量 docking
- 上传 `receptor.pdbqt` 并直接调用本机 `vina.exe`
- 支持切换到本地官方 `Vina-GPU` 后端，直接使用本机 GPU，不依赖云服务器
- 输入 `PDB ID` 从 RCSB 下载结构，按链清洗受体，移除水 / 氢 / 杂原子
- 使用 `Meeko` 直接执行 `PDB/PQR -> receptor.pdbqt`
- 使用 `Meeko` 直接执行 `SDF/MOL2/MOL -> ligand.pdbqt`
- 支持导入包含 `SMILES` 列的 `CSV`，直接批量生成候选配体 `PDBQT`
- 支持导入可选配置文件，自动回填 `center_*`、`size_*`、`num_modes`、`energy_range`、`seed` 等参数
- 选择当前配体或晶体配体自动估算 docking box
- 为每次成功运行保存 `run.json`，支持任务历史回看
- 使用 `3Dmol.js` 查看受体、配体和搜索盒
- 下载 `poses_out.pdbqt`、`docking.log`、`vina.conf` 和准备后的 `PDBQT`
- 提供 `PyInstaller` 打包脚本，生成 Windows 桌面版

## 技术选型

- GUI: `Streamlit`
- Docking engine: 本机 `AutoDock Vina` / 官方 `Vina-GPU`
- 结构下载: `RCSB PDB` 官方接口
- 结构准备: `Meeko`
- 3D 可视化: `3Dmol.js`
- 桌面打包: `PyInstaller`

## 快速开始

### 1. 安装依赖

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### 2. 准备本地 Docking 后端

CPU 模式先确认你能在终端运行 `vina.exe`：

```powershell
vina.exe --help
```

GPU 模式需要准备本地官方 `Vina-GPU` 可执行文件：

- `vina_gpu.exe` 或 `Vina-GPU.exe`
- `vina_gpu_k.exe` 或 `Vina-GPU-K.exe`

应用支持两种本地后端：

- `AutoDock Vina (CPU)`
- `Official Vina-GPU (本地 GPU)`

可以通过环境变量预填路径：

- `VINA_EXE`
- `VINA_GPU_EXE`
- `VINA_GPU_K_EXE`
- `VINA_GPU_THREAD`
- `VINA_GPU_SEARCH_DEPTH`

示例：

```powershell
$env:VINA_EXE="C:\tools\vina_1.2.7_windows_x86_64\bin\vina.exe"
$env:VINA_GPU_EXE="C:\Tools\Vina-GPU\vina_gpu.exe"
$env:VINA_GPU_K_EXE="C:\Tools\Vina-GPU\vina_gpu_k.exe"
streamlit run app.py
```

### 3. 启动 Web 界面

```powershell
streamlit run app.py
```

## 使用流程

### A. 直接 docking

1. 选择“上传 `receptor.pdbqt`”
2. 上传一个或多个 `ligand.pdbqt`
3. 在左侧边栏选择本地 CPU 或 GPU 后端
4. 手工输入搜索盒，或点击“根据当前配体估算 box”
5. CPU 模式可设置 `exhaustiveness`、`num_modes`、`energy_range`
6. GPU 模式使用侧边栏中的 `GPU Thread`、`Search Depth`
7. 点击开始运行

### A1. 本地 GPU 模式说明

- 当前接入的是官方 `Vina-GPU`
- 计算全部在本机执行，不需要部署云服务器
- 程序会自动为每次 GPU 任务生成 `vina-gpu.conf`
- 如果检测不到 `Kernel2_Opt.bin`，程序会自动切到 `vina_gpu_k.exe`
- 一旦在 `Vina-GPU` 目录生成了 `Kernel2_Opt.bin`，后续会优先使用 `vina_gpu.exe`
- 历史任务、批量任务、结果预览和下载功能对 CPU / GPU 后端都可用

### B. 用 `PDB ID helper + Meeko` 准备受体

1. 在受体区域切换到 `PDB ID helper + Meeko`
2. 输入如 `1STP` 的 `PDB ID` 并下载结构
3. 选择需要保留的链、辅因子以及是否去水 / 去氢
4. 如有晶体配体，可点击“根据晶体配体估算 box”
5. 在下方点击“用 Meeko 生成 `receptor.pdbqt`”
6. 如需要，下载清洗后的 `PDB` 或准备后的 `receptor.pdbqt`
7. 继续准备配体并执行 docking

### C. 用 `Meeko` 准备受体和配体

- 受体：
  - 切换到 `上传 PDB/PQR + Meeko`
  - 上传 `PDB` 或 `PQR`
  - 选择电荷模型和 `allow_bad_res`
  - 点击“用 Meeko 生成 `receptor.pdbqt`”
- 配体：
  - 切换到 `上传 SDF/MOL2/MOL + Meeko`
  - 上传一个或多个原始配体文件
  - 选择电荷模型
  - 点击“用 Meeko 准备配体”

说明：

- 当前 `SDF/MOL2/MOL` 模式支持单个文件内包含多个 ligand，会自动拆分并分别生成 `PDBQT`
- `PQR` 输入时可为受体选择 `read` 电荷模型
- 准备成功后，界面会直接把生成的 `PDBQT` 作为 CPU / GPU docking 输入

## 批量 docking

- 直接上传多个 `ligand.pdbqt`，或批量准备多个原始配体
- 每个配体会生成自己的运行目录
- `Batch Summary` 标签页会汇总所有成功和失败任务
- `Current Run` 默认展示本轮中得分最好的任务

## 历史任务

- 每次成功运行都会在 `runs/` 下创建独立目录
- 目录内包含：
  - `vina.conf`
  - `poses_out.pdbqt`
  - `docking.log`
  - 运行时使用的受体 / 配体 `PDBQT`
  - `run.json` 任务摘要
  - GPU 模式下会额外保存 `vina-gpu.conf` 和 GPU 输出目录中的结果文件
- `History` 标签页会读取这些 `run.json` 并允许重新查看 pose

### D. 导入 `CSV(SMILES)` 候选化合物

1. 在配体区域切换到 `导入 CSV(SMILES) + Meeko`
2. 上传一个包含 `SMILES` 列的 `CSV`
3. 选择 `SMILES` 列，以及可选的名称列
4. 选择电荷模型（`gasteiger` 或 `zero`）
5. 点击“用 Meeko 准备 CSV 配体”
6. 程序会把每一行候选化合物转换为 `ligand.pdbqt`，并继续支持批量 docking、预览和下载

说明：

- 默认会优先识别 `smiles`、`smile`、`canonical_smiles`、`isomeric_smiles`
- 如果没有名称列，程序会自动按 `文件名-row-N` 命名
- 空白 `SMILES` 行会自动跳过

### E. 导入可选配置文件

1. 在左侧边栏点击“导入可选配置文件”
2. 选择一个 `vina.conf`、`vina-gpu.conf` 或自定义 `key = value` 文本
3. 点击“应用配置文件”
4. 程序会自动回填搜索盒、CPU/GPU 数值参数，并在识别到 GPU 专属字段时自动切换后端

当前支持的常用键包括：

- `center_x`、`center_y`、`center_z`
- `size_x`、`size_y`、`size_z`
- `exhaustiveness`、`cpu`
- `thread`、`search_depth`
- `num_modes`、`energy_range`、`seed`

## Windows 桌面打包

### 1. 安装桌面打包依赖

```powershell
.\.venv\Scripts\Activate.ps1
pip install -r requirements-desktop.txt
```

### 2. 执行打包

```powershell
powershell -ExecutionPolicy Bypass -File .\scripts\build_desktop.ps1
```

默认会生成：

- `dist\VinaDockStudio\VinaDockStudio.exe`

这是一个 Windows 桌面启动器，内部会直接启动 `Streamlit` 应用。当前为了保证稳定启动，打包结果会保留一个控制台窗口用于输出运行日志。打包后的程序不会自带 `vina.exe` 或官方 `Vina-GPU`，你仍需要在目标机器上放置本地 CPU / GPU docking 程序，并在界面中配置路径。

## 目录结构

- `app.py`：主界面
- `dockgui/models.py`：数据结构
- `dockgui/meeko_tools.py`：Meeko 封装
- `dockgui/parsing.py`：PDBQT 转换、pose 解析、box 估算
- `dockgui/rcsb.py`：PDB ID 下载、受体清洗、晶体配体提取
- `dockgui/vina.py`：CPU / GPU docking 命令行封装
- `dockgui/history.py`：历史任务清单与回读
- `dockgui/viewer.py`：3Dmol 三维显示
- `launcher.py`：桌面版入口
- `VinaDockStudio.spec`：PyInstaller 打包配置
- `scripts/build_desktop.ps1`：Windows 打包脚本

## 后续可扩展

- 优化多分子 `SDF` 的命名、进度展示和失败回报
- 增加口袋预测 / 自动 box 建议
- 增加任务队列和后台执行
- 增加一键下载 `PDB ID` 后直接准备配体
