from __future__ import annotations

import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from uuid import uuid4

from dockgui.history import save_run_manifest
from dockgui.models import DockingBackendConfig, DockingBox, DockingJob, DockingRun
from dockgui.parsing import parse_vina_poses, sanitize_filename

CPU_CANDIDATES = ("vina", "vina.exe")
OFFICIAL_GPU_CANDIDATES = ("vina_gpu.exe", "Vina-GPU.exe", "vina_gpu", "Vina-GPU")
OFFICIAL_GPU_KERNEL_CANDIDATES = ("vina_gpu_k.exe", "Vina-GPU-K.exe", "vina_gpu_k", "Vina-GPU-K")
GPU_21_CANDIDATES = (
    "AutoDock-Vina-GPU-2-1.exe",
    "AutoDock-Vina-GPU-2.1.exe",
    "Vina-GPU-2.1.exe",
)
KERNEL_BIN_NAME = "Kernel2_Opt.bin"


def _resolve_executable(explicit_path: str | None, candidates: tuple[str, ...], error_message: str) -> str:
    all_candidates = [explicit_path, *[shutil.which(candidate) for candidate in candidates]]
    for candidate in all_candidates:
        if not candidate:
            continue
        candidate_path = Path(candidate)
        if candidate_path.exists() or shutil.which(str(candidate_path)):
            return str(candidate_path)
    raise FileNotFoundError(error_message)


def resolve_vina_executable(explicit_path: str | None) -> str:
    return _resolve_executable(
        explicit_path,
        CPU_CANDIDATES,
        "未找到 AutoDock Vina 可执行文件，请在侧边栏指定 `vina.exe` 路径。",
    )


def resolve_official_vina_gpu_executable(explicit_path: str | None) -> str:
    return _resolve_executable(
        explicit_path,
        OFFICIAL_GPU_CANDIDATES,
        "未找到官方 Vina-GPU 可执行文件，请在侧边栏指定 `vina_gpu.exe` 或 `Vina-GPU.exe` 路径。",
    )


def resolve_official_vina_gpu_kernel_executable(explicit_path: str | None) -> str:
    return _resolve_executable(
        explicit_path,
        OFFICIAL_GPU_KERNEL_CANDIDATES,
        "未找到官方 Vina-GPU-K 可执行文件，请在侧边栏指定 `vina_gpu_k.exe` 或 `Vina-GPU-K.exe` 路径。",
    )


def resolve_vina_gpu_21_executable(explicit_path: str | None) -> str:
    return _resolve_executable(
        explicit_path,
        GPU_21_CANDIDATES,
        "未找到 Vina-GPU 2.1 可执行文件，请在侧边栏指定本地路径。",
    )


def resolve_backend_executable(backend: DockingBackendConfig) -> str:
    if backend.kind == "vina_cpu":
        return resolve_vina_executable(backend.executable)
    if backend.kind == "vina_gpu_official":
        return resolve_official_vina_gpu_executable(backend.executable)
    if backend.kind == "vina_gpu_21":
        return resolve_vina_gpu_21_executable(backend.executable)
    raise ValueError(f"不支持的 docking 后端：{backend.kind}")


def create_run_dir(base_dir: Path) -> Path:
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    run_dir = base_dir / f"run-{stamp}-{uuid4().hex[:8]}"
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def format_box(box: DockingBox) -> str:
    return (
        f"center_x = {box.center_x:.3f}\n"
        f"center_y = {box.center_y:.3f}\n"
        f"center_z = {box.center_z:.3f}\n"
        f"size_x = {box.size_x:.3f}\n"
        f"size_y = {box.size_y:.3f}\n"
        f"size_z = {box.size_z:.3f}\n"
    )


def write_cpu_config(job: DockingJob, run_dir: Path, receptor_name: str, ligand_name: str) -> Path:
    config_path = run_dir / "vina.conf"
    config_lines = [
        f"receptor = {receptor_name}",
        f"ligand = {ligand_name}",
        "",
        format_box(job.docking_box).strip(),
        "",
        f"exhaustiveness = {job.exhaustiveness}",
        f"num_modes = {job.num_modes}",
        f"energy_range = {job.energy_range:.3f}",
        f"cpu = {job.cpu}",
    ]
    if job.seed is not None:
        config_lines.append(f"seed = {job.seed}")
    config_path.write_text("\n".join(config_lines) + "\n", encoding="utf-8")
    return config_path


def write_official_gpu_config(
    job: DockingJob,
    run_dir: Path,
    receptor_path: Path,
    ligand_path: Path,
    output_path: Path,
    log_path: Path,
    backend: DockingBackendConfig,
) -> Path:
    config_path = run_dir / "vina-gpu.conf"
    config_lines = [
        f"receptor = {receptor_path.resolve()}",
        f"ligand = {ligand_path.resolve()}",
        "",
        format_box(job.docking_box).strip(),
        "",
        f"thread = {backend.thread}",
        f"num_modes = {job.num_modes}",
        f"energy_range = {job.energy_range:.3f}",
        f"out = {output_path.resolve()}",
        f"log = {log_path.resolve()}",
    ]
    if backend.search_depth is not None:
        config_lines.append(f"search_depth = {backend.search_depth}")
    if job.seed is not None:
        config_lines.append(f"seed = {job.seed}")
    config_path.write_text("\n".join(config_lines) + "\n", encoding="utf-8")
    return config_path


def write_gpu_21_config(
    job: DockingJob,
    run_dir: Path,
    receptor_name: str,
    ligand_name: str,
    backend: DockingBackendConfig,
    executable_path: Path,
) -> tuple[Path, Path, Path]:
    ligand_dir = run_dir / "ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)
    output_dir = run_dir / "gpu_out"
    output_dir.mkdir(parents=True, exist_ok=True)

    receptor_path = run_dir / receptor_name
    ligand_path = ligand_dir / ligand_name
    if not receptor_path.exists():
        raise FileNotFoundError(f"缺少受体文件：{receptor_path}")
    if not ligand_path.exists():
        raise FileNotFoundError(f"缺少配体文件：{ligand_path}")

    opencl_path = backend.opencl_binary_path.strip() or str(executable_path.resolve().parent)
    config_path = run_dir / "vina-gpu.conf"
    config_lines = [
        f"receptor = {receptor_path.name}",
        f"ligand_directory = {ligand_dir.name}",
        f"output_directory = {output_dir.name}",
        "",
        format_box(job.docking_box).strip(),
        "",
        f"thread = {backend.thread}",
        f"opencl_binary_path = {opencl_path}",
        f"rilc_bfgs = {1 if backend.rilc_bfgs else 0}",
    ]
    if backend.search_depth is not None:
        config_lines.append(f"search_depth = {backend.search_depth}")
    config_path.write_text("\n".join(config_lines) + "\n", encoding="utf-8")
    return config_path, ligand_dir, output_dir


def _find_gpu_output(output_dir: Path, ligand_name: str) -> Path | None:
    pdbqt_files = sorted(output_dir.glob("*.pdbqt"))
    if not pdbqt_files:
        return None
    ligand_stem = Path(ligand_name).stem.lower()
    preferred = [path for path in pdbqt_files if ligand_stem in path.stem.lower()]
    return preferred[0] if preferred else pdbqt_files[0]


def _build_run(
    job: DockingJob,
    run_dir: Path,
    receptor_path: Path,
    ligand_path: Path,
    config_path: Path,
    output_path: Path,
    log_path: Path,
    stdout: str,
    stderr: str,
    command: list[str],
) -> DockingRun:
    poses = parse_vina_poses(output_path.read_text(encoding="utf-8", errors="ignore"))
    run = DockingRun(
        run_dir=run_dir,
        receptor_path=receptor_path,
        ligand_path=ligand_path,
        config_path=config_path,
        output_path=output_path,
        log_path=log_path,
        stdout=stdout,
        stderr=stderr,
        command=command,
        docking_box=job.docking_box,
        poses=poses,
        created_at=datetime.now().isoformat(timespec="seconds"),
        receptor_label=job.receptor_name,
        ligand_label=job.ligand_name,
        source_label=job.source_label,
        engine_name=job.engine_name,
    )
    save_run_manifest(run)
    return run


def _run_cpu(job: DockingJob, backend: DockingBackendConfig, base_dir: Path) -> DockingRun:
    vina_path = resolve_backend_executable(backend)
    run_dir = create_run_dir(base_dir)

    receptor_name = sanitize_filename(job.receptor_name)
    ligand_name = sanitize_filename(job.ligand_name)
    receptor_path = run_dir / receptor_name
    ligand_path = run_dir / ligand_name
    receptor_path.write_bytes(job.receptor_bytes)
    ligand_path.write_bytes(job.ligand_bytes)

    config_path = write_cpu_config(job, run_dir, receptor_name, ligand_name)
    output_path = run_dir / "poses_out.pdbqt"
    log_path = run_dir / "docking.log"
    command = [
        vina_path,
        "--config",
        str(config_path.name),
        "--out",
        str(output_path.name),
        "--log",
        str(log_path.name),
    ]

    completed = subprocess.run(command, cwd=run_dir, capture_output=True, text=True, check=False)
    stdout = completed.stdout or ""
    stderr = completed.stderr or ""

    if completed.returncode != 0:
        raise RuntimeError((stdout + "\n" + stderr).strip() or "Vina 运行失败。")
    if not output_path.exists():
        raise RuntimeError("Vina 运行完成，但未生成 `poses_out.pdbqt`。")

    return _build_run(job, run_dir, receptor_path, ligand_path, config_path, output_path, log_path, stdout, stderr, command)


def _select_official_gpu_runner(backend: DockingBackendConfig, executable_path: Path) -> tuple[Path, str]:
    kernel_bin = executable_path.resolve().parent / KERNEL_BIN_NAME
    if kernel_bin.exists():
        return executable_path, "vina_gpu"
    kernel_executable = Path(resolve_official_vina_gpu_kernel_executable(backend.kernel_executable or None))
    return kernel_executable, "vina_gpu_k"


def _run_official_gpu(job: DockingJob, backend: DockingBackendConfig, base_dir: Path) -> DockingRun:
    executable_path = Path(resolve_official_vina_gpu_executable(backend.executable))
    runner_path, runner_mode = _select_official_gpu_runner(backend, executable_path)
    run_dir = create_run_dir(base_dir)

    receptor_name = sanitize_filename(job.receptor_name)
    ligand_name = sanitize_filename(job.ligand_name)
    receptor_path = run_dir / receptor_name
    ligand_path = run_dir / ligand_name
    receptor_path.write_bytes(job.receptor_bytes)
    ligand_path.write_bytes(job.ligand_bytes)

    output_path = run_dir / "poses_out.pdbqt"
    log_path = run_dir / "docking.log"
    config_path = write_official_gpu_config(job, run_dir, receptor_path, ligand_path, output_path, log_path, backend)
    command = [
        str(runner_path),
        "--config",
        str(config_path.resolve()),
    ]

    completed = subprocess.run(
        command,
        cwd=executable_path.resolve().parent,
        capture_output=True,
        text=True,
        check=False,
    )
    stdout = completed.stdout or ""
    stderr = completed.stderr or ""

    if completed.returncode != 0:
        hint = ""
        if runner_mode == "vina_gpu":
            hint = "\n提示：若官方 `vina_gpu.exe` 首次运行失败，请确认可执行文件目录下已有 `Kernel2_Opt.bin`。"
        raise RuntimeError(((stdout + "\n" + stderr).strip() + hint).strip() or "官方 Vina-GPU 运行失败。")
    if not output_path.exists():
        raise RuntimeError("官方 Vina-GPU 运行完成，但未生成 `poses_out.pdbqt`。")

    return _build_run(job, run_dir, receptor_path, ligand_path, config_path, output_path, log_path, stdout, stderr, command)


def _run_gpu_21(job: DockingJob, backend: DockingBackendConfig, base_dir: Path) -> DockingRun:
    executable = Path(resolve_vina_gpu_21_executable(backend.executable))
    run_dir = create_run_dir(base_dir)

    receptor_name = sanitize_filename(job.receptor_name)
    ligand_name = sanitize_filename(job.ligand_name)
    receptor_path = run_dir / receptor_name
    receptor_path.write_bytes(job.receptor_bytes)

    ligand_dir = run_dir / "ligands"
    ligand_dir.mkdir(parents=True, exist_ok=True)
    ligand_path = ligand_dir / ligand_name
    ligand_path.write_bytes(job.ligand_bytes)

    config_path, _, output_dir = write_gpu_21_config(job, run_dir, receptor_name, ligand_name, backend, executable)
    log_path = run_dir / "docking.log"
    command = [str(executable), "--config", str(config_path.name)]

    completed = subprocess.run(command, cwd=run_dir, capture_output=True, text=True, check=False)
    stdout = completed.stdout or ""
    stderr = completed.stderr or ""
    log_path.write_text(((stdout + "\n" + stderr).strip() + "\n") if (stdout or stderr) else "", encoding="utf-8")

    if completed.returncode != 0:
        raise RuntimeError((stdout + "\n" + stderr).strip() or "Vina-GPU 2.1 运行失败。")

    output_path = _find_gpu_output(output_dir, ligand_name)
    if output_path is None or not output_path.exists():
        raise RuntimeError("Vina-GPU 2.1 运行完成，但未在输出目录中找到 `PDBQT` 结果。")

    return _build_run(job, run_dir, receptor_path, ligand_path, config_path, output_path, log_path, stdout, stderr, command)


def run_docking(job: DockingJob, backend: DockingBackendConfig, base_dir: Path) -> DockingRun:
    if backend.kind == "vina_cpu":
        return _run_cpu(job, backend, base_dir)
    if backend.kind == "vina_gpu_official":
        return _run_official_gpu(job, backend, base_dir)
    if backend.kind == "vina_gpu_21":
        return _run_gpu_21(job, backend, base_dir)
    raise ValueError(f"不支持的 docking 后端：{backend.kind}")


def run_vina(job: DockingJob, vina_executable: str | None, base_dir: Path) -> DockingRun:
    backend = DockingBackendConfig(kind="vina_cpu", executable=vina_executable or "vina")
    return run_docking(job, backend, base_dir)
