from __future__ import annotations

import hashlib
import os
import subprocess
from pathlib import Path

import pandas as pd
import streamlit as st

from dockgui.history import list_history, load_run_from_manifest
from dockgui.meeko_tools import prepare_ligands_pdbqt, prepare_receptor_pdbqt, probe_meeko
from dockgui.models import DockingBackendConfig, DockingBox, DockingJob, DockingRun, PdbEntryMetadata, StructureInput
from dockgui.parsing import pdbqt_to_pdb, structure_has_atoms, suggest_box_from_structure
from dockgui.rcsb import clean_pdb_text, extract_residue_pdb, fetch_pdb_entry
from dockgui.vina import resolve_backend_executable, run_docking
from dockgui.viewer import render_structure_viewer

DEFAULT_BOX = {
    "center_x": 0.0,
    "center_y": 0.0,
    "center_z": 0.0,
    "size_x": 20.0,
    "size_y": 20.0,
    "size_z": 20.0,
}

DIRECT_RECEPTOR_MODE = "上传 receptor.pdbqt"
HELPER_RECEPTOR_MODE = "PDB ID helper + Meeko"
RAW_RECEPTOR_MODE = "上传 PDB/PQR + Meeko"
DIRECT_LIGAND_MODE = "上传 ligand.pdbqt"
MEEKO_LIGAND_MODE = "上传 SDF/MOL2/MOL + Meeko"

CPU_BACKEND = "AutoDock Vina (CPU)"
GPU_BACKEND = "Official Vina-GPU (本地 GPU)"


def init_state() -> None:
    for key, value in DEFAULT_BOX.items():
        st.session_state.setdefault(key, value)
    st.session_state.setdefault("last_run", None)
    st.session_state.setdefault("batch_runs", [])
    st.session_state.setdefault("batch_failures", [])
    st.session_state.setdefault("pdb_helper_metadata", None)
    st.session_state.setdefault("pdb_helper_pdb_text", None)
    st.session_state.setdefault("helper_prepared_receptor", None)
    st.session_state.setdefault("helper_prepared_receptor_log", "")
    st.session_state.setdefault("helper_prepared_receptor_signature", "")
    st.session_state.setdefault("raw_prepared_receptor", None)
    st.session_state.setdefault("raw_prepared_receptor_log", "")
    st.session_state.setdefault("raw_prepared_receptor_signature", "")
    st.session_state.setdefault("prepared_ligands", [])
    st.session_state.setdefault("prepared_ligand_logs", {})
    st.session_state.setdefault("prepared_ligand_failures", [])
    st.session_state.setdefault("prepared_ligands_signature", "")


def fingerprint_parts(*parts: object) -> str:
    digest = hashlib.sha1()
    for part in parts:
        if isinstance(part, bytes):
            payload = part
        else:
            payload = str(part).encode("utf-8", errors="ignore")
        digest.update(len(payload).to_bytes(8, "big", signed=False))
        digest.update(payload)
    return digest.hexdigest()


def sync_signature_state(signature_key: str, signature: str, resets: dict[str, object]) -> None:
    if st.session_state.get(signature_key) == signature:
        return
    st.session_state[signature_key] = signature
    for key, value in resets.items():
        st.session_state[key] = value


def upload_to_structure(upload, *, file_format: str, source_label: str) -> StructureInput:
    return StructureInput(
        name=upload.name,
        data=upload.getvalue(),
        file_format=file_format,
        source_label=source_label,
    )


def combined_log_text(stdout: str, stderr: str) -> str:
    parts: list[str] = []
    if stdout.strip():
        parts.append(f"[stdout]\n{stdout.strip()}")
    if stderr.strip():
        parts.append(f"[stderr]\n{stderr.strip()}")
    return "\n\n".join(parts)


def stored_structure(key: str) -> StructureInput | None:
    value = st.session_state.get(key)
    return value if isinstance(value, StructureInput) else None


def stored_structure_list(key: str) -> list[StructureInput]:
    value = st.session_state.get(key, [])
    if not isinstance(value, list):
        return []
    return [item for item in value if isinstance(item, StructureInput)]


def parse_optional_int(value: str) -> int | None:
    stripped = value.strip()
    if not stripped:
        return None
    return int(stripped)


def preferred_path(*paths: str) -> str:
    for candidate in paths:
        if candidate and Path(candidate).exists():
            return candidate
    return paths[0] if paths else ""


def apply_box(box: DockingBox) -> None:
    st.session_state["center_x"] = round(box.center_x, 3)
    st.session_state["center_y"] = round(box.center_y, 3)
    st.session_state["center_z"] = round(box.center_z, 3)
    st.session_state["size_x"] = round(box.size_x, 3)
    st.session_state["size_y"] = round(box.size_y, 3)
    st.session_state["size_z"] = round(box.size_z, 3)


def current_box() -> DockingBox:
    return DockingBox(
        center_x=float(st.session_state["center_x"]),
        center_y=float(st.session_state["center_y"]),
        center_z=float(st.session_state["center_z"]),
        size_x=float(st.session_state["size_x"]),
        size_y=float(st.session_state["size_y"]),
        size_z=float(st.session_state["size_z"]),
    )


def box_valid(box: DockingBox) -> bool:
    return all(value > 0 for value in (box.size_x, box.size_y, box.size_z))


def format_chain_label(chain_id: str) -> str:
    return chain_id or "(blank)"


def show_environment_panel() -> tuple[DockingBackendConfig, Path, bool]:
    meeko_available, meeko_detail = probe_meeko()
    with st.sidebar:
        st.header("运行环境")
        backend_label = st.selectbox(
            "Docking 后端",
            options=[CPU_BACKEND, GPU_BACKEND],
            help="全部在本机执行；选择 GPU 时不会接入云服务器。",
        )

        if backend_label == CPU_BACKEND:
            executable = st.text_input(
                "Vina 可执行文件",
                value=os.environ.get("VINA_EXE", "vina"),
                help="可填写 `vina`、`vina.exe`，或绝对路径。",
            )
            backend = DockingBackendConfig(kind="vina_cpu", executable=executable)
        else:
            default_gpu = preferred_path(
                os.environ.get("VINA_GPU_EXE", ""),
                r"C:\Tools\Vina-GPU\vina_gpu.exe",
                r"C:\Tools\Vina-GPU\Vina-GPU.exe",
                "vina_gpu.exe",
            )
            default_gpu_k = preferred_path(
                os.environ.get("VINA_GPU_K_EXE", ""),
                r"C:\Tools\Vina-GPU\vina_gpu_k.exe",
                r"C:\Tools\Vina-GPU\Vina-GPU-K.exe",
                "vina_gpu_k.exe",
            )
            executable = st.text_input(
                "Vina-GPU 可执行文件",
                value=default_gpu,
                help="填写本地 `vina_gpu.exe` 或 `Vina-GPU.exe` 的完整路径。",
            )
            kernel_executable = st.text_input(
                "Vina-GPU-K 可执行文件",
                value=default_gpu_k,
                help="首次运行或缺少 `Kernel2_Opt.bin` 时，会自动使用 `vina_gpu_k.exe`。",
            )
            thread = int(
                st.number_input(
                    "GPU Thread",
                    min_value=1,
                    value=int(os.environ.get("VINA_GPU_THREAD", "1000")),
                    step=100,
                    help="官方 README 默认值是 `1000`，一般建议小于 `10000`。",
                )
            )
            search_depth_raw = st.text_input(
                "Search Depth",
                value=os.environ.get("VINA_GPU_SEARCH_DEPTH", ""),
                help="留空则交给官方 `Vina-GPU` 自动决定。",
            )
            try:
                search_depth = parse_optional_int(search_depth_raw)
            except ValueError:
                search_depth = None
                st.warning("`Search Depth` 必须是整数，当前已忽略。")
            backend = DockingBackendConfig(
                kind="vina_gpu_official",
                executable=executable,
                kernel_executable=kernel_executable,
                thread=thread,
                search_depth=search_depth,
            )

        runs_dir_raw = st.text_input(
            "结果目录",
            value="runs",
            help="每个单体或批量任务都会创建独立运行目录。",
        )
        runs_dir = Path(runs_dir_raw).expanduser()

        if st.button("检查当前后端", use_container_width=True):
            try:
                resolved = resolve_backend_executable(backend)
                st.success(f"已找到：{resolved}")
                if backend.kind == "vina_gpu_official":
                    kernel_dir = Path(resolved).resolve().parent
                    kernel_bin = kernel_dir / "Kernel2_Opt.bin"
                    st.caption(f"官方 Vina-GPU 目录：{kernel_dir}")
                    if kernel_bin.exists():
                        st.caption(f"已检测到 `{kernel_bin.name}`，将优先运行 `vina_gpu.exe`。")
                    else:
                        st.caption("未检测到 `Kernel2_Opt.bin`，程序会自动切到 `vina_gpu_k.exe`。")
            except Exception as exc:
                st.error(str(exc))

        st.markdown("---")
        st.subheader("Meeko")
        if meeko_available:
            st.success(meeko_detail)
            st.caption("支持 `PDB/PQR -> receptor.pdbqt` 和 `SDF/MOL2/MOL -> ligand.pdbqt`。")
        else:
            st.warning(f"Meeko 不可用：{meeko_detail}")
            st.caption("安装 `requirements.txt` 后，可直接在界面中准备受体和配体。")

        st.markdown("---")
        if backend.kind == "vina_gpu_official":
            st.caption("当前使用本地官方 `Vina-GPU` 执行 docking，不需要云服务器。")
        else:
            st.caption("当前使用本地 `AutoDock Vina` CPU 后端。")
        st.caption("PDB ID helper 会从 RCSB 下载结构并做受体清洗。")
    return backend, runs_dir, meeko_available


def show_pdb_helper_panel(meeko_available: bool) -> dict[str, object] | None:
    input_col, action_col, clear_col = st.columns([1.4, 0.8, 0.6])
    with input_col:
        pdb_id = st.text_input("PDB ID", value="", max_chars=4, placeholder="例如 1STP", key="pdb_helper_id")
    with action_col:
        download_clicked = st.button("下载结构", use_container_width=True, key="download_pdb_entry")
    with clear_col:
        clear_clicked = st.button("清空", use_container_width=True, key="clear_pdb_entry")

    if clear_clicked:
        st.session_state["pdb_helper_metadata"] = None
        st.session_state["pdb_helper_pdb_text"] = None
        st.session_state["helper_prepared_receptor"] = None
        st.session_state["helper_prepared_receptor_log"] = ""
        st.session_state["helper_prepared_receptor_signature"] = ""
        st.rerun()

    if download_clicked:
        try:
            metadata, pdb_text = fetch_pdb_entry(pdb_id)
            st.session_state["pdb_helper_metadata"] = metadata
            st.session_state["pdb_helper_pdb_text"] = pdb_text
            st.success(f"已下载 {metadata.pdb_id}")
        except Exception as exc:
            st.error(str(exc))

    metadata = st.session_state.get("pdb_helper_metadata")
    pdb_text = st.session_state.get("pdb_helper_pdb_text")
    if not isinstance(metadata, PdbEntryMetadata) or not isinstance(pdb_text, str):
        st.caption("输入 PDB ID 后可直接从 RCSB 下载结构，并生成清洗后的受体 PDB。")
        return None

    metric_col1, metric_col2, metric_col3 = st.columns(3)
    with metric_col1:
        st.metric("PDB", metadata.pdb_id)
    with metric_col2:
        st.metric("Chains", len(metadata.chains))
    with metric_col3:
        resolution_text = f"{metadata.resolution:.2f} Å" if metadata.resolution is not None else "n/a"
        st.metric("Resolution", resolution_text)

    st.caption(metadata.title)
    if metadata.experimental_method:
        st.caption(f"Method: {metadata.experimental_method}")

    chain_key = f"pdb_helper_chains_{metadata.pdb_id}"
    keep_key = f"pdb_helper_keep_{metadata.pdb_id}"
    reference_key = f"pdb_helper_reference_{metadata.pdb_id}"
    remove_water_key = f"pdb_helper_remove_water_{metadata.pdb_id}"
    remove_h_key = f"pdb_helper_remove_h_{metadata.pdb_id}"

    selected_chains = st.multiselect(
        "保留链",
        options=metadata.chains,
        default=metadata.chains,
        format_func=format_chain_label,
        key=chain_key,
    )
    ligand_map = {ligand.token: ligand for ligand in metadata.ligands}
    ligand_tokens = [ligand.token for ligand in metadata.ligands]
    kept_residue_tokens = st.multiselect(
        "保留辅因子/配体到清洗后的受体",
        options=ligand_tokens,
        default=[],
        format_func=lambda token: ligand_map[token].label,
        key=keep_key,
    )
    reference_options = [""] + ligand_tokens
    reference_token = st.selectbox(
        "用于自动 box 的晶体配体",
        options=reference_options,
        format_func=lambda token: "不使用" if token == "" else ligand_map[token].label,
        key=reference_key,
    )

    toggle_col1, toggle_col2 = st.columns(2)
    with toggle_col1:
        remove_water = st.checkbox("移除水分子", value=True, key=remove_water_key)
    with toggle_col2:
        remove_hydrogen = st.checkbox("移除氢原子", value=True, key=remove_h_key)

    cleaned_pdb_text, stats = clean_pdb_text(
        pdb_text,
        selected_chains,
        set(kept_residue_tokens),
        remove_water=remove_water,
        remove_hydrogen=remove_hydrogen,
    )

    reference_ligand_pdb = None
    reference_ligand_label = None
    if reference_token:
        try:
            selected_residue = ligand_map[reference_token]
            reference_ligand_pdb = extract_residue_pdb(
                pdb_text,
                selected_residue,
                remove_hydrogen=remove_hydrogen,
            )
            reference_ligand_label = selected_residue.label
        except Exception as exc:
            st.warning(str(exc))

    action_box_col, action_download_col = st.columns(2)
    with action_box_col:
        if st.button(
            "根据晶体配体估算 box",
            use_container_width=True,
            disabled=reference_ligand_pdb is None,
            key=f"estimate_box_from_pdb_{metadata.pdb_id}",
        ):
            try:
                apply_box(suggest_box_from_structure(reference_ligand_pdb or ""))
                st.rerun()
            except Exception as exc:
                st.error(f"无法根据晶体配体估算搜索盒：{exc}")
    with action_download_col:
        st.download_button(
            "下载清洗后的受体 PDB",
            data=cleaned_pdb_text.encode("utf-8"),
            file_name=f"{metadata.pdb_id}_cleaned_receptor.pdb",
            use_container_width=True,
        )

    st.caption(
        "这个 helper 适合下载结构、去水、去杂原子、保留指定辅因子，并用晶体配体自动生成初始搜索盒。"
    )
    if meeko_available:
        st.info("清洗完成后，可直接在下方用 Meeko 生成 `receptor.pdbqt`。")
    else:
        st.info("当前未检测到 Meeko。你仍可下载清洗后的 `PDB`，稍后手工转换为 `receptor.pdbqt`。")

    with st.expander("查看受体清洗摘要"):
        summary = pd.DataFrame(
            [
                {"Item": "Kept protein atoms", "Count": stats["kept_protein_atoms"]},
                {"Item": "Kept cofactor atoms", "Count": stats["kept_cofactor_atoms"]},
                {"Item": "Removed chain atoms", "Count": stats["removed_chain_atoms"]},
                {"Item": "Removed water atoms", "Count": stats["removed_water_atoms"]},
                {"Item": "Removed hetero atoms", "Count": stats["removed_hetero_atoms"]},
                {"Item": "Removed hydrogen atoms", "Count": stats["removed_hydrogen_atoms"]},
                {"Item": "Removed altloc atoms", "Count": stats["removed_altloc_atoms"]},
            ]
        )
        st.dataframe(summary, use_container_width=True, hide_index=True)

    if metadata.ligands:
        with st.expander("查看下载结构中的异源残基"):
            ligand_rows = [
                {
                    "Token": ligand.token,
                    "Residue": ligand.res_name,
                    "Chain": format_chain_label(ligand.chain_id),
                    "Residue ID": f"{ligand.res_seq}{ligand.insertion_code}",
                    "Atoms": ligand.atom_count,
                }
                for ligand in metadata.ligands
            ]
            st.dataframe(pd.DataFrame(ligand_rows), use_container_width=True, hide_index=True)

    return {
        "metadata": metadata,
        "cleaned_pdb_text": cleaned_pdb_text,
        "reference_ligand_pdb": reference_ligand_pdb,
        "reference_ligand_label": reference_ligand_label,
    }


def show_helper_receptor_prepare_panel(helper_data: dict[str, object] | None) -> StructureInput | None:
    if not helper_data:
        return None

    metadata = helper_data.get("metadata")
    cleaned_pdb_text = helper_data.get("cleaned_pdb_text")
    if not isinstance(metadata, PdbEntryMetadata) or not isinstance(cleaned_pdb_text, str):
        return None

    st.markdown("#### Meeko 受体准备")
    st.caption("将当前清洗后的受体 `PDB` 直接转换为 `receptor.pdbqt`，并把结果作为本次 docking 输入。")

    charge_model = st.selectbox(
        "受体电荷模型",
        options=["gasteiger", "zero"],
        key="helper_receptor_charge_model",
    )
    allow_bad_res = st.checkbox(
        "删除缺失原子的残基而不是报错",
        value=True,
        key="helper_receptor_allow_bad_res",
    )

    signature = fingerprint_parts(metadata.pdb_id, cleaned_pdb_text.encode("utf-8"), charge_model, allow_bad_res)
    sync_signature_state(
        "helper_prepared_receptor_signature",
        signature,
        {
            "helper_prepared_receptor": None,
            "helper_prepared_receptor_log": "",
        },
    )

    if st.button("用 Meeko 生成 receptor.pdbqt", use_container_width=True, key="prepare_helper_receptor"):
        try:
            prepared = prepare_receptor_pdbqt(
                f"{metadata.pdb_id}_cleaned_receptor.pdb",
                cleaned_pdb_text.encode("utf-8"),
                input_format="pdb",
                allow_bad_res=allow_bad_res,
                charge_model=charge_model,
            )
            st.session_state["helper_prepared_receptor"] = StructureInput(
                name=prepared.name,
                data=prepared.data,
                file_format="pdbqt",
                source_label=f"Meeko receptor from PDB helper {metadata.pdb_id}",
            )
            st.session_state["helper_prepared_receptor_log"] = combined_log_text(prepared.stdout, prepared.stderr)
            st.success(f"已生成 {prepared.name}")
        except Exception as exc:
            st.error(f"Meeko 受体准备失败：{exc}")

    prepared_receptor = stored_structure("helper_prepared_receptor")
    if prepared_receptor is None:
        return None

    st.success(f"当前用于 docking 的受体：`{prepared_receptor.name}`")
    st.download_button(
        "下载 receptor.pdbqt",
        data=prepared_receptor.data,
        file_name=prepared_receptor.name,
        use_container_width=True,
        key="download_helper_prepared_receptor",
    )
    log_text = str(st.session_state.get("helper_prepared_receptor_log", ""))
    if log_text:
        with st.expander("查看 Meeko 受体准备日志"):
            st.text(log_text)
    return prepared_receptor


def show_raw_receptor_prepare_panel(meeko_available: bool) -> StructureInput | None:
    if not meeko_available:
        st.warning("当前未检测到 Meeko，无法直接准备受体。")
        return None

    raw_upload = st.file_uploader(
        "上传原始受体 `PDB` 或 `PQR`",
        type=["pdb", "pqr"],
        key="raw_receptor_upload",
    )
    if raw_upload is None:
        st.caption("支持直接上传 `PDB` 或 `PQR`，并在界面中转换为 `receptor.pdbqt`。")
        return None

    input_format = Path(raw_upload.name).suffix.lower().lstrip(".")
    charge_options = ["gasteiger", "zero"] + (["read"] if input_format == "pqr" else [])
    charge_model = st.selectbox(
        "受体电荷模型",
        options=charge_options,
        key="raw_receptor_charge_model",
        help="`read` 仅在输入为 `PQR` 时可用。",
    )
    allow_bad_res = st.checkbox(
        "删除缺失原子的残基而不是报错",
        value=True,
        key="raw_receptor_allow_bad_res",
    )

    signature = fingerprint_parts(raw_upload.name, raw_upload.getvalue(), charge_model, allow_bad_res)
    sync_signature_state(
        "raw_prepared_receptor_signature",
        signature,
        {
            "raw_prepared_receptor": None,
            "raw_prepared_receptor_log": "",
        },
    )

    if st.button("用 Meeko 生成 receptor.pdbqt", use_container_width=True, key="prepare_raw_receptor"):
        try:
            prepared = prepare_receptor_pdbqt(
                raw_upload.name,
                raw_upload.getvalue(),
                input_format=input_format,
                allow_bad_res=allow_bad_res,
                charge_model=charge_model,
            )
            st.session_state["raw_prepared_receptor"] = StructureInput(
                name=prepared.name,
                data=prepared.data,
                file_format="pdbqt",
                source_label=f"Meeko receptor from uploaded {input_format.upper()}",
            )
            st.session_state["raw_prepared_receptor_log"] = combined_log_text(prepared.stdout, prepared.stderr)
            st.success(f"已生成 {prepared.name}")
        except Exception as exc:
            st.error(f"Meeko 受体准备失败：{exc}")

    prepared_receptor = stored_structure("raw_prepared_receptor")
    if prepared_receptor is None:
        return None

    st.success(f"当前用于 docking 的受体：`{prepared_receptor.name}`")
    st.download_button(
        "下载 receptor.pdbqt",
        data=prepared_receptor.data,
        file_name=prepared_receptor.name,
        use_container_width=True,
        key="download_raw_prepared_receptor",
    )
    log_text = str(st.session_state.get("raw_prepared_receptor_log", ""))
    if log_text:
        with st.expander("查看 Meeko 受体准备日志"):
            st.text(log_text)
    return prepared_receptor


def show_receptor_panel(meeko_available: bool) -> tuple[StructureInput | None, dict[str, object] | None]:
    st.subheader("1. Receptor")
    mode = st.radio(
        "受体输入方式",
        options=[DIRECT_RECEPTOR_MODE, HELPER_RECEPTOR_MODE, RAW_RECEPTOR_MODE],
        horizontal=True,
        key="receptor_input_mode",
    )

    if mode == DIRECT_RECEPTOR_MODE:
        receptor_upload = st.file_uploader(
            "上传受体 `receptor.pdbqt`",
            type=["pdbqt"],
            key="receptor_upload",
        )
        receptor_input = None
        if receptor_upload is not None:
            receptor_input = upload_to_structure(
                receptor_upload,
                file_format="pdbqt",
                source_label="Uploaded receptor PDBQT",
            )
        st.caption("执行 Vina 时直接使用这个文件。")
        return receptor_input, None

    if mode == HELPER_RECEPTOR_MODE:
        helper_data = show_pdb_helper_panel(meeko_available)
        receptor_input = show_helper_receptor_prepare_panel(helper_data) if meeko_available else None
        return receptor_input, helper_data

    return show_raw_receptor_prepare_panel(meeko_available), None


def show_direct_ligand_panel() -> tuple[list[StructureInput], StructureInput | None]:
    ligand_uploads = st.file_uploader(
        "上传一个或多个配体 `PDBQT`",
        type=["pdbqt"],
        accept_multiple_files=True,
        key="ligand_uploads",
    )
    if not ligand_uploads:
        st.caption("上传一个配体会执行单任务 docking；上传多个配体会执行批量 docking。")
        return [], None

    ligand_inputs = [
        upload_to_structure(upload, file_format="pdbqt", source_label="Uploaded ligand PDBQT")
        for upload in ligand_uploads
    ]
    if len(ligand_inputs) > 1:
        st.info(f"已检测到 {len(ligand_inputs)} 个配体，将按顺序执行批量 docking。")
    else:
        st.caption("当前是单任务 docking。")

    preview_index = st.selectbox(
        "预览配体",
        options=list(range(len(ligand_inputs))),
        format_func=lambda index: ligand_inputs[index].name,
        key="ligand_preview_index",
    )
    return ligand_inputs, ligand_inputs[preview_index]


def show_meeko_ligand_panel(meeko_available: bool) -> tuple[list[StructureInput], StructureInput | None]:
    if not meeko_available:
        st.warning("当前未检测到 Meeko，无法直接准备配体。")
        return [], None

    raw_uploads = st.file_uploader(
        "上传一个或多个原始配体 `SDF/MOL2/MOL`",
        type=["sdf", "mol2", "mol"],
        accept_multiple_files=True,
        key="raw_ligand_uploads",
    )
    if not raw_uploads:
        st.caption("该模式会先把原始配体转换为 `ligand.pdbqt`，再执行 docking。")
        return [], None

    charge_model = st.selectbox(
        "配体电荷模型",
        options=["gasteiger", "zero", "read"],
        key="meeko_ligand_charge_model",
        help="如果 `SDF/MOL2` 已经带有部分电荷，可尝试 `read`。",
    )
    if charge_model == "read":
        st.info("`read` 只适用于已经带部分电荷的 `SDF/MOL2`。普通化合物库通常应使用 `gasteiger`。")
    signature_parts: list[object] = [charge_model]
    for upload in raw_uploads:
        signature_parts.extend([upload.name, upload.getvalue()])
    signature = fingerprint_parts(*signature_parts)
    sync_signature_state(
        "prepared_ligands_signature",
        signature,
        {
            "prepared_ligands": [],
            "prepared_ligand_logs": {},
            "prepared_ligand_failures": [],
        },
    )

    if st.button("用 Meeko 准备配体", use_container_width=True, key="prepare_ligands"):
        successes: list[StructureInput] = []
        failures: list[dict[str, str]] = []
        logs: dict[str, str] = {}
        progress = st.progress(0.0)
        status = st.empty()

        for index, upload in enumerate(raw_uploads):
            status.write(f"正在准备 {upload.name} ({index + 1}/{len(raw_uploads)})")
            try:
                prepared_files = prepare_ligands_pdbqt(
                    upload.name,
                    upload.getvalue(),
                    charge_model=charge_model,
                )
                for prepared_file in prepared_files:
                    successes.append(
                        StructureInput(
                            name=prepared_file.name,
                            data=prepared_file.data,
                            file_format="pdbqt",
                            source_label=f"Meeko ligand from {upload.name}",
                        )
                    )
                    logs[prepared_file.name] = combined_log_text(prepared_file.stdout, prepared_file.stderr)
            except Exception as exc:
                failures.append({"Ligand": upload.name, "Error": str(exc)})
            progress.progress((index + 1) / len(raw_uploads))

        status.empty()
        progress.empty()
        st.session_state["prepared_ligands"] = successes
        st.session_state["prepared_ligand_logs"] = logs
        st.session_state["prepared_ligand_failures"] = failures

        if successes:
            st.success(f"已准备 {len(successes)} 个配体。")
        if failures:
            st.warning(f"有 {len(failures)} 个配体准备失败。")

    prepared_ligands = stored_structure_list("prepared_ligands")
    ligand_failures = st.session_state.get("prepared_ligand_failures", [])
    ligand_logs = st.session_state.get("prepared_ligand_logs", {})
    if not isinstance(ligand_logs, dict):
        ligand_logs = {}

    if not prepared_ligands:
        return [], None

    if len(prepared_ligands) > 1:
        st.info(f"当前有 {len(prepared_ligands)} 个已准备好的配体，将按顺序执行批量 docking。")
    else:
        st.caption("当前是单任务 docking。")

    st.dataframe(
        pd.DataFrame(
            [
                {
                    "Ligand": ligand.name,
                    "Source": ligand.source_label,
                    "Size (bytes)": len(ligand.data),
                }
                for ligand in prepared_ligands
            ]
        ),
        use_container_width=True,
        hide_index=True,
    )

    preview_index = st.selectbox(
        "预览已准备好的配体",
        options=list(range(len(prepared_ligands))),
        format_func=lambda index: prepared_ligands[index].name,
        key="prepared_ligand_preview_index",
    )
    selected_ligand = prepared_ligands[preview_index]

    with st.expander("下载已准备好的配体"):
        for ligand in prepared_ligands:
            st.download_button(
                f"下载 {ligand.name}",
                data=ligand.data,
                file_name=ligand.name,
                use_container_width=True,
                key=f"download_prepared_ligand_{ligand.name}",
            )

    log_text = str(ligand_logs.get(selected_ligand.name, ""))
    if log_text:
        with st.expander(f"查看 {selected_ligand.name} 的 Meeko 日志"):
            st.text(log_text)

    if ligand_failures:
        with st.expander("查看配体准备失败记录"):
            st.dataframe(pd.DataFrame(ligand_failures), use_container_width=True, hide_index=True)

    return prepared_ligands, selected_ligand


def show_ligand_panel(meeko_available: bool) -> tuple[list[StructureInput], StructureInput | None]:
    st.subheader("2. Ligands")
    mode = st.radio(
        "配体输入方式",
        options=[DIRECT_LIGAND_MODE, MEEKO_LIGAND_MODE],
        horizontal=True,
        key="ligand_input_mode",
    )
    if mode == DIRECT_LIGAND_MODE:
        return show_direct_ligand_panel()
    return show_meeko_ligand_panel(meeko_available)


def show_box_panel(uploaded_ligand_text: str | None, reference_ligand_pdb: str | None) -> DockingBox:
    st.subheader("3. Binding Site")
    action_col1, action_col2, hint_col = st.columns([1, 1, 1.6])
    with action_col1:
        if st.button(
            "根据上传配体估算 box",
            use_container_width=True,
            disabled=uploaded_ligand_text is None,
            key="estimate_box_from_upload",
        ):
            try:
                apply_box(suggest_box_from_structure(uploaded_ligand_text or ""))
                st.rerun()
            except Exception as exc:
                st.error(f"无法根据上传配体估算搜索盒：{exc}")
    with action_col2:
        if st.button(
            "根据晶体配体估算 box",
            use_container_width=True,
            disabled=reference_ligand_pdb is None,
            key="estimate_box_from_reference",
        ):
            try:
                apply_box(suggest_box_from_structure(reference_ligand_pdb or ""))
                st.rerun()
            except Exception as exc:
                st.error(f"无法根据晶体配体估算搜索盒：{exc}")
    with hint_col:
        st.caption("类似 SwissDock 的 workflow：支持手工输入，或根据上传配体/晶体配体自动生成初始 box。")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.number_input("Center X", key="center_x", step=0.5, format="%.3f")
        st.number_input("Size X", key="size_x", min_value=1.0, step=1.0, format="%.3f")
    with col2:
        st.number_input("Center Y", key="center_y", step=0.5, format="%.3f")
        st.number_input("Size Y", key="size_y", min_value=1.0, step=1.0, format="%.3f")
    with col3:
        st.number_input("Center Z", key="center_z", step=0.5, format="%.3f")
        st.number_input("Size Z", key="size_z", min_value=1.0, step=1.0, format="%.3f")
    return current_box()


def show_options_panel(backend: DockingBackendConfig) -> tuple[int, int, float, int, int | None]:
    st.subheader("4. Docking Options")
    if backend.kind != "vina_cpu":
        st.info("当前使用 GPU 后端。该模式主要读取侧边栏中的 `GPU Thread`、`Search Depth`，以及这里的输出控制参数。")
        col1, col2 = st.columns(2)
        with col1:
            num_modes = st.number_input("Num Modes", min_value=1, value=9, step=1, key="gpu_num_modes")
        with col2:
            energy_range = st.number_input(
                "Energy Range",
                min_value=0.5,
                value=3.0,
                step=0.5,
                format="%.1f",
                key="gpu_energy_range",
            )
        seed_raw = st.text_input("Random Seed", value="", help="留空则由官方 `Vina-GPU` 自动生成。", key="gpu_seed")
        seed = None
        if seed_raw.strip():
            try:
                seed = int(seed_raw)
            except ValueError:
                st.warning("Random Seed 必须是整数，当前将忽略该值。")
        return 8, int(num_modes), float(energy_range), 0, seed

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        exhaustiveness = st.number_input("Exhaustiveness", min_value=1, value=8, step=1)
    with col2:
        num_modes = st.number_input("Num Modes", min_value=1, value=9, step=1)
    with col3:
        energy_range = st.number_input("Energy Range", min_value=0.5, value=3.0, step=0.5, format="%.1f")
    with col4:
        cpu = st.number_input("CPU", min_value=0, value=0, step=1, help="0 表示交给 Vina 自动决定。")

    seed_raw = st.text_input("Random Seed", value="", help="留空则由 Vina 自动生成。")
    seed = None
    if seed_raw.strip():
        try:
            seed = int(seed_raw)
        except ValueError:
            st.warning("Random Seed 必须是整数，当前将忽略该值。")
    return int(exhaustiveness), int(num_modes), float(energy_range), int(cpu), seed


def build_source_label(receptor_input: StructureInput, ligand_input: StructureInput) -> str:
    receptor_source = receptor_input.source_label or receptor_input.name
    ligand_source = ligand_input.source_label or ligand_input.name
    return f"{receptor_source}; ligand: {ligand_source}"


def preview_panel(
    receptor_input: StructureInput | None,
    helper_data: dict[str, object] | None,
    ligand_preview: StructureInput | None,
    box: DockingBox,
) -> None:
    st.subheader("5. Structure Preview")

    receptor_pdb = None
    receptor_label = "无"
    if helper_data and isinstance(helper_data.get("cleaned_pdb_text"), str):
        receptor_pdb = str(helper_data["cleaned_pdb_text"])
        metadata = helper_data.get("metadata")
        receptor_label = f"PDB helper: {metadata.pdb_id}" if isinstance(metadata, PdbEntryMetadata) else "PDB helper"
    elif receptor_input is not None and structure_has_atoms(receptor_input.text):
        receptor_pdb = pdbqt_to_pdb(receptor_input.text)
        receptor_label = receptor_input.name

    ligand_pdb = None
    ligand_label = "无"
    if ligand_preview is not None and structure_has_atoms(ligand_preview.text):
        ligand_pdb = pdbqt_to_pdb(ligand_preview.text)
        ligand_label = ligand_preview.name
    elif helper_data and isinstance(helper_data.get("reference_ligand_pdb"), str):
        ligand_pdb = str(helper_data["reference_ligand_pdb"])
        ligand_label = str(helper_data.get("reference_ligand_label") or "Crystal ligand")

    if receptor_pdb is None and ligand_pdb is None:
        st.info("上传结构、下载 PDB 或完成 Meeko 准备后，会在这里显示三维预览。")
        return

    st.caption(f"Receptor: {receptor_label} | Ligand: {ligand_label}")
    render_structure_viewer(receptor_pdb, ligand_pdb, box)


def run_panel(
    receptor_input: StructureInput | None,
    ligand_inputs: list[StructureInput],
    box: DockingBox,
    backend: DockingBackendConfig,
    runs_dir: Path,
    options: tuple[int, int, float, int, int | None],
    helper_data: dict[str, object] | None,
) -> None:
    st.subheader("6. Execute")
    exhaustiveness, num_modes, energy_range, cpu, seed = options
    ligand_count = len(ligand_inputs)
    run_label = "开始 Docking" if ligand_count <= 1 else f"开始批量 Docking（{ligand_count} 个配体）"
    disabled = receptor_input is None or ligand_count == 0 or not box_valid(box)

    receptor_mode = st.session_state.get("receptor_input_mode", DIRECT_RECEPTOR_MODE)
    ligand_mode = st.session_state.get("ligand_input_mode", DIRECT_LIGAND_MODE)
    if helper_data and receptor_input is None:
        st.info("PDB ID helper 已经帮你清洗受体并估算 box。继续用 Meeko 生成 `receptor.pdbqt` 后即可运行 Vina。")
    elif receptor_mode != DIRECT_RECEPTOR_MODE and receptor_input is None:
        st.info("先完成 Meeko 受体准备，再执行 Vina。")

    if ligand_mode != DIRECT_LIGAND_MODE and ligand_count == 0:
        st.caption("先用 Meeko 准备配体，再开始 docking。")
    elif receptor_input is not None and ligand_count == 0:
        st.caption("上传或准备配体后即可开始 docking。")

    if st.button(run_label, type="primary", use_container_width=True, disabled=disabled):
        runs_dir.mkdir(parents=True, exist_ok=True)
        successes: list[DockingRun] = []
        failures: list[dict[str, str]] = []
        progress = st.progress(0.0)
        status = st.empty()

        for index, ligand_input in enumerate(ligand_inputs):
            status.write(f"正在运行 {ligand_input.name} ({index + 1}/{ligand_count})")
            job = DockingJob(
                receptor_name=receptor_input.name,
                receptor_bytes=receptor_input.data,
                ligand_name=ligand_input.name,
                ligand_bytes=ligand_input.data,
                docking_box=box,
                exhaustiveness=exhaustiveness,
                num_modes=num_modes,
                energy_range=energy_range,
                cpu=cpu,
                seed=seed,
                source_label=build_source_label(receptor_input, ligand_input),
                engine_name=backend.label,
            )
            try:
                successes.append(run_docking(job, backend=backend, base_dir=runs_dir))
            except Exception as exc:
                failures.append({"Ligand": ligand_input.name, "Error": str(exc)})
            progress.progress((index + 1) / ligand_count)

        status.empty()
        progress.empty()
        st.session_state["batch_runs"] = successes
        st.session_state["batch_failures"] = failures

        if successes:
            st.session_state["last_run"] = min(
                successes,
                key=lambda run: run.best_affinity if run.best_affinity is not None else float("inf"),
            )
            if len(successes) == 1 and not failures:
                st.success(f"Docking 完成，结果保存在 `{successes[0].run_dir}`")
            else:
                st.success(f"批量 docking 完成：成功 {len(successes)} 个，失败 {len(failures)} 个。")
        else:
            st.session_state["last_run"] = None
            st.error("本次批量任务全部失败，请检查当前本地 docking 后端路径和输入文件。")

        if failures:
            st.warning("有部分配体运行失败，可在 Batch Summary 中查看错误信息。")


def render_run_details(run: DockingRun | None, *, context_key: str) -> None:
    if run is None:
        st.caption("这里会显示最近一次成功运行的 docking 结果。")
        return

    if not run.poses:
        st.warning("该任务没有解析到 pose，请检查 `docking.log`。")
        return

    metric_col1, metric_col2, metric_col3 = st.columns(3)
    with metric_col1:
        best_affinity = f"{run.best_affinity:.2f} kcal/mol" if run.best_affinity is not None else "n/a"
        st.metric("最佳打分", best_affinity)
    with metric_col2:
        st.metric("Pose 数量", len(run.poses))
    with metric_col3:
        st.metric("运行时间", run.created_at or "n/a")

    st.caption(f"Engine: {run.engine_name} | Source: {run.source_label or 'n/a'} | Run dir: {run.run_dir}")
    pose_table = pd.DataFrame([pose.as_row() for pose in run.poses])
    st.dataframe(pose_table, use_container_width=True, hide_index=True)

    pose_index = st.selectbox(
        "查看 pose",
        options=list(range(len(run.poses))),
        format_func=lambda index: f"Mode {run.poses[index].mode} | {run.poses[index].affinity:.2f} kcal/mol",
        key=f"pose_selector_{context_key}_{run.run_dir.name}",
    )
    selected_pose = run.poses[pose_index]
    receptor_text = run.receptor_path.read_text(encoding="utf-8", errors="ignore") if run.receptor_path.exists() else ""
    receptor_pdb = pdbqt_to_pdb(receptor_text) if structure_has_atoms(receptor_text) else None
    ligand_pdb = pdbqt_to_pdb(selected_pose.pdbqt_text)
    render_structure_viewer(receptor_pdb, ligand_pdb, run.docking_box)

    download_col1, download_col2, download_col3 = st.columns(3)
    with download_col1:
        if run.output_path.exists():
            st.download_button(
                "下载 poses_out.pdbqt",
                data=run.output_path.read_bytes(),
                file_name=run.output_path.name,
                use_container_width=True,
                key=f"download_pose_{context_key}_{run.run_dir.name}",
            )
    with download_col2:
        if run.log_path.exists():
            st.download_button(
                "下载 docking.log",
                data=run.log_path.read_bytes(),
                file_name=run.log_path.name,
                use_container_width=True,
                key=f"download_log_{context_key}_{run.run_dir.name}",
            )
    with download_col3:
        if run.config_path.exists():
            st.download_button(
                "下载 vina.conf",
                data=run.config_path.read_bytes(),
                file_name=run.config_path.name,
                use_container_width=True,
                key=f"download_config_{context_key}_{run.run_dir.name}",
            )

    with st.expander("查看命令与日志"):
        if run.command:
            st.code(subprocess.list2cmdline(run.command), language="bash")
        st.text(run.stdout or "(stdout 为空)")
        if run.stderr:
            st.text(run.stderr)


def render_batch_summary() -> None:
    batch_runs: list[DockingRun] = st.session_state.get("batch_runs", [])
    batch_failures: list[dict[str, str]] = st.session_state.get("batch_failures", [])

    if not batch_runs and not batch_failures:
        st.caption("批量 docking 后，这里会显示所有配体的汇总结果。")
        return

    if batch_runs:
        summary_rows = [run.as_batch_row() for run in batch_runs]
        st.dataframe(pd.DataFrame(summary_rows), use_container_width=True, hide_index=True)
        inspect_index = st.selectbox(
            "查看成功任务",
            options=list(range(len(batch_runs))),
            format_func=lambda index: batch_runs[index].ligand_label or batch_runs[index].ligand_path.name,
            key="batch_inspect_index",
        )
        render_run_details(batch_runs[inspect_index], context_key="batch")

    if batch_failures:
        with st.expander("查看失败任务"):
            st.dataframe(pd.DataFrame(batch_failures), use_container_width=True, hide_index=True)


def render_history_panel(runs_dir: Path) -> None:
    records = list_history(runs_dir)
    if not records:
        st.caption("`runs/` 下还没有可加载的历史任务。")
        return

    st.dataframe(
        pd.DataFrame([record.as_row() for record in records]),
        use_container_width=True,
        hide_index=True,
    )
    selected_index = st.selectbox(
        "加载历史任务",
        options=list(range(len(records))),
        format_func=lambda index: records[index].label,
        key="history_run_index",
    )
    try:
        selected_run = load_run_from_manifest(records[selected_index].manifest_path)
    except Exception as exc:
        st.error(f"无法读取历史任务：{exc}")
        return

    render_run_details(selected_run, context_key="history")


def main() -> None:
    st.set_page_config(page_title="VinaDock Studio", layout="wide")
    init_state()

    st.title("VinaDock Studio")
    st.caption(
        "参考 SwissDock workflow 的本地可视化 AutoDock Vina / Vina-GPU 工具：支持 PDB ID 下载与受体清洗、Meeko 自动准备、本地 CPU/GPU docking、任务历史和结果三维查看。"
    )

    backend, runs_dir, meeko_available = show_environment_panel()
    left, right = st.columns([1.15, 0.85])

    with left:
        receptor_input, helper_data = show_receptor_panel(meeko_available)
        ligand_inputs, ligand_preview = show_ligand_panel(meeko_available)
        reference_ligand_pdb = None
        if helper_data and isinstance(helper_data.get("reference_ligand_pdb"), str):
            reference_ligand_pdb = str(helper_data["reference_ligand_pdb"])
        ligand_preview_text = ligand_preview.text if ligand_preview is not None else None
        box = show_box_panel(ligand_preview_text, reference_ligand_pdb)
        options = show_options_panel(backend)
        run_panel(receptor_input, ligand_inputs, box, backend, runs_dir, options, helper_data)

    with right:
        preview_panel(receptor_input, helper_data, ligand_preview, box)

    st.markdown("---")
    current_tab, batch_tab, history_tab = st.tabs(["Current Run", "Batch Summary", "History"])
    with current_tab:
        render_run_details(st.session_state.get("last_run"), context_key="current")
    with batch_tab:
        render_batch_summary()
    with history_tab:
        render_history_panel(runs_dir)


if __name__ == "__main__":
    main()
