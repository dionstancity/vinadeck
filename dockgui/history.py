from __future__ import annotations

import json
from pathlib import Path

from dockgui.models import DockingBox, DockingRun, RunHistoryRecord
from dockgui.parsing import parse_vina_poses

MANIFEST_NAME = "run.json"


def save_run_manifest(run: DockingRun) -> Path:
    manifest_path = run.run_dir / MANIFEST_NAME
    manifest = {
        "schema_version": 1,
        "created_at": run.created_at,
        "status": run.status,
        "receptor_label": run.receptor_label or run.receptor_path.name,
        "ligand_label": run.ligand_label or run.ligand_path.name,
        "source_label": run.source_label,
        "engine_name": run.engine_name,
        "best_affinity": run.best_affinity,
        "pose_count": len(run.poses),
        "command": run.command,
        "stdout": run.stdout,
        "stderr": run.stderr,
        "docking_box": run.docking_box.as_dict(),
        "files": {
            "receptor": str(run.receptor_path.relative_to(run.run_dir)),
            "ligand": str(run.ligand_path.relative_to(run.run_dir)),
            "config": str(run.config_path.relative_to(run.run_dir)),
            "output": str(run.output_path.relative_to(run.run_dir)),
            "log": str(run.log_path.relative_to(run.run_dir)),
        },
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
    run.manifest_path = manifest_path
    return manifest_path


def list_history(base_dir: Path) -> list[RunHistoryRecord]:
    if not base_dir.exists():
        return []

    records: list[RunHistoryRecord] = []
    for manifest_path in base_dir.glob(f"run-*/{MANIFEST_NAME}"):
        try:
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue

        records.append(
            RunHistoryRecord(
                run_dir=manifest_path.parent,
                manifest_path=manifest_path,
                created_at=str(data.get("created_at", "")),
                receptor_label=str(data.get("receptor_label", "")),
                ligand_label=str(data.get("ligand_label", "")),
                best_affinity=float(data["best_affinity"]) if data.get("best_affinity") is not None else None,
                pose_count=int(data.get("pose_count", 0)),
                source_label=str(data.get("source_label", "")),
                engine_name=str(data.get("engine_name", "AutoDock Vina (CPU)")),
                status=str(data.get("status", "success")),
            )
        )

    return sorted(records, key=lambda record: record.created_at, reverse=True)


def load_run_from_manifest(manifest_path: Path) -> DockingRun:
    data = json.loads(manifest_path.read_text(encoding="utf-8"))
    run_dir = manifest_path.parent
    files = data.get("files", {})
    receptor_path = run_dir / str(files.get("receptor", "receptor.pdbqt"))
    ligand_path = run_dir / str(files.get("ligand", "ligand.pdbqt"))
    config_path = run_dir / str(files.get("config", "vina.conf"))
    output_path = run_dir / str(files.get("output", "poses_out.pdbqt"))
    log_path = run_dir / str(files.get("log", "docking.log"))
    poses = []
    if output_path.exists():
        poses = parse_vina_poses(output_path.read_text(encoding="utf-8", errors="ignore"))

    return DockingRun(
        run_dir=run_dir,
        receptor_path=receptor_path,
        ligand_path=ligand_path,
        config_path=config_path,
        output_path=output_path,
        log_path=log_path,
        stdout=str(data.get("stdout", "")),
        stderr=str(data.get("stderr", "")),
        command=[str(item) for item in data.get("command", [])],
        docking_box=DockingBox.from_dict(data.get("docking_box", {})),
        poses=poses,
        created_at=str(data.get("created_at", "")),
        receptor_label=str(data.get("receptor_label", receptor_path.name)),
        ligand_label=str(data.get("ligand_label", ligand_path.name)),
        source_label=str(data.get("source_label", "")),
        engine_name=str(data.get("engine_name", "AutoDock Vina (CPU)")),
        status=str(data.get("status", "success")),
        manifest_path=manifest_path,
    )
