from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True, frozen=True)
class DockingBackendConfig:
    kind: str
    executable: str
    kernel_executable: str = ""
    thread: int = 8000
    search_depth: int | None = None
    opencl_binary_path: str = ""
    rilc_bfgs: bool = True

    @property
    def label(self) -> str:
        labels = {
            "vina_cpu": "AutoDock Vina (CPU)",
            "vina_gpu_official": "Official Vina-GPU",
            "vina_gpu": "Vina-GPU",
            "vina_gpu_21": "Vina-GPU 2.1",
        }
        return labels.get(self.kind, self.kind)


@dataclass(slots=True)
class DockingBox:
    center_x: float
    center_y: float
    center_z: float
    size_x: float
    size_y: float
    size_z: float

    def as_dict(self) -> dict[str, float]:
        return {
            "center_x": self.center_x,
            "center_y": self.center_y,
            "center_z": self.center_z,
            "size_x": self.size_x,
            "size_y": self.size_y,
            "size_z": self.size_z,
        }

    @classmethod
    def from_dict(cls, values: dict[str, float]) -> "DockingBox":
        return cls(
            center_x=float(values["center_x"]),
            center_y=float(values["center_y"]),
            center_z=float(values["center_z"]),
            size_x=float(values["size_x"]),
            size_y=float(values["size_y"]),
            size_z=float(values["size_z"]),
        )


@dataclass(slots=True, frozen=True)
class StructureInput:
    name: str
    data: bytes
    file_format: str
    source_label: str = ""

    @property
    def text(self) -> str:
        return self.data.decode("utf-8", errors="ignore")


@dataclass(slots=True)
class DockingPose:
    mode: int
    affinity: float
    rmsd_lb: float
    rmsd_ub: float
    pdbqt_text: str

    def as_row(self) -> dict[str, float | int]:
        return {
            "Mode": self.mode,
            "Affinity (kcal/mol)": self.affinity,
            "RMSD l.b.": self.rmsd_lb,
            "RMSD u.b.": self.rmsd_ub,
        }


@dataclass(slots=True)
class DockingJob:
    receptor_name: str
    receptor_bytes: bytes
    ligand_name: str
    ligand_bytes: bytes
    docking_box: DockingBox
    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 3.0
    cpu: int = 0
    seed: int | None = None
    source_label: str = "Uploaded receptor PDBQT"
    engine_name: str = "AutoDock Vina (CPU)"


@dataclass(slots=True)
class DockingRun:
    run_dir: Path
    receptor_path: Path
    ligand_path: Path
    config_path: Path
    output_path: Path
    log_path: Path
    stdout: str
    stderr: str
    command: list[str]
    docking_box: DockingBox
    poses: list[DockingPose]
    created_at: str = ""
    receptor_label: str = ""
    ligand_label: str = ""
    source_label: str = ""
    engine_name: str = "AutoDock Vina (CPU)"
    status: str = "success"
    manifest_path: Path | None = None

    @property
    def best_affinity(self) -> float | None:
        return self.poses[0].affinity if self.poses else None

    def as_batch_row(self) -> dict[str, str | float | int]:
        return {
            "Ligand": self.ligand_label or self.ligand_path.name,
            "Best affinity (kcal/mol)": self.best_affinity if self.best_affinity is not None else "",
            "Poses": len(self.poses),
            "Engine": self.engine_name,
            "Run directory": str(self.run_dir),
            "Created": self.created_at,
        }


@dataclass(slots=True, frozen=True)
class PdbResidue:
    res_name: str
    chain_id: str
    res_seq: int
    insertion_code: str = ""
    atom_count: int = 0

    @property
    def token(self) -> str:
        chain_id = self.chain_id or "_"
        insertion = self.insertion_code or ""
        return f"{self.res_name}:{chain_id}:{self.res_seq}{insertion}"

    @property
    def label(self) -> str:
        chain_label = self.chain_id or "(blank)"
        insertion = self.insertion_code or ""
        return f"{self.res_name} | Chain {chain_label} | Residue {self.res_seq}{insertion} | {self.atom_count} atoms"


@dataclass(slots=True)
class PdbEntryMetadata:
    pdb_id: str
    title: str
    experimental_method: str
    resolution: float | None
    chains: list[str]
    ligands: list[PdbResidue]


@dataclass(slots=True)
class RunHistoryRecord:
    run_dir: Path
    manifest_path: Path
    created_at: str
    receptor_label: str
    ligand_label: str
    best_affinity: float | None
    pose_count: int
    source_label: str
    engine_name: str
    status: str

    @property
    def label(self) -> str:
        score = f"{self.best_affinity:.2f} kcal/mol" if self.best_affinity is not None else "n/a"
        return f"{self.created_at} | {self.ligand_label} | {score}"

    def as_row(self) -> dict[str, str | float | int]:
        return {
            "Created": self.created_at,
            "Ligand": self.ligand_label,
            "Receptor": self.receptor_label,
            "Best affinity (kcal/mol)": self.best_affinity if self.best_affinity is not None else "",
            "Poses": self.pose_count,
            "Engine": self.engine_name,
            "Source": self.source_label,
            "Status": self.status,
            "Run directory": str(self.run_dir),
        }
