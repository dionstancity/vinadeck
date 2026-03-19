from __future__ import annotations

from pathlib import Path

from dockgui.models import DockingBox, DockingPose

ATOM_RECORDS = ("ATOM", "HETATM")
TWO_LETTER_ELEMENTS = {"BR", "CL", "FE", "MG", "ZN", "NA", "CA", "MN", "CU", "CO"}


def sanitize_filename(name: str) -> str:
    sanitized = "".join(char if char.isalnum() or char in {"-", "_", "."} else "_" for char in Path(name).name)
    return sanitized or "structure.pdbqt"


def iter_atom_lines(text: str):
    for line in text.splitlines():
        if line.startswith(ATOM_RECORDS):
            yield line.rstrip("\n")


def _parse_coordinate(line: str, start: int, end: int) -> float:
    return float(line[start:end].strip())


def _parse_int(line: str, start: int, end: int, default: int = 1) -> int:
    try:
        return int(line[start:end].strip())
    except ValueError:
        return default


def _infer_element(atom_name: str) -> str:
    letters = "".join(char for char in atom_name if char.isalpha()).upper()
    if not letters:
        return "C"
    if len(letters) >= 2 and letters[:2] in TWO_LETTER_ELEMENTS:
        return letters[:2].title()
    return letters[0].upper()


def pdbqt_to_pdb(text: str) -> str:
    pdb_lines: list[str] = []
    for line in iter_atom_lines(text):
        record = line[0:6].strip() or "ATOM"
        serial = _parse_int(line, 6, 11)
        atom_name = line[12:16].strip() or "C"
        residue_name = (line[17:20].strip() or "LIG")[:3]
        chain_id = (line[21:22].strip() or "A")[:1]
        residue_id = _parse_int(line, 22, 26)
        x = _parse_coordinate(line, 30, 38)
        y = _parse_coordinate(line, 38, 46)
        z = _parse_coordinate(line, 46, 54)
        element = _infer_element(atom_name)
        pdb_lines.append(
            f"{record:<6}{serial:>5} {atom_name:<4} {residue_name:>3} {chain_id:1}{residue_id:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.00:>6.2f}{0.00:>6.2f}          {element:>2}"
        )
    pdb_lines.append("END")
    return "\n".join(pdb_lines) + "\n"


def structure_has_atoms(text: str) -> bool:
    return any(True for _ in iter_atom_lines(text))


def compute_bounds(text: str) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    points: list[tuple[float, float, float]] = []
    for line in iter_atom_lines(text):
        points.append(
            (
                _parse_coordinate(line, 30, 38),
                _parse_coordinate(line, 38, 46),
                _parse_coordinate(line, 46, 54),
            )
        )
    if not points:
        raise ValueError("未在结构文件中找到原子坐标。")

    xs = [point[0] for point in points]
    ys = [point[1] for point in points]
    zs = [point[2] for point in points]
    return (min(xs), min(ys), min(zs)), (max(xs), max(ys), max(zs))


def suggest_box_from_structure(text: str, padding: float = 8.0, minimum_size: float = 12.0) -> DockingBox:
    lower, upper = compute_bounds(text)
    size_x = max((upper[0] - lower[0]) + padding, minimum_size)
    size_y = max((upper[1] - lower[1]) + padding, minimum_size)
    size_z = max((upper[2] - lower[2]) + padding, minimum_size)
    return DockingBox(
        center_x=(lower[0] + upper[0]) / 2.0,
        center_y=(lower[1] + upper[1]) / 2.0,
        center_z=(lower[2] + upper[2]) / 2.0,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
    )


def parse_vina_poses(text: str) -> list[DockingPose]:
    poses: list[DockingPose] = []
    current_lines: list[str] = []
    current_mode: int | None = None
    affinity = 0.0
    rmsd_lb = 0.0
    rmsd_ub = 0.0

    def flush_current() -> None:
        nonlocal affinity, current_lines, current_mode, rmsd_lb, rmsd_ub
        if current_lines and current_mode is not None:
            poses.append(
                DockingPose(
                    mode=current_mode,
                    affinity=affinity,
                    rmsd_lb=rmsd_lb,
                    rmsd_ub=rmsd_ub,
                    pdbqt_text="\n".join(current_lines) + "\n",
                )
            )
        current_lines = []
        current_mode = None
        affinity = 0.0
        rmsd_lb = 0.0
        rmsd_ub = 0.0

    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        if line.startswith("MODEL"):
            flush_current()
            current_lines = [line]
            parts = line.split()
            current_mode = int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else len(poses) + 1
            continue

        if current_mode is None:
            continue

        current_lines.append(line)
        if line.startswith("REMARK VINA RESULT:"):
            parts = line.split()
            if len(parts) >= 6:
                affinity = float(parts[3])
                rmsd_lb = float(parts[4])
                rmsd_ub = float(parts[5])
        elif line.startswith("ENDMDL"):
            flush_current()

    flush_current()
    return poses
