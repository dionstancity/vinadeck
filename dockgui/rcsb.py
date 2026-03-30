from __future__ import annotations

import json
from collections import OrderedDict
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from dockgui.models import PdbEntryMetadata, PdbResidue

ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
STRUCTURE_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
ATOM_RECORDS = {"ATOM", "HETATM"}
KEEP_ALT_LOCS = {"", "A", "1"}
WATER_RESIDUES = {"HOH", "WAT", "DOD"}
HEADER_RECORDS = ("HEADER", "TITLE ", "COMPND", "SOURCE", "KEYWDS", "EXPDTA", "AUTHOR")


def normalize_pdb_id(value: str) -> str:
    pdb_id = value.strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise ValueError("PDB ID 必须是 4 位字母或数字，例如 `1STP`。")
    return pdb_id


def _open_url(url: str):
    request = Request(url, headers={"User-Agent": "VinaDock Studio/1.0"})
    return urlopen(request, timeout=20)


def _fetch_json(url: str) -> dict:
    try:
        with _open_url(url) as response:
            return json.load(response)
    except HTTPError as exc:
        if exc.code == 404:
            raise ValueError("指定的 PDB ID 不存在。") from exc
        raise RuntimeError(f"访问 RCSB 失败：HTTP {exc.code}") from exc
    except URLError as exc:
        raise RuntimeError("无法连接 RCSB，请检查网络连接。") from exc


def _fetch_text(url: str) -> str:
    try:
        with _open_url(url) as response:
            return response.read().decode("utf-8", errors="ignore")
    except HTTPError as exc:
        if exc.code == 404:
            raise ValueError("指定的 PDB 结构文件不存在。") from exc
        raise RuntimeError(f"下载结构文件失败：HTTP {exc.code}") from exc
    except URLError as exc:
        raise RuntimeError("无法下载结构文件，请检查网络连接。") from exc


def _safe_int(value: str, default: int = 0) -> int:
    try:
        return int(value.strip())
    except ValueError:
        return default


def _safe_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value.strip())
    except ValueError:
        return default


def _infer_element(atom_name: str, element: str) -> str:
    if element:
        return element.upper()
    letters = "".join(char for char in atom_name if char.isalpha()).upper()
    if not letters:
        return ""
    return letters[:2].strip()


def iter_primary_model_lines(text: str):
    has_model = False
    inside_first_model = False
    for raw_line in text.splitlines():
        line = raw_line.rstrip("\n")
        record = line[0:6].strip()
        if record == "MODEL":
            if not has_model:
                has_model = True
                inside_first_model = True
            continue
        if record == "ENDMDL":
            if inside_first_model:
                break
            continue
        if has_model and not inside_first_model:
            continue
        yield line


def parse_atom_line(line: str) -> dict[str, str | int | float] | None:
    record = line[0:6].strip()
    if record not in ATOM_RECORDS:
        return None
    atom_name = line[12:16].strip()
    alt_loc = line[16:17].strip()
    res_name = line[17:20].strip()
    chain_id = line[21:22].strip()
    res_seq = _safe_int(line[22:26], 0)
    insertion_code = line[26:27].strip()
    x = _safe_float(line[30:38])
    y = _safe_float(line[38:46])
    z = _safe_float(line[46:54])
    element = _infer_element(atom_name, line[76:78].strip())
    return {
        "record": record,
        "atom_name": atom_name,
        "alt_loc": alt_loc,
        "res_name": res_name,
        "chain_id": chain_id,
        "res_seq": res_seq,
        "insertion_code": insertion_code,
        "x": x,
        "y": y,
        "z": z,
        "element": element,
        "line": line,
    }


def is_hydrogen(atom: dict[str, str | int | float]) -> bool:
    return str(atom["element"]).upper() == "H" or str(atom["atom_name"]).upper().startswith("H")


def residue_token(res_name: str, chain_id: str, res_seq: int, insertion_code: str = "") -> str:
    chain_value = chain_id or "_"
    insertion = insertion_code or ""
    return f"{res_name}:{chain_value}:{res_seq}{insertion}"


def summarize_pdb_structure(text: str) -> tuple[list[str], list[PdbResidue]]:
    chains: OrderedDict[str, None] = OrderedDict()
    residues: OrderedDict[str, PdbResidue] = OrderedDict()
    counts: dict[str, int] = {}

    for line in iter_primary_model_lines(text):
        atom = parse_atom_line(line)
        if atom is None:
            continue
        if str(atom["alt_loc"]) not in KEEP_ALT_LOCS:
            continue

        chain_id = str(atom["chain_id"])
        chains.setdefault(chain_id, None)

        if atom["record"] != "HETATM":
            continue
        if str(atom["res_name"]).upper() in WATER_RESIDUES:
            continue

        token = residue_token(
            str(atom["res_name"]),
            chain_id,
            int(atom["res_seq"]),
            str(atom["insertion_code"]),
        )
        counts[token] = counts.get(token, 0) + 1
        if token not in residues:
            residues[token] = PdbResidue(
                res_name=str(atom["res_name"]),
                chain_id=chain_id,
                res_seq=int(atom["res_seq"]),
                insertion_code=str(atom["insertion_code"]),
                atom_count=0,
            )

    ligand_list = [
        PdbResidue(
            res_name=residue.res_name,
            chain_id=residue.chain_id,
            res_seq=residue.res_seq,
            insertion_code=residue.insertion_code,
            atom_count=counts[residue.token],
        )
        for residue in residues.values()
    ]
    ligand_list.sort(key=lambda residue: (residue.chain_id, residue.res_name, residue.res_seq, residue.insertion_code))
    return list(chains.keys()), ligand_list


def fetch_pdb_entry(pdb_id: str) -> tuple[PdbEntryMetadata, str]:
    normalized_id = normalize_pdb_id(pdb_id)
    entry_data = _fetch_json(ENTRY_URL.format(pdb_id=normalized_id))
    pdb_text = _fetch_text(STRUCTURE_URL.format(pdb_id=normalized_id))
    chains, ligands = summarize_pdb_structure(pdb_text)

    entry_info = entry_data.get("rcsb_entry_info", {})
    refine_data = entry_data.get("refine", [{}])
    resolution = None
    resolution_values = entry_info.get("resolution_combined") or []
    if resolution_values:
        resolution = float(resolution_values[0])
    elif refine_data and refine_data[0].get("ls_dres_high") is not None:
        resolution = float(refine_data[0]["ls_dres_high"])

    metadata = PdbEntryMetadata(
        pdb_id=normalized_id,
        title=str(entry_data.get("struct", {}).get("title", normalized_id)),
        experimental_method=str(
            entry_info.get("experimental_method")
            or (entry_data.get("exptl", [{}])[0].get("method") if entry_data.get("exptl") else "")
        ),
        resolution=resolution,
        chains=chains,
        ligands=ligands,
    )
    return metadata, pdb_text


def clean_pdb_text(
    text: str,
    selected_chains: list[str] | None,
    kept_residue_tokens: set[str] | None,
    *,
    remove_water: bool = True,
    remove_hydrogen: bool = True,
) -> tuple[str, dict[str, int]]:
    chain_filter = set(selected_chains) if selected_chains else None
    kept_tokens = kept_residue_tokens or set()
    header_lines: list[str] = []
    atom_lines: list[str] = []
    stats = {
        "kept_protein_atoms": 0,
        "kept_cofactor_atoms": 0,
        "removed_chain_atoms": 0,
        "removed_water_atoms": 0,
        "removed_hetero_atoms": 0,
        "removed_hydrogen_atoms": 0,
        "removed_altloc_atoms": 0,
    }

    for line in iter_primary_model_lines(text):
        if line.startswith(HEADER_RECORDS):
            header_lines.append(line)
            continue

        atom = parse_atom_line(line)
        if atom is None:
            continue

        if str(atom["alt_loc"]) not in KEEP_ALT_LOCS:
            stats["removed_altloc_atoms"] += 1
            continue

        if chain_filter is not None and str(atom["chain_id"]) not in chain_filter:
            stats["removed_chain_atoms"] += 1
            continue

        if remove_hydrogen and is_hydrogen(atom):
            stats["removed_hydrogen_atoms"] += 1
            continue

        if atom["record"] == "HETATM":
            current_residue_token = residue_token(
                str(atom["res_name"]),
                str(atom["chain_id"]),
                int(atom["res_seq"]),
                str(atom["insertion_code"]),
            )
            if str(atom["res_name"]).upper() in WATER_RESIDUES and remove_water:
                stats["removed_water_atoms"] += 1
                continue
            if current_residue_token not in kept_tokens:
                stats["removed_hetero_atoms"] += 1
                continue
            stats["kept_cofactor_atoms"] += 1
        else:
            stats["kept_protein_atoms"] += 1

        atom_lines.append(str(atom["line"]))

    cleaned_lines = header_lines + atom_lines + ["END"]
    return "\n".join(cleaned_lines) + "\n", stats


def extract_residue_pdb(text: str, residue: PdbResidue, *, remove_hydrogen: bool = True) -> str:
    residue_lines: list[str] = []
    target_token = residue.token
    for line in iter_primary_model_lines(text):
        atom = parse_atom_line(line)
        if atom is None:
            continue
        if str(atom["alt_loc"]) not in KEEP_ALT_LOCS:
            continue
        current_token = residue_token(
            str(atom["res_name"]),
            str(atom["chain_id"]),
            int(atom["res_seq"]),
            str(atom["insertion_code"]),
        )
        if current_token != target_token:
            continue
        if remove_hydrogen and is_hydrogen(atom):
            continue
        residue_lines.append(str(atom["line"]))

    if not residue_lines:
        raise ValueError("未在结构中找到所选配体。")
    return "\n".join(residue_lines) + "\nEND\n"
