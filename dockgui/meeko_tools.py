from __future__ import annotations

import io
import os
import sys
from contextlib import redirect_stderr, redirect_stdout
from dataclasses import dataclass
from functools import lru_cache
from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from tempfile import TemporaryDirectory

from dockgui.parsing import sanitize_filename

LIGAND_FORMATS = {"sdf", "mol2", "mol"}
RECEPTOR_FORMATS = {"pdb", "pqr"}
LIGAND_CHARGE_MODELS = {"gasteiger", "zero", "read"}
RECEPTOR_CHARGE_MODELS = {"gasteiger", "zero", "read"}


@dataclass(slots=True, frozen=True)
class MeekoPreparedFile:
    name: str
    data: bytes
    stdout: str
    stderr: str


@lru_cache(maxsize=1)
def probe_meeko() -> tuple[bool, str]:
    try:
        import_module("meeko.cli.mk_prepare_ligand")
        import_module("meeko.cli.mk_prepare_receptor")
        try:
            meeko_version = version("meeko")
        except PackageNotFoundError:
            meeko_version = "unknown"
        return True, f"Meeko {meeko_version} is available"
    except Exception as exc:
        return False, f"{type(exc).__name__}: {exc}"


def _normalize_exit_code(code: object) -> int:
    if code is None:
        return 0
    if isinstance(code, int):
        return code
    return 1


def _run_cli(main_func, argv: list[str], *, cwd: Path) -> tuple[int, str, str]:
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    original_argv = sys.argv[:]
    original_cwd = Path.cwd()
    exit_code = 0

    try:
        sys.argv = argv
        os.chdir(cwd)
        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            try:
                main_func()
            except SystemExit as exc:
                exit_code = _normalize_exit_code(exc.code)
            except Exception as exc:
                print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
                exit_code = 1
    finally:
        sys.argv = original_argv
        os.chdir(original_cwd)

    return exit_code, stdout_buffer.getvalue().strip(), stderr_buffer.getvalue().strip()


def _format_failure(stdout: str, stderr: str, fallback: str) -> str:
    message = "\n".join(part for part in (stderr, stdout) if part).strip()
    return message or fallback


def _looks_like_missing_partial_charges(message: str) -> bool:
    normalized = message.lower()
    return (
        "partialcharge contains none" in normalized
        or "the list of charges based on atom property name partialcharge contains none" in normalized
    )


def _build_ligand_name(
    *,
    output_stem: str,
    molecule_name: str,
    index: int,
    setup_index: int,
    total_setups: int,
    seen_names: set[str],
) -> str:
    candidate_stem = sanitize_filename(molecule_name).strip()
    if candidate_stem.lower().endswith(".pdbqt"):
        candidate_stem = Path(candidate_stem).stem
    if not candidate_stem:
        candidate_stem = f"{output_stem}-{index}"
    if total_setups > 1:
        candidate_stem = f"{candidate_stem}-setup-{setup_index}"

    candidate_name = f"{candidate_stem}.pdbqt"
    counter = 2
    while candidate_name.lower() in seen_names:
        candidate_name = f"{candidate_stem}-{counter}.pdbqt"
        counter += 1
    seen_names.add(candidate_name.lower())
    return candidate_name


def _summarize_failures(failures: list[str], *, prepared_count: int) -> str:
    if not failures:
        return ""
    lines = [f"Prepared {prepared_count} ligands; skipped {len(failures)} molecules."]
    lines.extend(failures[:20])
    if len(failures) > 20:
        lines.append(f"... {len(failures) - 20} more failures omitted.")
    return "\n".join(lines)


def _prepare_sdf_ligands_pdbqt(
    input_path: Path,
    *,
    output_stem: str,
    charge_model: str,
) -> list[MeekoPreparedFile]:
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.MolStandardize import rdMolStandardize

    charge_atom_prop = "PartialCharge" if charge_model == "read" else None
    preparator = MoleculePreparation(charge_model=charge_model, charge_atom_prop=charge_atom_prop)
    fragment_chooser = rdMolStandardize.LargestFragmentChooser()
    supplier = Chem.SDMolSupplier(str(input_path), removeHs=False, sanitize=True, strictParsing=False)
    molecules = list(supplier)
    del supplier

    prepared_files: list[MeekoPreparedFile] = []
    failures: list[str] = []
    seen_names: set[str] = set()

    for index, molecule in enumerate(molecules, start=1):
        if molecule is None:
            failures.append(f"Molecule {index}: unreadable structure in input SDF.")
            continue

        molecule_name = ""
        if molecule.HasProp("_Name"):
            molecule_name = molecule.GetProp("_Name").strip()
        if not molecule_name:
            molecule_name = f"{output_stem}-{index}"

        try:
            working_molecule = fragment_chooser.choose(Chem.Mol(molecule))
            has_3d = working_molecule.GetNumConformers() > 0 and working_molecule.GetConformer().Is3D()
            working_molecule = Chem.AddHs(working_molecule, addCoords=has_3d)

            if working_molecule.GetNumConformers() == 0 or not working_molecule.GetConformer().Is3D():
                params = AllChem.ETKDGv3()
                params.randomSeed = 42 + index
                status = AllChem.EmbedMolecule(working_molecule, params)
                if status != 0:
                    status = AllChem.EmbedMolecule(
                        working_molecule,
                        maxAttempts=1000,
                        randomSeed=42 + index,
                        useRandomCoords=True,
                        ignoreSmoothingFailures=True,
                    )
                if status != 0:
                    raise RuntimeError("3D embedding failed")
                try:
                    if AllChem.MMFFHasAllMoleculeParams(working_molecule):
                        AllChem.MMFFOptimizeMolecule(working_molecule)
                    else:
                        AllChem.UFFOptimizeMolecule(working_molecule)
                except Exception:
                    pass

            setups = preparator.prepare(working_molecule)
            if not setups:
                raise RuntimeError("Meeko did not return any ligand setup.")

            for setup_index, setup in enumerate(setups, start=1):
                pdbqt_string, success, error_message = PDBQTWriterLegacy.write_string(setup)
                if not success or not pdbqt_string.strip():
                    raise RuntimeError(error_message or "Failed to generate ligand PDBQT.")

                prepared_files.append(
                    MeekoPreparedFile(
                        name=_build_ligand_name(
                            output_stem=output_stem,
                            molecule_name=molecule_name,
                            index=index,
                            setup_index=setup_index,
                            total_setups=len(setups),
                            seen_names=seen_names,
                        ),
                        data=pdbqt_string.encode("utf-8"),
                        stdout="Prepared via RDKit + Meeko API.",
                        stderr="",
                    )
                )
        except Exception as exc:
            failures.append(f"Molecule {index} ({molecule_name}): {exc}")

    failure_summary = _summarize_failures(failures, prepared_count=len(prepared_files))
    if not prepared_files:
        if charge_model == "read" and _looks_like_missing_partial_charges(failure_summary):
            raise RuntimeError(
                "Charge model 'read' requires existing per-atom partial charges in the input SDF/MOL2. "
                "This file does not contain usable partial charges. Switch to 'gasteiger' for standard ligand libraries."
            )
        raise RuntimeError(failure_summary or "Meeko did not generate any ligand PDBQT files.")

    if not failure_summary:
        return prepared_files

    return [
        MeekoPreparedFile(
            name=prepared_file.name,
            data=prepared_file.data,
            stdout=prepared_file.stdout,
            stderr=failure_summary,
        )
        for prepared_file in prepared_files
    ]


def prepare_ligands_pdbqt(
    ligand_name: str,
    ligand_bytes: bytes,
    *,
    charge_model: str = "gasteiger",
) -> list[MeekoPreparedFile]:
    available, detail = probe_meeko()
    if not available:
        raise RuntimeError(f"Meeko is not available: {detail}")

    ligand_suffix = Path(ligand_name).suffix.lower().lstrip(".")
    if ligand_suffix not in LIGAND_FORMATS:
        supported = ", ".join(sorted(LIGAND_FORMATS))
        raise ValueError(f"Ligand format must be one of: {supported}")
    if charge_model not in LIGAND_CHARGE_MODELS:
        supported = ", ".join(sorted(LIGAND_CHARGE_MODELS))
        raise ValueError(f"Unsupported ligand charge model: {charge_model}. Choices: {supported}")

    safe_input_name = sanitize_filename(ligand_name)
    output_stem = Path(safe_input_name).stem or "ligand"

    with TemporaryDirectory(prefix="vinadock-ligand-") as temp_dir:
        workdir = Path(temp_dir)
        input_path = workdir / safe_input_name
        output_dir = workdir / "prepared_ligands"
        input_path.write_bytes(ligand_bytes)

        if ligand_suffix == "sdf":
            return _prepare_sdf_ligands_pdbqt(
                input_path,
                output_stem=output_stem,
                charge_model=charge_model,
            )

        module = import_module("meeko.cli.mk_prepare_ligand")
        command = [
            "mk_prepare_ligand.py",
            "-i",
            input_path.name,
            "--multimol_outdir",
            output_dir.name,
            "--multimol_prefix",
            output_stem,
            "--charge_model",
            charge_model,
        ]
        exit_code, stdout, stderr = _run_cli(module.main, command, cwd=workdir)

        output_paths = sorted(output_dir.glob("*.pdbqt"))
        if exit_code != 0 and not output_paths:
            failure_message = _format_failure(stdout, stderr, "Meeko ligand preparation failed.")
            if charge_model == "read" and _looks_like_missing_partial_charges(failure_message):
                raise RuntimeError(
                    "Charge model 'read' requires existing per-atom partial charges in the input SDF/MOL2. "
                    "This file does not contain usable partial charges. Switch to 'gasteiger' for standard ligand libraries."
                )
            raise RuntimeError(failure_message)
        if not output_paths:
            raise RuntimeError("Meeko did not generate any ligand PDBQT files.")

        prepared_files: list[MeekoPreparedFile] = []
        for output_path in output_paths:
            display_name = output_path.name
            if len(output_paths) == 1:
                display_name = f"{output_stem}.pdbqt"
            prepared_files.append(
                MeekoPreparedFile(
                    name=display_name,
                    data=output_path.read_bytes(),
                    stdout=stdout,
                    stderr=stderr,
                )
            )

        return prepared_files


def prepare_ligand_pdbqt(
    ligand_name: str,
    ligand_bytes: bytes,
    *,
    charge_model: str = "gasteiger",
) -> MeekoPreparedFile:
    prepared_files = prepare_ligands_pdbqt(
        ligand_name,
        ligand_bytes,
        charge_model=charge_model,
    )
    if len(prepared_files) != 1:
        raise RuntimeError("Input contains multiple ligands; use the batch ligand preparation flow instead.")
    return prepared_files[0]


def prepare_receptor_pdbqt(
    receptor_name: str,
    receptor_bytes: bytes,
    *,
    input_format: str,
    allow_bad_res: bool = True,
    charge_model: str = "gasteiger",
) -> MeekoPreparedFile:
    available, detail = probe_meeko()
    if not available:
        raise RuntimeError(f"Meeko is not available: {detail}")

    normalized_format = input_format.lower().lstrip(".")
    if normalized_format not in RECEPTOR_FORMATS:
        supported = ", ".join(sorted(RECEPTOR_FORMATS))
        raise ValueError(f"Receptor format must be one of: {supported}")
    if charge_model not in RECEPTOR_CHARGE_MODELS:
        supported = ", ".join(sorted(RECEPTOR_CHARGE_MODELS))
        raise ValueError(f"Unsupported receptor charge model: {charge_model}. Choices: {supported}")
    if charge_model == "read" and normalized_format != "pqr":
        raise ValueError("Receptor charge model 'read' only works with PQR input.")

    module = import_module("meeko.cli.mk_prepare_receptor")
    safe_input_name = sanitize_filename(receptor_name)
    output_name = f"{Path(safe_input_name).stem}.pdbqt"
    read_flag = "--read_pdb" if normalized_format == "pdb" else "--read_pqr"

    with TemporaryDirectory(prefix="vinadock-receptor-") as temp_dir:
        workdir = Path(temp_dir)
        input_path = workdir / safe_input_name
        output_path = workdir / output_name
        input_path.write_bytes(receptor_bytes)

        command = [
            "mk_prepare_receptor.py",
            read_flag,
            input_path.name,
            "-p",
            output_path.name,
            "--charge_model",
            charge_model,
        ]
        if allow_bad_res:
            command.append("-a")

        exit_code, stdout, stderr = _run_cli(module.main, command, cwd=workdir)

        if exit_code != 0:
            raise RuntimeError(_format_failure(stdout, stderr, "Meeko receptor preparation failed."))
        if not output_path.exists():
            raise RuntimeError("Meeko did not generate a receptor PDBQT file.")

        return MeekoPreparedFile(
            name=output_name,
            data=output_path.read_bytes(),
            stdout=stdout,
            stderr=stderr,
        )
