#!/usr/bin/env python3
"""Annotate E6 engineering features without scoring or ranking.

This script merges optional engineering feature annotations into the E6 input
contract. It does not calculate engineering scores or rank candidates.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd


DEFAULT_E6_INPUT_CSV = (
    "data/yrbp_branch/E6_engineering/"
    "E6_engineering_scoring_input_v0.1.csv"
)
DEFAULT_OUT_CSV = (
    "data/yrbp_branch/E6_engineering/"
    "E6_engineering_feature_annotated_v0.1.csv"
)

E6_REQUIRED_FIELDS = [
    "candidate_id",
    "accession",
    "receptor_id",
    "receptor_format",
    "receptor_type",
    "docking_result_status",
    "docking_hit_count",
    "docking_top_mode_rank",
    "docking_top_affinity_kcal_mol",
    "docking_top_rmsd_lb",
    "docking_top_rmsd_ub",
    "docking_affinity_mean",
    "docking_affinity_min",
    "docking_affinity_max",
    "docking_affinity_sd",
    "docking_evidence_ready",
    "expression_feature_status",
    "chaperone_feature_status",
    "immobilization_feature_status",
    "expression_feasibility_score",
    "chaperone_independence_score",
    "immobilization_feasibility_score",
    "engineering_scoring_ready",
    "recommended_next_step",
    "notes",
]

EXPRESSION_COLUMNS = [
    "sequence_length_aa",
    "molecular_weight_kda",
    "cysteine_count",
    "predicted_tm_helix_count",
    "signal_peptide_predicted",
    "low_complexity_fraction",
    "predicted_solubility_label",
    "expression_risk_label",
    "expression_annotation_source",
    "expression_notes",
]

CHAPERONE_COLUMNS = [
    "known_chaperone_required",
    "chaperone_evidence_source",
    "complex_assembly_flag",
    "tail_fiber_or_rbp_complex_flag",
    "oligomerization_risk_label",
    "chaperone_dependency_risk_label",
    "chaperone_notes",
]

IMMOBILIZATION_COLUMNS = [
    "lysine_count",
    "cysteine_handle_count",
    "n_terminal_handle_available",
    "c_terminal_handle_available",
    "predicted_pi",
    "surface_lysine_proxy",
    "immobilization_handle_label",
    "immobilization_risk_label",
    "immobilization_notes",
]

ANNOTATION_COLUMNS = (
    ["candidate_id", "accession"]
    + EXPRESSION_COLUMNS
    + CHAPERONE_COLUMNS
    + IMMOBILIZATION_COLUMNS
)

EXPRESSION_CORE_FIELDS = [
    "sequence_length_aa",
    "molecular_weight_kda",
    "predicted_solubility_label",
    "expression_risk_label",
]

CHAPERONE_CORE_FIELDS = [
    "known_chaperone_required",
    "complex_assembly_flag",
    "chaperone_dependency_risk_label",
]

IMMOBILIZATION_CORE_FIELDS = [
    "lysine_count",
    "cysteine_handle_count",
    "immobilization_handle_label",
    "immobilization_risk_label",
]

BOOLEAN_COLUMNS = [
    "signal_peptide_predicted",
    "known_chaperone_required",
    "complex_assembly_flag",
    "tail_fiber_or_rbp_complex_flag",
    "n_terminal_handle_available",
    "c_terminal_handle_available",
]

NUMERIC_COLUMNS = [
    "sequence_length_aa",
    "molecular_weight_kda",
    "cysteine_count",
    "predicted_tm_helix_count",
    "low_complexity_fraction",
    "lysine_count",
    "cysteine_handle_count",
    "predicted_pi",
    "surface_lysine_proxy",
]

FEATURE_COLUMNS = EXPRESSION_COLUMNS + CHAPERONE_COLUMNS + IMMOBILIZATION_COLUMNS
OUTPUT_FIELDS = E6_REQUIRED_FIELDS + [
    field for field in FEATURE_COLUMNS if field not in E6_REQUIRED_FIELDS
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge E6 engineering feature annotations."
    )
    parser.add_argument(
        "--e6_input_csv",
        default=DEFAULT_E6_INPUT_CSV,
        help="Input E6 engineering scoring contract CSV.",
    )
    parser.add_argument(
        "--annotation_csv",
        default="",
        help="Optional engineering feature annotation CSV.",
    )
    parser.add_argument(
        "--out_csv",
        default=DEFAULT_OUT_CSV,
        help="Output E6 engineering feature annotated CSV.",
    )
    return parser.parse_args()


def read_csv_with_fields(path: Path, required_fields: List[str]) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=required_fields)

    df = pd.read_csv(path, dtype=str, keep_default_na=False)
    for field in required_fields:
        if field not in df.columns:
            df[field] = ""
    return df.fillna("")


def read_annotation(path_text: str) -> pd.DataFrame:
    if not path_text.strip():
        return pd.DataFrame(columns=ANNOTATION_COLUMNS)
    return read_csv_with_fields(Path(path_text), ANNOTATION_COLUMNS)


def parse_bool_annotation(value: object) -> str:
    text = str(value).strip()
    lowered = text.lower()
    if not text:
        return ""
    if lowered in {"true", "1", "yes", "y"}:
        return "True"
    if lowered in {"false", "0", "no", "n"}:
        return "False"
    return text


def parse_float(value: object) -> Optional[float]:
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def format_float(value: Optional[float]) -> str:
    if value is None:
        return ""
    return f"{value:g}"


def normalize_numeric(value: object) -> str:
    text = str(value).strip()
    if not text:
        return ""
    parsed = parse_float(text)
    if parsed is None:
        return text
    return format_float(parsed)


def is_nonempty(value: object) -> bool:
    return bool(str(value).strip())


def append_notes(*parts: object) -> str:
    clean_parts = [str(part).strip() for part in parts if str(part).strip()]
    return "; ".join(clean_parts)


def normalize_annotation_row(row: pd.Series) -> Dict[str, str]:
    normalized: Dict[str, str] = {}
    for field in ANNOTATION_COLUMNS:
        value = row.get(field, "")
        if field in BOOLEAN_COLUMNS:
            normalized[field] = parse_bool_annotation(value)
        elif field in NUMERIC_COLUMNS:
            normalized[field] = normalize_numeric(value)
        else:
            normalized[field] = str(value).strip()
    return normalized


def build_annotation_indexes(
    annotation_df: pd.DataFrame,
) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]], set, set]:
    candidate_index: Dict[str, Dict[str, str]] = {}
    accession_index: Dict[str, Dict[str, str]] = {}
    duplicate_candidates = set()
    duplicate_accessions = set()

    for _, row in annotation_df.iterrows():
        normalized = normalize_annotation_row(row)
        candidate_id = normalized.get("candidate_id", "")
        accession = normalized.get("accession", "")

        if candidate_id:
            if candidate_id in candidate_index:
                duplicate_candidates.add(candidate_id)
            else:
                candidate_index[candidate_id] = normalized

        if accession:
            if accession in accession_index:
                duplicate_accessions.add(accession)
            else:
                accession_index[accession] = normalized

    return candidate_index, accession_index, duplicate_candidates, duplicate_accessions


def find_annotation(
    row: pd.Series,
    candidate_index: Dict[str, Dict[str, str]],
    accession_index: Dict[str, Dict[str, str]],
    duplicate_candidates: set,
    duplicate_accessions: set,
) -> Tuple[Dict[str, str], str]:
    candidate_id = str(row.get("candidate_id", "")).strip()
    accession = str(row.get("accession", "")).strip()

    if candidate_id and candidate_id in candidate_index:
        warning = "duplicate_annotation_warning" if candidate_id in duplicate_candidates else ""
        return candidate_index[candidate_id], warning

    if accession and accession in accession_index:
        warning = "duplicate_annotation_warning" if accession in duplicate_accessions else ""
        return accession_index[accession], warning

    return {}, ""


def feature_status(row: Dict[str, str], fields: List[str], core_fields: List[str]) -> str:
    has_any = any(is_nonempty(row.get(field, "")) for field in fields)
    has_all_core = all(is_nonempty(row.get(field, "")) for field in core_fields)
    if has_all_core:
        return "annotated"
    if has_any:
        return "partial"
    return "not_annotated"


def make_output_row(
    row: pd.Series,
    annotation: Dict[str, str],
    annotation_warning: str,
) -> Dict[str, str]:
    output = {field: str(row.get(field, "")).strip() for field in E6_REQUIRED_FIELDS}

    for field in FEATURE_COLUMNS:
        output[field] = annotation.get(field, "")

    expression_status = feature_status(
        output, EXPRESSION_COLUMNS, EXPRESSION_CORE_FIELDS
    )
    chaperone_status = feature_status(
        output, CHAPERONE_COLUMNS, CHAPERONE_CORE_FIELDS
    )
    immobilization_status = feature_status(
        output, IMMOBILIZATION_COLUMNS, IMMOBILIZATION_CORE_FIELDS
    )
    engineering_ready = (
        expression_status == "annotated"
        and chaperone_status == "annotated"
        and immobilization_status == "annotated"
    )

    output["expression_feature_status"] = expression_status
    output["chaperone_feature_status"] = chaperone_status
    output["immobilization_feature_status"] = immobilization_status
    output["expression_feasibility_score"] = str(
        row.get("expression_feasibility_score", "")
    ).strip()
    output["chaperone_independence_score"] = str(
        row.get("chaperone_independence_score", "")
    ).strip()
    output["immobilization_feasibility_score"] = str(
        row.get("immobilization_feasibility_score", "")
    ).strip()
    output["engineering_scoring_ready"] = str(engineering_ready)
    if engineering_ready:
        output["recommended_next_step"] = "run_e6_engineering_scoring"
    else:
        output["recommended_next_step"] = (
            "complete_expression_chaperone_immobilization_annotations"
        )

    annotation_notes = append_notes(
        annotation.get("expression_notes", ""),
        annotation.get("chaperone_notes", ""),
        annotation.get("immobilization_notes", ""),
        annotation_warning,
    )
    output["notes"] = append_notes(row.get("notes", ""), annotation_notes)
    return output


def write_output(rows: List[Dict[str, str]], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(rows, columns=OUTPUT_FIELDS)
    out_df = out_df.fillna("")
    out_df.to_csv(out_csv, index=False, encoding="utf-8", lineterminator="\n")


def main() -> int:
    args = parse_args()
    e6_input_csv = Path(args.e6_input_csv)
    out_csv = Path(args.out_csv)

    e6_df = read_csv_with_fields(e6_input_csv, E6_REQUIRED_FIELDS)
    annotation_df = read_annotation(args.annotation_csv)
    (
        candidate_index,
        accession_index,
        duplicate_candidates,
        duplicate_accessions,
    ) = build_annotation_indexes(annotation_df)

    output_rows: List[Dict[str, str]] = []
    fully_annotated_count = 0
    for _, row in e6_df.iterrows():
        annotation, warning = find_annotation(
            row,
            candidate_index,
            accession_index,
            duplicate_candidates,
            duplicate_accessions,
        )
        output_row = make_output_row(row, annotation, warning)
        if output_row["engineering_scoring_ready"] == "True":
            fully_annotated_count += 1
        output_rows.append(output_row)

    write_output(output_rows, out_csv)

    print(f"Read {len(e6_df)} E6 input rows.")
    print(f"Read {len(annotation_df)} annotation rows.")
    print(f"Wrote {len(output_rows)} annotated engineering rows.")
    print(f"Fully annotated rows: {fully_annotated_count}.")
    print(f"Output path: {out_csv}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
