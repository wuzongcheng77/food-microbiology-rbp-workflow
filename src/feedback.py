"""Deterministic feedback; extracted from archived v5.4.9g, registry-driven."""
import csv
import hashlib
import math
from collections import Counter, defaultdict
from pathlib import Path

def read_csv(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def file_sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest().upper()


def seq_sha(seq):
    return hashlib.sha256(seq.encode("utf-8")).hexdigest()


def num(v, default=None):
    if v in (None, ""):
        return default
    try:
        return float(v)
    except ValueError:
        return default


def fmt(v):
    return "" if v is None else f"{float(v):.6f}"


def b(v):
    return "true" if v else "false"


def scalar_sim(a, c):
    return None if a is None or c is None else max(0.0, min(1.0, 1.0 - abs(a - c)))


def minmax_sim(a, c, vals):
    if a is None or c is None:
        return None
    vals = [float(v) for v in vals if v is not None]
    if not vals:
        return None
    span = max(vals) - min(vals)
    return 1.0 if span == 0 else max(0.0, min(1.0, 1.0 - abs(a - c) / span))


def mean(vals):
    vals = [v for v in vals if v is not None]
    return None if not vals else sum(vals) / len(vals)


def weighted(parts, fallback=0.5):
    got = [(w, v) for _, w, v in parts if v is not None]
    if not got:
        return fallback, True, "low"
    score = sum(w * v for w, v in got) / sum(w for w, _ in got)
    norm = len(got) != len(parts)
    conf = "high" if not norm else ("medium" if len(got) >= 3 else "low")
    return score, norm, conf


def jaccard(a, c):
    return None if not a or not c else len(a & c) / len(a | c)


def register_verified(evidence):
    grouped = defaultdict(list)
    for row in evidence:
        if row['evidence_layer'] != 'C_verified':
            raise ValueError('Unexpected evidence layer')
        grouped[row['accession']].append(row)
    result = {}
    for candidate, records in grouped.items():
        if not all(r['evidence_status'] in ('pass','strong_pass') for r in records):
            raise ValueError('Unsupported evidence status; author review required')
        result[candidate] = 1.0
    return result

def reproduce_feedback(root):
    TOP50 = root/'inputs/frozen_top50.csv'
    E6 = root/'inputs/feedback/e6.csv'
    DOMAIN = root/'inputs/feedback/domains.csv'
    CONTACT = root/'inputs/feedback/contacts.csv'
    B0 = root/'inputs/feedback/sequences.csv'
    B1 = root/'inputs/feedback/no_optional_b1.csv'
    top50 = read_csv(TOP50)
    accessions = [r['source_accession'] for r in top50]
    evidence = read_csv(root/'inputs/feedback/registry.csv')
    status_scores = register_verified(evidence)
    if len(status_scores) != 1:
        raise ValueError('This frozen single-anchor policy requires exactly one verified candidate')
    ANCHOR = next(iter(status_scores))
    record_ids = {r["candidate_record_id"] for r in top50}
    seqs = {}
    for row in read_csv(B0):
        acc = row.get("source_accession")
        rec = row.get("candidate_record_id")
        seq = (row.get("sequence_normalized") or "").strip().upper()
        if seq and (acc in accessions or rec in record_ids):
            seqs[acc] = {
                "seq": seq,
                "source": B0.relative_to(root).as_posix(),
                "sha": row.get("sequence_sha256") or seq_sha(seq),
                "status": "recovered_from_B0_standard_candidate_records",
                "notes": f"candidate_record_id={rec}; identity_status={row.get('identity_status','')}",
            }
    if B1.exists():
        for row in read_csv(B1):
            acc = row.get("representative_accession")
            seq = (row.get("sequence_normalized") or "").strip().upper()
            if seq and acc in accessions and acc not in seqs:
                seqs[acc] = {
                    "seq": seq,
                    "source": B1.relative_to(root).as_posix(),
                    "sha": row.get("sequence_sha256") or seq_sha(seq),
                    "status": "recovered_from_B1_clean_candidate_universe",
                    "notes": f"clean_candidate_id={row.get('clean_candidate_id','')}; dedup_status={row.get('dedup_status','')}",
                }
    if ANCHOR not in seqs:
        raise SystemExit("STOP: verified anchor sequence cannot be recovered")
    
    seq_manifest = []
    for row in top50:
        acc = row["source_accession"]
        src = seqs.get(acc)
        seq_manifest.append({
            "accession": acc,
            "sequence_recovered": b(bool(src)),
            "sequence_length": str(len(src["seq"])) if src else "",
            "sequence_source_file": src["source"] if src else "",
            "sequence_sha256": seq_sha(src["seq"]) if src else "",
            "recovery_status": src["status"] if src else "missing_local_sequence_body",
            "notes": src["notes"] if src else "not found in searched local B0/B1 sequence stores",
        })
    
    e6_by_acc = {r["source_accession"]: r for r in read_csv(E6) if r.get("source_accession") in accessions}
    contact_by_acc = {r["source_accession"]: r for r in read_csv(CONTACT) if r.get("source_accession") in accessions}
    domains = defaultdict(set)
    for row in read_csv(DOMAIN):
        acc = row.get("accession")
        if acc not in accessions or row.get("evidence_status") != "hit":
            continue
        db = (row.get("signature_database") or "").lower()
        sig = row.get("signature_accession") or ""
        if db in {"mobidblite", "coils"} or sig in {"mobidb-lite", "Coil"}:
            continue
        token = row.get("interpro_accession") or sig
        if token:
            domains[acc].add(token)
    
    AA = "ACDEFGHIKLMNPQRSTVWY"
    HYD = set("AILMFWYV")
    CHG = set("DEKRH")
    ACI = set("DE")
    BAS = set("KRH")
    POL = set("STNQCY")
    ARO = set("FWY")
    
    
    def frac(seq, chars):
        return sum(1 for x in seq if x in chars) / len(seq) if seq else None
    
    
    def net_charge(seq):
        return (sum(1 for x in seq if x in BAS) - sum(1 for x in seq if x in ACI)) / len(seq) if seq else None
    
    
    def low_complexity(seq):
        return max(Counter(seq).values()) / len(seq) if seq else None
    
    
    def vector(seq):
        counts = Counter(seq)
        n = len(seq)
        vals = [counts.get(x, 0) / n for x in AA]
        vals += [frac(seq, HYD), frac(seq, CHG), frac(seq, ACI), frac(seq, BAS), frac(seq, POL), frac(seq, ARO),
                 frac(seq, {"C"}), frac(seq, {"P"}), low_complexity(seq), (net_charge(seq) + 1.0) / 2.0]
        return vals
    
    
    def cosine(a, c):
        dot = sum(x * y for x, y in zip(a, c))
        na = math.sqrt(sum(x * x for x in a))
        nc = math.sqrt(sum(y * y for y in c))
        return None if na == 0 or nc == 0 else max(0.0, min(1.0, dot / (na * nc)))
    
    
    def close(a, c, scale):
        return None if a is None or c is None else max(0.0, min(1.0, 1.0 - abs(a - c) / scale))
    
    
    anchor_seq = seqs[ANCHOR]["seq"]
    anchor_vec = vector(anchor_seq)
    
    
    def phys_sim(acc):
        if acc not in seqs:
            return None
        seq = seqs[acc]["seq"]
        cos = cosine(vector(seq), anchor_vec)
        closeness = mean([
            1.0 - abs(len(seq) - len(anchor_seq)) / max(len(seq), len(anchor_seq)),
            close(frac(seq, HYD), frac(anchor_seq, HYD), 1.0),
            close(frac(seq, {"C"}), frac(anchor_seq, {"C"}), 0.10),
            close(net_charge(seq), net_charge(anchor_seq), 0.50),
        ])
        if cos is None:
            return closeness
        if closeness is None:
            return cos
        return max(0.0, min(1.0, 0.70 * cos + 0.30 * closeness))
    
    
    contact_metrics = ["contact_pair_count", "close_contact_pair_count", "polar_contact_pair_count",
                       "contact_residue_count", "receptor_contact_atom_count", "min_ligand_receptor_distance_angstrom"]
    contact_vals = {m: [num(contact_by_acc.get(a, {}).get(m)) for a in accessions] for m in contact_metrics}
    broad_vals = [num(contact_by_acc.get(a, {}).get("broad_vina_score")) for a in accessions]
    anchor_top = next(r for r in top50 if r["source_accession"] == ANCHOR)
    
    
    def domain_sim(acc):
        return jaccard(domains.get(acc), domains.get(ANCHOR))
    
    
    def receptor_sim(acc):
        cur, anc = contact_by_acc.get(acc), contact_by_acc.get(ANCHOR)
        if not cur or not anc:
            return None
        return mean([minmax_sim(num(cur.get(m)), num(anc.get(m)), contact_vals[m]) for m in contact_metrics])
    
    
    def structure_sim(acc):
        cur, anc = e6_by_acc.get(acc), e6_by_acc.get(ANCHOR)
        if not cur or not anc:
            return None
        return scalar_sim(num(cur.get("E6_2_structure_reliability_score")), num(anc.get("E6_2_structure_reliability_score")))
    
    
    def target_dock_sim(acc, row):
        return mean([
            scalar_sim(num(row.get("target_context_score")), num(anchor_top.get("target_context_score"))),
            scalar_sim(num(row.get("docking_support_score")), num(anchor_top.get("docking_support_score"))),
            minmax_sim(num(contact_by_acc.get(acc, {}).get("broad_vina_score")),
                       num(contact_by_acc.get(ANCHOR, {}).get("broad_vina_score")), broad_vals),
        ])
    
    
    post, comps = [], []
    for row in top50:
        acc = row["source_accession"]
        e6 = e6_by_acc.get(acc, {})
        e64, e65, e66 = (num(e6.get("E6_4_constructability_expression_score")),
                         num(e6.get("E6_5_solubility_aggregation_chaperone_score")),
                         num(e6.get("E6_6_application_engineering_score")))
        dep, dep_norm, dep_conf = weighted([
            ("E6_4_constructability_expression_score", 0.40, e64),
            ("E6_5_solubility_aggregation_chaperone_score", 0.40, e65),
            ("E6_6_application_engineering_score", 0.20, e66),
        ])
        if not dep_norm and all(e6.get(k) == "limited_evidence" for k in ["E6_4_score_status", "E6_5_score_status", "E6_6_score_status"]):
            dep_conf = "medium"
        ds, rs, ss, ps, ts = domain_sim(acc), receptor_sim(acc), structure_sim(acc), phys_sim(acc), target_dock_sim(acc, row)
        anchor_sim, sim_norm, sim_conf = weighted([
            ("domain_architecture_similarity", 0.25, ds),
            ("receptor_binding_region_similarity", 0.20, rs),
            ("structure_similarity", 0.20, ss),
            ("physicochemical_similarity", 0.15, ps),
            ("target_context_or_docking_similarity", 0.20, ts),
        ])
        status_score = status_scores.get(acc, 0.5)
        cver = 0.70 * status_score + 0.30 * anchor_sim
        gpre = num(row["G_pre1_target_engineering_composite_score"])
        gpro = 0.60 * gpre + 0.10 * dep + 0.30 * cver
        common = {
            "G_pro1_rank": "", "accession": acc,
            "G_pre1_rank": row["G_pre1_target_engineering_composite_rank"], "G_pre1_score": fmt(gpre),
            "deployment_readiness_score": fmt(dep), "C_verified_status_score": fmt(status_score),
            "domain_architecture_similarity": fmt(ds), "receptor_binding_region_similarity": fmt(rs),
            "structure_similarity": fmt(ss), "physicochemical_similarity": fmt(ps),
            "target_context_or_docking_similarity": fmt(ts), "C_verified_anchor_similarity_score": fmt(anchor_sim),
            "C_verified_score": fmt(cver), "G_pro1_score": fmt(gpro),
            "wet_lab_status": "verified_positive" if acc in status_scores else "not_tested",
            "evidence_registry_status": "registered_direct_evidence" if acc in status_scores else "no_direct_wetlab_evidence_anchor_similarity_only",
            "sequence_recovered": b(acc in seqs), "similarity_confidence": sim_conf,
            "deployment_readiness_confidence": dep_conf,
            "forced_rank_status": row.get("forced_rank_status") or e6.get("forced_rank_status") or "not_forced",
            "hard_coded_priority_status": row.get("hard_coded_priority_status") or e6.get("hard_coded_priority_status") or "not_hard_coded",
            "manual_override_status": row.get("manual_override_status") or "none",
            "gold_standard_bonus_status": row.get("gold_standard_bonus_status") or e6.get("gold_standard_bonus_status") or "no_gold_standard_bonus",
            "candidate_policy": row.get("ordinary_candidate_policy") or e6.get("ordinary_candidate_policy") or "ordinary_candidate",
        }
        post.append(common)
        comps.append({
            **common, "source_pair_id": row.get("pair_id", ""), "source_candidate_record_id": row.get("candidate_record_id", ""),
            "E6_4_constructability_expression_score": fmt(e64), "E6_5_solubility_aggregation_chaperone_score": fmt(e65),
            "E6_6_application_engineering_score": fmt(e66),
            "deployment_score_normalized_by_available_fields": b(dep_norm),
            "score_normalized_by_available_fields": b(sim_norm),
            "domain_component_status": "computed" if ds is not None else "missing",
            "receptor_binding_region_component_status": "computed_from_contact_profile" if rs is not None else "missing",
            "structure_component_status": "computed_from_E6_2_structure_reliability_score" if ss is not None else "missing",
            "physicochemical_component_status": "computed_from_recovered_sequence" if ps is not None else "missing_sequence_or_sequence_features",
            "target_context_or_docking_component_status": "computed_from_target_context_docking_and_broad_vina" if ts is not None else "missing",
            "sequence_sha256": seqs.get(acc, {}).get("sha", ""), "sequence_source_file": seqs.get(acc, {}).get("source", ""),
        })
    
    post.sort(key=lambda r: (-float(r["G_pro1_score"]), int(r["G_pre1_rank"])))
    for i, row in enumerate(post, 1):
        row["G_pro1_rank"] = str(i)
    rank_by_acc = {r["accession"]: r["G_pro1_rank"] for r in post}
    for row in comps:
        row["G_pro1_rank"] = rank_by_acc[row["accession"]]
    comps.sort(key=lambda r: int(r["G_pro1_rank"]))
    return post, comps
