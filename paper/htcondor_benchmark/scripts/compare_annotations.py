#!/usr/bin/env python3
"""
compare_annotations.py - Compute sensitivity/specificity/F1 between two GFF3 files.

Following Eilbeck et al. (2009), Quantitative measures for the management and
comparison of annotated genomes, BMC Bioinformatics 10:67. This is the
canonical methodology for evaluating gene prediction tools.

Computes metrics at three levels:
1. Gene level   - whole gene match (with three stringency thresholds)
2. Exon level   - individual exon match
3. Nucleotide level - per-base coding sequence agreement

Usage:
    compare_annotations.py PREDICTED.gff3 REFERENCE.gff3 [--output report.json]

Example:
    compare_annotations.py myconote_out/annotated.gff3 SGD_R64-1-1.gff3
"""

import argparse
import json
import sys
from collections import defaultdict
from typing import Dict, List, Set, Tuple


# ─────────────────────────────────────────────────────────────────────────────
# GFF3 parsing
# ─────────────────────────────────────────────────────────────────────────────

class Feature:
    """A single GFF3 feature (gene, mRNA, CDS, exon, etc.)."""
    __slots__ = ('seqid', 'source', 'ftype', 'start', 'end', 'strand',
                 'phase', 'attributes', 'children')

    def __init__(self, seqid, source, ftype, start, end, strand, phase, attributes):
        self.seqid = seqid
        self.source = source
        self.ftype = ftype
        self.start = start
        self.end = end
        self.strand = strand
        self.phase = phase
        self.attributes = attributes
        self.children: List['Feature'] = []

    def get_id(self) -> str:
        return self.attributes.get('ID', '')

    def get_parent(self) -> str:
        return self.attributes.get('Parent', '')

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def load_rename_map(path: str) -> Dict[str, str]:
    """Load a `myconote-cli sort --rename-table` TSV and return new_id -> original_id.

    The sort step renames contigs to scaffold_N (longest first) for NCBI-clean
    output, so the predicted GFF3 uses scaffold_N seqids while the reference GFF3
    keeps the original accessions (e.g. NC_001133.9). This inverts the table
    (columns: original_id, new_id, length) so predicted seqids can be lifted back
    to the reference namespace before comparison.
    """
    mapping: Dict[str, str] = {}
    with open(path) as f:
        header = f.readline()  # original_id\tnew_id\tlength
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            original_id, new_id = parts[0], parts[1]
            mapping[new_id] = original_id
    return mapping


def parse_gff3(path: str, seqid_map: Dict[str, str] = None) -> Dict[str, List[Feature]]:
    """Parse a GFF3 file. Returns a dict of seqid -> list of gene features.

    If seqid_map is given, each column-1 seqid is translated through it (a seqid
    absent from the map is left unchanged), so a renamed prediction can be
    compared against a reference in the original namespace.
    """
    features = []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if seqid_map:
                seqid = seqid_map.get(seqid, seqid)
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            phase_val = int(phase) if phase != '.' else None

            attr_dict = {}
            for kv in attrs.split(';'):
                if '=' in kv:
                    k, v = kv.split('=', 1)
                    attr_dict[k] = v

            features.append(Feature(
                seqid=seqid, source=source, ftype=ftype,
                start=start, end=end, strand=strand,
                phase=phase_val, attributes=attr_dict,
            ))

    # Build feature ID index for parent lookups
    by_id: Dict[str, Feature] = {}
    for feat in features:
        fid = feat.get_id()
        if fid:
            by_id[fid] = feat

    # Attach children to parents
    for feat in features:
        parent_id = feat.get_parent()
        if parent_id and parent_id in by_id:
            by_id[parent_id].children.append(feat)

    # Group genes by seqid
    genes_by_seqid: Dict[str, List[Feature]] = defaultdict(list)
    for feat in features:
        if feat.ftype == 'gene':
            genes_by_seqid[feat.seqid].append(feat)

    # Sort each seqid's genes by start position
    for seqid in genes_by_seqid:
        genes_by_seqid[seqid].sort(key=lambda g: g.start)

    return dict(genes_by_seqid)


def get_cds_features(gene: Feature) -> List[Feature]:
    """Recursively collect all CDS features under a gene."""
    cds_list = []
    stack = list(gene.children)
    while stack:
        feat = stack.pop()
        if feat.ftype == 'CDS':
            cds_list.append(feat)
        stack.extend(feat.children)
    return sorted(cds_list, key=lambda c: c.start)


def get_exon_features(gene: Feature) -> List[Feature]:
    """Recursively collect all exon features under a gene."""
    exon_list = []
    stack = list(gene.children)
    while stack:
        feat = stack.pop()
        if feat.ftype == 'exon':
            exon_list.append(feat)
        stack.extend(feat.children)
    return sorted(exon_list, key=lambda e: e.start)


# ─────────────────────────────────────────────────────────────────────────────
# Overlap calculations
# ─────────────────────────────────────────────────────────────────────────────

def overlap_length(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Length of overlap between two intervals (1-based, inclusive)."""
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def features_overlap(f1: Feature, f2: Feature) -> bool:
    """Two features overlap on the same strand."""
    return (f1.strand == f2.strand
            and overlap_length(f1.start, f1.end, f2.start, f2.end) > 0)


def overlap_fraction(query: Feature, target: Feature) -> float:
    """Fraction of query feature length that overlaps target."""
    if query.strand != target.strand:
        return 0.0
    ov = overlap_length(query.start, query.end, target.start, target.end)
    return ov / query.length


def features_match_strict(f1: Feature, f2: Feature) -> bool:
    """Two features match exactly (same start, end, strand)."""
    return (f1.strand == f2.strand
            and f1.start == f2.start
            and f1.end == f2.end)


# ─────────────────────────────────────────────────────────────────────────────
# Gene-level metrics
# ─────────────────────────────────────────────────────────────────────────────

def compute_gene_metrics(predicted_genes: Dict[str, List[Feature]],
                          reference_genes: Dict[str, List[Feature]]) -> Dict:
    """
    Compute gene-level sensitivity, specificity, F1 at three stringency levels:
        loose:    any overlap on the same strand
        moderate: >=50% reciprocal overlap
        strict:   exact start AND stop coordinates match
    """
    results = {}

    for stringency in ['loose', 'moderate', 'strict']:
        tp = 0  # predicted gene matches a reference gene
        fp = 0  # predicted gene with no reference match
        fn = 0  # reference gene with no predicted match

        # Track which reference genes have been matched
        matched_refs: Dict[str, Set[int]] = defaultdict(set)

        all_seqids = set(predicted_genes) | set(reference_genes)

        for seqid in all_seqids:
            preds = predicted_genes.get(seqid, [])
            refs = reference_genes.get(seqid, [])

            for pred in preds:
                found_match = False
                for i, ref in enumerate(refs):
                    if pred.strand != ref.strand:
                        continue
                    if pred.end < ref.start:
                        break  # sorted, no further overlaps possible
                    if pred.start > ref.end:
                        continue

                    if stringency == 'loose':
                        if features_overlap(pred, ref):
                            found_match = True
                            matched_refs[seqid].add(i)
                            break
                    elif stringency == 'moderate':
                        ov = overlap_length(pred.start, pred.end, ref.start, ref.end)
                        if ov >= 0.5 * min(pred.length, ref.length):
                            found_match = True
                            matched_refs[seqid].add(i)
                            break
                    elif stringency == 'strict':
                        if features_match_strict(pred, ref):
                            found_match = True
                            matched_refs[seqid].add(i)
                            break

                if found_match:
                    tp += 1
                else:
                    fp += 1

            # Count false negatives (unmatched reference genes)
            for i, ref in enumerate(refs):
                if i not in matched_refs[seqid]:
                    fn += 1

        sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        spec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        f1 = (2 * sens * spec / (sens + spec)) if (sens + spec) > 0 else 0.0

        results[stringency] = {
            'TP': tp,
            'FP': fp,
            'FN': fn,
            'sensitivity': round(sens, 4),
            'specificity': round(spec, 4),
            'F1': round(f1, 4),
        }

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Exon-level metrics
# ─────────────────────────────────────────────────────────────────────────────

def compute_exon_metrics(predicted_genes: Dict[str, List[Feature]],
                          reference_genes: Dict[str, List[Feature]]) -> Dict:
    """Compute exon-level sensitivity and specificity."""
    pred_exons: Dict[str, Set[Tuple[int, int, str]]] = defaultdict(set)
    ref_exons: Dict[str, Set[Tuple[int, int, str]]] = defaultdict(set)

    for seqid, genes in predicted_genes.items():
        for gene in genes:
            for exon in get_exon_features(gene) or get_cds_features(gene):
                pred_exons[seqid].add((exon.start, exon.end, exon.strand))

    for seqid, genes in reference_genes.items():
        for gene in genes:
            for exon in get_exon_features(gene) or get_cds_features(gene):
                ref_exons[seqid].add((exon.start, exon.end, exon.strand))

    tp_exact = 0
    fp_exact = 0
    fn_exact = 0
    tp_overlap = 0
    fp_overlap = 0
    fn_overlap = 0

    all_seqids = set(pred_exons) | set(ref_exons)
    for seqid in all_seqids:
        preds = pred_exons.get(seqid, set())
        refs = ref_exons.get(seqid, set())

        # Exact matches
        exact_matches = preds & refs
        tp_exact += len(exact_matches)
        fp_exact += len(preds - refs)
        fn_exact += len(refs - preds)

        # Overlap-based matching (more permissive)
        ref_list = sorted(refs)
        for pred in preds:
            ps, pe, pstrand = pred
            matched = False
            for rs, re, rstrand in ref_list:
                if pstrand != rstrand:
                    continue
                if re < ps:
                    continue
                if rs > pe:
                    break
                if overlap_length(ps, pe, rs, re) > 0:
                    matched = True
                    break
            if matched:
                tp_overlap += 1
            else:
                fp_overlap += 1

        pred_list = sorted(preds)
        for ref in refs:
            rs, re, rstrand = ref
            matched = False
            for ps, pe, pstrand in pred_list:
                if rstrand != pstrand:
                    continue
                if pe < rs:
                    continue
                if ps > re:
                    break
                if overlap_length(ps, pe, rs, re) > 0:
                    matched = True
                    break
            if not matched:
                fn_overlap += 1

    def metric(tp, fp, fn):
        sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        spec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        f1 = (2 * sens * spec / (sens + spec)) if (sens + spec) > 0 else 0.0
        return {
            'TP': tp, 'FP': fp, 'FN': fn,
            'sensitivity': round(sens, 4),
            'specificity': round(spec, 4),
            'F1': round(f1, 4),
        }

    return {
        'exact_match': metric(tp_exact, fp_exact, fn_exact),
        'overlap_match': metric(tp_overlap, fp_overlap, fn_overlap),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Nucleotide-level metrics
# ─────────────────────────────────────────────────────────────────────────────

def compute_nucleotide_metrics(predicted_genes: Dict[str, List[Feature]],
                                reference_genes: Dict[str, List[Feature]]) -> Dict:
    """
    Per-nucleotide coding region accuracy.
    For each seqid, build a set of coding positions from each annotation,
    then compute the intersection.
    """
    def coding_positions(genes_by_seqid):
        positions: Dict[str, Set[int]] = defaultdict(set)
        for seqid, genes in genes_by_seqid.items():
            for gene in genes:
                for cds in get_cds_features(gene):
                    for pos in range(cds.start, cds.end + 1):
                        positions[seqid].add(pos)
        return positions

    pred_coding = coding_positions(predicted_genes)
    ref_coding = coding_positions(reference_genes)

    tp = 0
    fp = 0
    fn = 0
    all_seqids = set(pred_coding) | set(ref_coding)
    for seqid in all_seqids:
        p = pred_coding.get(seqid, set())
        r = ref_coding.get(seqid, set())
        tp += len(p & r)
        fp += len(p - r)
        fn += len(r - p)

    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = (2 * sens * spec / (sens + spec)) if (sens + spec) > 0 else 0.0

    return {
        'TP_nt': tp,
        'FP_nt': fp,
        'FN_nt': fn,
        'sensitivity': round(sens, 4),
        'specificity': round(spec, 4),
        'F1': round(f1, 4),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Summary statistics
# ─────────────────────────────────────────────────────────────────────────────

def count_features(genes_by_seqid: Dict[str, List[Feature]]) -> Dict[str, int]:
    """Count basic features in an annotation."""
    n_genes = 0
    n_exons = 0
    n_cds = 0
    total_gene_length = 0
    n_single_exon = 0

    for seqid, genes in genes_by_seqid.items():
        for gene in genes:
            n_genes += 1
            total_gene_length += gene.length
            exons = get_exon_features(gene) or get_cds_features(gene)
            n_exons += len(exons)
            n_cds += len(get_cds_features(gene))
            if len(exons) == 1:
                n_single_exon += 1

    return {
        'n_genes': n_genes,
        'n_exons': n_exons,
        'n_cds': n_cds,
        'mean_gene_length': round(total_gene_length / n_genes, 1) if n_genes else 0,
        'n_single_exon_genes': n_single_exon,
        'mean_exons_per_gene': round(n_exons / n_genes, 2) if n_genes else 0,
    }


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Compare two GFF3 annotations and compute accuracy metrics',
        epilog='Following Eilbeck et al. (2009) BMC Bioinformatics 10:67.',
    )
    parser.add_argument('predicted', help='Predicted annotation GFF3')
    parser.add_argument('reference', help='Reference (gold standard) annotation GFF3')
    parser.add_argument('--output', '-o', help='Output JSON file (default: stdout)')
    parser.add_argument('--label', help='Label for this comparison (e.g., "myconote_vs_SGD")')
    parser.add_argument('--rename-table',
                        help='sort --rename-table TSV (original_id, new_id, length); '
                             'lifts predicted seqids back to the reference namespace')
    args = parser.parse_args()

    seqid_map = None
    if args.rename_table:
        seqid_map = load_rename_map(args.rename_table)
        print(f'Loaded rename table: {len(seqid_map)} contigs '
              f'(predicted seqids lifted to reference namespace)', file=sys.stderr)

    print(f'Parsing predicted annotation: {args.predicted}', file=sys.stderr)
    pred = parse_gff3(args.predicted, seqid_map=seqid_map)
    pred_summary = count_features(pred)

    print(f'Parsing reference annotation: {args.reference}', file=sys.stderr)
    ref = parse_gff3(args.reference)
    ref_summary = count_features(ref)

    print(f'Predicted: {pred_summary["n_genes"]:,} genes', file=sys.stderr)
    print(f'Reference: {ref_summary["n_genes"]:,} genes', file=sys.stderr)
    print('Computing metrics...', file=sys.stderr)

    gene_metrics = compute_gene_metrics(pred, ref)
    exon_metrics = compute_exon_metrics(pred, ref)
    nt_metrics = compute_nucleotide_metrics(pred, ref)

    report = {
        'label': args.label or '',
        'predicted_file': args.predicted,
        'reference_file': args.reference,
        'predicted_summary': pred_summary,
        'reference_summary': ref_summary,
        'gene_level': gene_metrics,
        'exon_level': exon_metrics,
        'nucleotide_level': nt_metrics,
    }

    # Pretty-print summary to stderr
    print('', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print('ACCURACY METRICS', file=sys.stderr)
    print('=' * 70, file=sys.stderr)
    print('', file=sys.stderr)
    print('Gene-level:', file=sys.stderr)
    for s in ['loose', 'moderate', 'strict']:
        m = gene_metrics[s]
        print(f"  {s:10s}  Sens={m['sensitivity']:.3f}  "
              f"Spec={m['specificity']:.3f}  F1={m['F1']:.3f}  "
              f"(TP={m['TP']:,} FP={m['FP']:,} FN={m['FN']:,})",
              file=sys.stderr)
    print('', file=sys.stderr)
    print('Exon-level:', file=sys.stderr)
    for s in ['exact_match', 'overlap_match']:
        m = exon_metrics[s]
        print(f"  {s:14s}  Sens={m['sensitivity']:.3f}  "
              f"Spec={m['specificity']:.3f}  F1={m['F1']:.3f}",
              file=sys.stderr)
    print('', file=sys.stderr)
    print('Nucleotide-level:', file=sys.stderr)
    print(f"  Sens={nt_metrics['sensitivity']:.3f}  "
          f"Spec={nt_metrics['specificity']:.3f}  F1={nt_metrics['F1']:.3f}",
          file=sys.stderr)
    print('=' * 70, file=sys.stderr)

    output_text = json.dumps(report, indent=2)
    if args.output:
        with open(args.output, 'w') as f:
            f.write(output_text)
        print(f'Wrote: {args.output}', file=sys.stderr)
    else:
        print(output_text)


if __name__ == '__main__':
    main()
