#!/usr/bin/env python3
"""
extract_protein_metrics.py
Compute functional annotation quality metrics from an annotated GFF3.

Measures:
  - fraction of genes with a non-empty product description
  - fraction with a Swiss-Prot hit (from `Dbxref=SwissProt:` or `db_xref=`)
  - fraction with a Pfam domain (from `Dbxref=Pfam:` or `Ontology_term=Pfam:`)
  - fraction with at least one GO term (from `Ontology_term=GO:`)
  - fraction with an EC number (from `Dbxref=EC:`)
  - mean/median product description length

Also validates the companion protein FASTA (if provided) by counting:
  - proteins with internal stop codons (invalid)
  - proteins shorter than a minimum length (default 30 aa)

Usage:
    python3 extract_protein_metrics.py ANNOTATED.gff3 \
        [--proteins proteins.fa] \
        [--min-protein-len 30] \
        --output metrics.json \
        --label LABEL
"""

import argparse
import json
import re
import statistics
import sys
from pathlib import Path


GENE_TYPES = {'gene', 'pseudogene', 'ncRNA_gene'}


def parse_attributes(attr_str: str) -> dict:
    """Parse GFF3 column 9 into a dict. Multi-value keys are kept as lists."""
    out = {}
    for kv in attr_str.strip().split(';'):
        if '=' not in kv:
            continue
        k, v = kv.split('=', 1)
        k = k.strip()
        v = v.strip()
        if not k:
            continue
        if k in out:
            if isinstance(out[k], list):
                out[k].append(v)
            else:
                out[k] = [out[k], v]
        else:
            out[k] = v
    return out


def attr_values(attrs: dict, key: str) -> list:
    """Return all comma-separated values for a GFF3 attribute key."""
    v = attrs.get(key)
    if v is None:
        return []
    if isinstance(v, list):
        flat = []
        for item in v:
            flat.extend([x.strip() for x in item.split(',') if x.strip()])
        return flat
    return [x.strip() for x in v.split(',') if x.strip()]


def iterate_gene_attrs(gff_path: Path):
    """Yield (gene_id, attrs_dict) for every gene feature in a GFF3."""
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype not in GENE_TYPES:
                continue
            attrs = parse_attributes(parts[8])
            gene_id = attrs.get('ID') or attrs.get('Name') or '?'
            yield gene_id, attrs


def has_hit(attrs: dict, needle: str) -> bool:
    """True if any Dbxref/db_xref/Ontology_term value starts with `needle:`."""
    needle_lc = needle.lower() + ':'
    for key in ('Dbxref', 'db_xref', 'Ontology_term'):
        for v in attr_values(attrs, key):
            if v.lower().startswith(needle_lc):
                return True
    return False


def parse_fasta_proteins(fasta_path: Path):
    """Yield (header, sequence) tuples from a protein FASTA."""
    header = None
    chunks = []
    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, ''.join(chunks)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('gff', help='Annotated GFF3 to evaluate')
    ap.add_argument('--proteins', help='Optional protein FASTA for validity checks')
    ap.add_argument('--min-protein-len', type=int, default=30,
                    help='Minimum acceptable protein length (default 30)')
    ap.add_argument('--output', required=True, help='Output JSON path')
    ap.add_argument('--label', default='', help='Label to include in JSON')
    args = ap.parse_args()

    gff_path = Path(args.gff)
    if not gff_path.exists():
        print(f'ERROR: GFF3 not found: {gff_path}', file=sys.stderr)
        sys.exit(1)

    n_genes = 0
    n_with_product = 0
    n_with_swissprot = 0
    n_with_pfam = 0
    n_with_go = 0
    n_with_ec = 0
    product_lengths = []
    generic_patterns = re.compile(
        r'^(hypothetical protein|unknown|uncharacterized|predicted protein)\s*$',
        re.IGNORECASE,
    )
    n_generic_products = 0

    for _gene_id, attrs in iterate_gene_attrs(gff_path):
        n_genes += 1

        product = attrs.get('product') or attrs.get('Note') or ''
        if isinstance(product, list):
            product = product[0] if product else ''
        product = product.strip()
        if product:
            n_with_product += 1
            product_lengths.append(len(product))
            if generic_patterns.match(product):
                n_generic_products += 1

        if has_hit(attrs, 'SwissProt') or has_hit(attrs, 'UniProt'):
            n_with_swissprot += 1
        if has_hit(attrs, 'Pfam'):
            n_with_pfam += 1
        if has_hit(attrs, 'GO'):
            n_with_go += 1
        if has_hit(attrs, 'EC'):
            n_with_ec += 1

    def frac(x):
        return round(x / n_genes, 4) if n_genes else 0.0

    result = {
        'label': args.label,
        'gff3': str(gff_path),
        'n_genes': n_genes,
        'product': {
            'n_with_product': n_with_product,
            'fraction_with_product': frac(n_with_product),
            'n_generic_products': n_generic_products,
            'fraction_informative': frac(n_with_product - n_generic_products),
            'mean_product_length': round(statistics.mean(product_lengths), 1)
                if product_lengths else 0,
            'median_product_length': int(statistics.median(product_lengths))
                if product_lengths else 0,
        },
        'swissprot': {
            'n_with_hit': n_with_swissprot,
            'fraction': frac(n_with_swissprot),
        },
        'pfam': {
            'n_with_domain': n_with_pfam,
            'fraction': frac(n_with_pfam),
        },
        'go': {
            'n_with_term': n_with_go,
            'fraction': frac(n_with_go),
        },
        'ec': {
            'n_with_number': n_with_ec,
            'fraction': frac(n_with_ec),
        },
    }

    if args.proteins:
        prot_path = Path(args.proteins)
        if not prot_path.exists():
            print(f'Warning: protein FASTA not found: {prot_path}', file=sys.stderr)
        else:
            n_prot = 0
            n_internal_stop = 0
            n_too_short = 0
            lens = []
            for _hdr, seq in parse_fasta_proteins(prot_path):
                n_prot += 1
                seq_stripped = seq.rstrip('*')
                if '*' in seq_stripped:
                    n_internal_stop += 1
                if len(seq_stripped) < args.min_protein_len:
                    n_too_short += 1
                lens.append(len(seq_stripped))
            result['proteins'] = {
                'n_proteins': n_prot,
                'n_with_internal_stop': n_internal_stop,
                'fraction_invalid_internal_stop': round(n_internal_stop / n_prot, 4)
                    if n_prot else 0.0,
                'n_shorter_than_min': n_too_short,
                'min_length_threshold': args.min_protein_len,
                'mean_length_aa': round(statistics.mean(lens), 1) if lens else 0,
                'median_length_aa': int(statistics.median(lens)) if lens else 0,
            }

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(f'Wrote: {out_path}')


if __name__ == '__main__':
    main()
