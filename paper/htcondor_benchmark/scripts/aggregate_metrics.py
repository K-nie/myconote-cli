#!/usr/bin/env python3
"""
aggregate_metrics.py
Combine all per-job metrics JSON files into manuscript-ready tables.

Reads:
  results/{tool}/{genome_id}/rep{N}/metrics.json
  results/{tool}/{genome_id}/rep{N}/performance.json

Writes:
  results/aggregated/per_run_metrics.tsv      - one row per (tool, genome, rep)
  results/aggregated/per_tool_summary.tsv     - mean/SD per tool
  results/aggregated/per_genome_summary.tsv   - mean/SD per genome
  results/aggregated/comparison_table.tsv     - manuscript-ready table
  results/aggregated/statistical_tests.tsv    - paired t-tests, MycoNote vs others

Usage:
    python3 aggregate_metrics.py [results_dir]
"""

import json
import os
import sys
import statistics
from collections import defaultdict
from pathlib import Path


ALL_TOOLS = ['myconote', 'funannotate', 'maker', 'braker']
TOOLS = list(ALL_TOOLS)  # overwritten in main() to just the tools that ran
METRIC_FIELDS = [
    'gene_loose_sens', 'gene_loose_spec', 'gene_loose_F1',
    'gene_strict_sens', 'gene_strict_spec', 'gene_strict_F1',
    'exon_overlap_sens', 'exon_overlap_spec', 'exon_overlap_F1',
    'nt_sens', 'nt_spec', 'nt_F1',
    'wall_seconds', 'peak_rss_mb',
    'n_predicted_genes', 'n_reference_genes',
]


def flatten_metrics(metrics_json: dict, performance_json: dict) -> dict:
    """Convert nested metrics JSON into a flat dict for tabular output."""
    out = {}

    gene = metrics_json.get('gene_level', {})
    if 'loose' in gene:
        out['gene_loose_sens'] = gene['loose'].get('sensitivity', 0)
        out['gene_loose_spec'] = gene['loose'].get('specificity', 0)
        out['gene_loose_F1'] = gene['loose'].get('F1', 0)
    if 'strict' in gene:
        out['gene_strict_sens'] = gene['strict'].get('sensitivity', 0)
        out['gene_strict_spec'] = gene['strict'].get('specificity', 0)
        out['gene_strict_F1'] = gene['strict'].get('F1', 0)

    exon = metrics_json.get('exon_level', {}).get('overlap_match', {})
    out['exon_overlap_sens'] = exon.get('sensitivity', 0)
    out['exon_overlap_spec'] = exon.get('specificity', 0)
    out['exon_overlap_F1'] = exon.get('F1', 0)

    nt = metrics_json.get('nucleotide_level', {})
    out['nt_sens'] = nt.get('sensitivity', 0)
    out['nt_spec'] = nt.get('specificity', 0)
    out['nt_F1'] = nt.get('F1', 0)

    out['n_predicted_genes'] = metrics_json.get('predicted_summary', {}).get('n_genes', 0)
    out['n_reference_genes'] = metrics_json.get('reference_summary', {}).get('n_genes', 0)

    out['wall_seconds'] = performance_json.get('total_wall_seconds', 0)
    out['peak_rss_mb'] = performance_json.get('peak_rss_mb', 0)

    return out


def load_run(tool: str, genome_id: str, rep: int, results_dir: Path) -> dict:
    """Load metrics for a single run, returning empty dict if missing."""
    run_dir = results_dir / tool / genome_id / f'rep{rep}'
    metrics_path = run_dir / 'metrics.json'
    perf_path = run_dir / 'performance.json'

    if not metrics_path.exists():
        return {}

    try:
        with open(metrics_path) as f:
            metrics = json.load(f)
    except Exception as e:
        print(f'  Warning: failed to load {metrics_path}: {e}', file=sys.stderr)
        return {}

    perf = {}
    if perf_path.exists():
        try:
            with open(perf_path) as f:
                perf = json.load(f)
        except Exception as e:
            print(f'  Warning: failed to load {perf_path}: {e}', file=sys.stderr)

    flat = flatten_metrics(metrics, perf)
    flat['tool'] = tool
    flat['genome'] = genome_id
    flat['replicate'] = rep
    return flat


def load_genomes(config_path: Path) -> list:
    """Read the list of genome IDs from the config TSV."""
    genomes = []
    with open(config_path) as f:
        next(f)  # skip header
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue  # skip blank and comment lines
            parts = line.split('\t')
            if parts and parts[0]:
                genomes.append(parts[0])
    return genomes


def write_per_run(rows: list, out_path: Path):
    """Write the long-form per-run metrics table."""
    if not rows:
        print('No rows to write', file=sys.stderr)
        return

    fields = ['tool', 'genome', 'replicate'] + METRIC_FIELDS
    with open(out_path, 'w') as f:
        f.write('\t'.join(fields) + '\n')
        for row in rows:
            f.write('\t'.join(str(row.get(k, '')) for k in fields) + '\n')
    print(f'Wrote: {out_path}')


def write_per_tool_summary(rows: list, out_path: Path):
    """Write per-tool summary with mean/SD across all runs."""
    by_tool = defaultdict(list)
    for row in rows:
        if row.get('tool'):
            by_tool[row['tool']].append(row)

    with open(out_path, 'w') as f:
        f.write('tool\tn_runs')
        for field in METRIC_FIELDS:
            f.write(f'\t{field}_mean\t{field}_sd')
        f.write('\n')

        for tool in TOOLS:
            tool_rows = by_tool.get(tool, [])
            f.write(f'{tool}\t{len(tool_rows)}')
            for field in METRIC_FIELDS:
                vals = [r[field] for r in tool_rows
                        if isinstance(r.get(field), (int, float)) and r[field] > 0]
                if vals:
                    mean = round(statistics.mean(vals), 4)
                    sd = round(statistics.stdev(vals), 4) if len(vals) > 1 else 0.0
                else:
                    mean = sd = 0.0
                f.write(f'\t{mean}\t{sd}')
            f.write('\n')
    print(f'Wrote: {out_path}')


def write_comparison_table(rows: list, out_path: Path):
    """Write a manuscript-ready comparison table."""
    by_tool = defaultdict(list)
    for row in rows:
        if row.get('tool'):
            by_tool[row['tool']].append(row)

    metrics_to_show = [
        ('gene_loose_F1', 'Gene F1 (loose)'),
        ('gene_strict_F1', 'Gene F1 (strict)'),
        ('exon_overlap_F1', 'Exon F1'),
        ('nt_F1', 'Nucleotide F1'),
        ('wall_seconds', 'Runtime (s)'),
        ('peak_rss_mb', 'Peak memory (MB)'),
    ]

    with open(out_path, 'w') as f:
        # Header
        f.write('Metric')
        for tool in TOOLS:
            f.write(f'\t{tool}')
        f.write('\n')

        # Rows
        for field, label in metrics_to_show:
            f.write(label)
            for tool in TOOLS:
                tool_rows = by_tool.get(tool, [])
                vals = [r[field] for r in tool_rows
                        if isinstance(r.get(field), (int, float)) and r[field] > 0]
                if vals:
                    mean = statistics.mean(vals)
                    sd = statistics.stdev(vals) if len(vals) > 1 else 0
                    f.write(f'\t{mean:.3f} ± {sd:.3f}')
                else:
                    f.write('\tNA')
            f.write('\n')

    print(f'Wrote: {out_path}')


def write_statistical_tests(rows: list, out_path: Path):
    """Paired t-tests comparing MycoNote-CLI vs each other tool."""
    other_tools = [t for t in TOOLS if t != 'myconote']
    if not other_tools:
        print('Only myconote results present — skipping statistical tests '
              '(nothing to compare against)', file=sys.stderr)
        return
    try:
        from scipy import stats
    except ImportError:
        print('scipy not available — skipping statistical tests', file=sys.stderr)
        return

    by_key = defaultdict(dict)
    for row in rows:
        key = (row.get('genome'), row.get('replicate'))
        tool = row.get('tool')
        by_key[key][tool] = row

    test_metrics = [
        'gene_loose_F1', 'gene_strict_F1', 'exon_overlap_F1', 'nt_F1',
        'wall_seconds', 'peak_rss_mb',
    ]

    with open(out_path, 'w') as f:
        f.write('comparison\tmetric\tn_pairs\tmyconote_mean\tother_mean\t'
                'mean_difference\tt_stat\tp_value\teffect_size_d\n')

        for other in other_tools:
            for metric in test_metrics:
                myco_vals = []
                other_vals = []
                for key, tools in by_key.items():
                    if 'myconote' in tools and other in tools:
                        m = tools['myconote'].get(metric)
                        o = tools[other].get(metric)
                        if isinstance(m, (int, float)) and isinstance(o, (int, float)):
                            if m > 0 and o > 0:
                                myco_vals.append(m)
                                other_vals.append(o)

                n = len(myco_vals)
                if n < 2:
                    f.write(f'myconote_vs_{other}\t{metric}\t{n}\tNA\tNA\tNA\tNA\tNA\tNA\n')
                    continue

                t_stat, p_value = stats.ttest_rel(myco_vals, other_vals)
                myco_mean = statistics.mean(myco_vals)
                other_mean = statistics.mean(other_vals)
                diff = myco_mean - other_mean

                # Cohen's d for paired samples
                pooled_sd = statistics.stdev(
                    [m - o for m, o in zip(myco_vals, other_vals)]
                )
                cohens_d = diff / pooled_sd if pooled_sd > 0 else 0

                f.write(f'myconote_vs_{other}\t{metric}\t{n}\t'
                        f'{myco_mean:.4f}\t{other_mean:.4f}\t{diff:.4f}\t'
                        f'{t_stat:.3f}\t{p_value:.4f}\t{cohens_d:.3f}\n')

    print(f'Wrote: {out_path}')


def main():
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[1])
    else:
        script_dir = Path(__file__).parent
        results_dir = script_dir.parent / 'results'

    if not results_dir.exists():
        print(f'ERROR: Results directory not found: {results_dir}', file=sys.stderr)
        sys.exit(1)

    config_path = results_dir.parent / 'configs' / 'genomes.tsv'
    if not config_path.exists():
        print(f'ERROR: Config file not found: {config_path}', file=sys.stderr)
        sys.exit(1)

    genomes = load_genomes(config_path)
    print(f'Found {len(genomes)} genomes in config')
    print(f'Loading metrics from {results_dir}...')

    # Only aggregate tools whose result directory exists — supports
    # partial runs (e.g. --myconote-only) without emitting a flood of
    # NA rows for tools that were never executed.
    tools = [t for t in ALL_TOOLS if (results_dir / t).is_dir()]
    if not tools:
        print(f'ERROR: no tool result directories found under {results_dir}',
              file=sys.stderr)
        sys.exit(1)
    print(f'Tools with results: {", ".join(tools)}')

    # Expose the active tool set to the writer helpers via a module global
    # so they can use it for ordering columns.
    global TOOLS
    TOOLS = tools

    all_rows = []
    for tool in tools:
        for genome in genomes:
            for rep in [1, 2, 3]:
                row = load_run(tool, genome, rep, results_dir)
                if row:
                    all_rows.append(row)

    print(f'Loaded {len(all_rows)} runs (out of {len(tools) * len(genomes) * 3} possible)')

    # Create aggregated output directory
    agg_dir = results_dir / 'aggregated'
    agg_dir.mkdir(exist_ok=True)

    write_per_run(all_rows, agg_dir / 'per_run_metrics.tsv')
    write_per_tool_summary(all_rows, agg_dir / 'per_tool_summary.tsv')
    write_comparison_table(all_rows, agg_dir / 'comparison_table.tsv')
    write_statistical_tests(all_rows, agg_dir / 'statistical_tests.tsv')

    print('')
    print(f'Aggregation complete. Outputs in: {agg_dir}')


if __name__ == '__main__':
    main()
