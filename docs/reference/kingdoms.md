# Kingdom Defaults

MycoNote-CLI ships with sensible per-kingdom defaults so you don't need to
specify every parameter manually. Pass `--kingdom` to a pipeline command to
activate these defaults. Fungi is the primary, validated target; the other
kingdoms are an experimental scaffold (see the Supported Kingdoms note).

## Augustus species models

| Kingdom | Default species | Notes |
|---------|----------------|-------|
| Fungi | `saccharomyces_cerevisiae_S288C` | Ascomycota; see `--clade` below for a clade-aware default |
| Plant | `arabidopsis` | *Arabidopsis thaliana* |
| Animal | `human` | *Homo sapiens* |
| Insect | `fly` | *Drosophila melanogaster* |
| Protist | `toxoplasma` | *Toxoplasma gondii* |

Override with `--species <model>` to use any installed Augustus species.

!!! tip "Clade-aware fungal default (`--clade`)"
    For fungi, `saccharomyces_cerevisiae_S288C` is a Saccharomycotina model and
    a poor default off that clade. Instead of hard-coding it, pass
    `predict --clade saccharomycotina|pezizomycotina|basidiomycota` to pick a
    clade-appropriate model (Pezizomycotina → `aspergillus_nidulans`,
    Basidiomycota → `cryptococcus_neoformans_neoformans_JEC21`). A fungal run
    with neither `--species` nor `--clade` still uses S288C but prints a warning.

## RepeatMasker clades

These match the per-kingdom clade MycoNote passes to repeat masking:

| Kingdom | Clade |
|---------|-------|
| Fungi | Fungi |
| Plant | Viridiplantae |
| Animal | Metazoa |
| Insect | Insecta |
| Protist | Protista |

## BUSCO lineages (suggested)

Only `fungi_odb10` is wired as a default; for other kingdoms set the lineage
explicitly with `--busco-lineage`. Conventional choices:

| Kingdom | Suggested lineage |
|---------|-------------------|
| Fungi | `fungi_odb10` (default) |
| Plant | `embryophyta_odb10` |
| Animal | `metazoa_odb10` |
| Insect | `insecta_odb10` |
| Protist | `eukaryota_odb10` |

Override with `--busco-lineage` on `myconote-cli annotate`.
