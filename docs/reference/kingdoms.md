# Kingdom Defaults

Myconote_CLI ships with sensible per-kingdom defaults so you don't need to specify every parameter manually. Pass `--kingdom` to any pipeline command to activate these defaults.

## Augustus species models

| Kingdom | Default species | Notes |
|---------|----------------|-------|
| Fungi | `saccharomyces_cerevisiae_S288C` | Ascomycota; good general-purpose fungal model |
| Plant | `arabidopsis` | *Arabidopsis thaliana* |
| Animal | `human` | *Homo sapiens* |
| Protist | `tetrahymena` | *Tetrahymena thermophila* |

Override with `--species <model>` to use any installed Augustus species.

## RepeatMasker libraries

| Kingdom | Library |
|---------|---------|
| Fungi | Fungi clade from Dfam |
| Plant | Viridiplantae clade |
| Animal | Metazoa clade |
| Protist | Full Dfam library |

## BUSCO lineages

| Kingdom | Default lineage |
|---------|----------------|
| Fungi | `fungi_odb10` |
| Plant | `embryophyta_odb10` |
| Animal | `metazoa_odb10` |
| Protist | `eukaryota_odb10` |

Override with `--busco-lineage` on `myconote stats` or `myconote annotate`.
