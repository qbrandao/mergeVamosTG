# mergeVamosTG

Outil en Rust permettant de fusionner et d'annoter les STRs provenant de VaMoS et Tandem-Genotypes.
Le processus d'annotation génomique (Gènes, Exons, Introns, UTRs) est parallélisé avec Rayon.

## Prérequis
- Rust (cargo)

## Utilisation
```bash
cargo run --release <vamos_hap1.vcf> <vamos_hap2.vcf> <tg_records.tsv> <output_file>
