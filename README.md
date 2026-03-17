# Pea Neutron Deletion Lines

Code used to generate the VCF file for [this EVA submission](https://www.ebi.ac.uk/eva/?eva-study=PRJEB110166).

```bash
conda env create -f environment.yml -n vcf
conda activate vcf
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/pisum_sativum_gca964186695v1gb/dna/Pisum_sativum_gca964186695v1gb.JIC_Psat_v1.3.dna.toplevel.fa.gz
gunzip Pisum_sativum_gca964186695v1gb.JIC_Psat_v1.3.dna.toplevel.fa.gz
uv run main.py
uv run pytest
vcf_validator -i output.vcf
```