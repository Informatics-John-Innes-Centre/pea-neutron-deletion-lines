# Pea Neutron Deletion Lines

Code used to generate the VCF file for [this EVA submission](https://www.ebi.ac.uk/eva/?eva-study=PRJEB110166).

## Manual Installation

```bash
conda env create -f environment.yml -n vcf
conda activate vcf
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/pisum_sativum_gca964186695v1gb/dna/Pisum_sativum_gca964186695v1gb.JIC_Psat_v1.3.dna.toplevel.fa.gz
gunzip Pisum_sativum_gca964186695v1gb.JIC_Psat_v1.3.dna.toplevel.fa.gz
uv run main.py
uv run pytest
vcf_validator -i output.vcf
```

## Docker

```bash
sudo docker build -t pea-neutron-deletion-lines .
sudo docker run --rm -v $(pwd):/output pea-neutron-deletion-lines -o /output/output.vcf
```

## Singularity

```bash
sudo docker build -t pea-neutron-deletion-lines .
sudo singularity build pea-neutron-deletion-lines.sif docker-daemon://pea-neutron-deletion-lines:latest
```