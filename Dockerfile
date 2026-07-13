FROM debian@sha256:d07d1b51c39f51188e60be9b64e6bf769fa94e187f092bc32b91305cfa34ba5a
LABEL org.opencontainers.image.source=https://github.com/Informatics-John-Innes-Centre/pea-neutron-deletion-lines

WORKDIR /app
COPY . .

RUN apt-get update && \
    apt-get install -y curl=8.14.1-2+deb13u4 && \
    curl -L -o miniforge.sh \
      https://github.com/conda-forge/miniforge/releases/download/26.3.2-3/Miniforge3-26.3.2-3-Linux-x86_64.sh && \
    bash miniforge.sh -b -p /opt/miniforge3 && \
    rm miniforge.sh


ENV PATH="/opt/miniforge3/bin:/opt/miniforge3/condabin:${PATH}"

RUN rm /opt/miniforge3/.condarc && \
conda config --set channel_alias "https://repo.prefix.dev" && \
conda env create -f environment.yml -n vcf

RUN conda run -n vcf uv sync --frozen

RUN conda run -n vcf wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/pisum_sativum_gca964186695v1gb/dna/Pisum_sativum_gca964186695v1gb.JIC_Psat_v1.3.dna.toplevel.fa.gz && \
    conda run -n vcf gunzip Pisum_sativum_gca964186695v1gb.JIC_Psat_v1.3.dna.toplevel.fa.gz

ENTRYPOINT ["conda", "run", "-n", "vcf", "uv", "run", "main.py"]
