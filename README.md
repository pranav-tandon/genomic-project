# genomic-project

# Variant Calling Workflow – README

## Overview

This repository contains a fully reproducible **Snakemake** workflow that processes raw whole‑exome sequencing (WES) reads (paired‑end FASTQ) to high‑confidence genetic variants (VCF/BCF). The pipeline integrates:

| Stage               | Tool                          | Purpose                                  |
| ------------------- | ----------------------------- | ---------------------------------------- |
| **QC**              | FastQC                        | Assess raw read quality                  |
| **Alignment**       | Bowtie2                       | Map reads to the *hg38* reference genome |
| **Post‑process**    | SAMtools                      | Convert, sort & index alignments         |
| **Variant Calling** | BCFtools (+ SAMtools mpileup) | Detect SNPs & indels, produce VCF/BCF    |

All software dependencies are isolated in a dedicated **Conda** environment (`genomics_env`) to guarantee identical results across machines.

---

## Quick Start

```bash
# 1 Clone repository and enter project dir
$ git clone <repo‑url>
$ cd variant‑pipeline

# 2 Create Conda env (+ Mamba recommended for speed)
$ conda env create -n genomics_env -f env.yml   # or: mamba env create ...
$ conda activate genomics_env

# 3 Download / link raw FASTQs & reference index
$ mkdir data
$ cp /path/to/subset_SRR099957_*fastq.gz data/
$ cp -r /path/to/hg38_bowtie2_index reference/
$ cp /path/to/hg38.fa reference/

# 4 Dry‑run (shows DAG & commands)
$ snakemake -npr --use-conda

# 5 Execute workflow (single core, 4 GB RAM)
$ snakemake --cores 1 --resources mem_mb=4000 --use-conda
```

> **Tip:** Increase `--cores` to parallelise multi‑sample runs or use a Snakemake cluster profile for HPC/cloud execution.

---

## Input Files

| Path                               | Description            |
| ---------------------------------- | ---------------------- |
| `data/subset_SRR099957_1.fastq.gz` | Forward reads          |
| `data/subset_SRR099957_2.fastq.gz` | Reverse reads          |
| `reference/hg38.*.bt2`             | Bowtie2 index files    |
| `reference/hg38.fa`                | Reference genome FASTA |

All other paths are generated automatically by Snakemake as defined in `config/config.yaml`.

---

## Output Artifacts

After successful execution, you will find:

```
results/
├── fastqc/
│   ├── subset_SRR099957_1_fastqc.html
│   └── subset_SRR099957_2_fastqc.html
├── alignments/
│   ├── subset_SRR099957.sorted.bam
│   └── subset_SRR099957.sorted.bam.bai
└── variants/
    ├── subset_SRR099957.bcf
    └── subset_SRR099957.vcf
```

*FastQC* reports can be opened in any browser; BAM/VCF files can be inspected with IGV or command‑line tools like `samtools` and `bcftools`.

---

## Configuration

Edit `config/config.yaml` to adjust:

* **samples:** FASTQ filenames or pattern for multi‑sample runs.
* **reference:** path to Bowtie2 index & FASTA.
* **threads / memory:** per‑rule resource usage.

---

## Extending the Pipeline

The modular Snakefile makes it easy to:

1. **Add read trimming** – insert a rule using *fastp* or *Trimmomatic* before alignment.
2. **Switch aligner** – replace `align_reads_bowtie2` with BWA‑MEM2 or HISAT2.
3. **Add variant filtering / annotation** – append rules invoking *bcftools filter*, *ANNOVAR*, or *VEP*.
4. **Run joint calling** – create a rule that merges or jointly calls across multiple samples.

---

## Troubleshooting

| Issue                            | Fix                                                                   |
| -------------------------------- | --------------------------------------------------------------------- |
| Conda solve is slow or conflicts | Use *Mamba* (`conda install mamba -n base -c conda-forge`).           |
| `bowtie2` cannot find index      | Ensure paths in `config.yaml` point to all six `.bt2` files.          |
| Snakemake stops on missing files | Run with `-npr` to view expected filenames, then check `config.yaml`. |

---

## Citation

If you use this workflow, please cite:

* Köster & Rahmann, *Snakemake* (2012)
* Grüning et al., *Bioconda* (2018)
* Langmead & Salzberg, *Bowtie 2* (2012)
* Li et al., *SAMtools* (2009)
* Danecek et al., *BCFtools* (2021)
* Babraham Bioinformatics, *FastQC* (2019)

---

## License

MIT License – see `LICENSE` file.

---

## Author & Acknowledgements

Developed and maintained by **Pranav Dhruv Tandon** (UBC).
Based on community best practices in genomics data analysis.

Special thanks to the open‑source developers of FastQC, Bowtie2, SAMtools/BCFtools, Snakemake, and the Bioconda community.
