![nf-core/exogap](./docs/images/exogap_logo_light.png#gh-light-mode-only)
![nf-core/exogap](./docs/images/exogap_logo_dark.png#gh-dark-mode-only)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

EXOGAP (EXotic Organism Genome Annotation Pipeline) is a Nextflow DSL2 workflow for comprehensive annotation of genomes from non-model species. It integrates tools for repetitive element identification, protein-coding gene prediction and non-coding gene prediction. The pipeline can annotate multiple related genomes in a single run and formats outputs for deposition in public databases.

Main modules:
- Repetitive element annotation: RepeatModeler, RepeatMasker (uses FamDB).
- Protein-coding gene annotation: MAKER2, AUGUSTUS, SNAP, BLAST, InterProScan.
- Non-coding gene annotation: tRNAscan-SE, RNAmmer, Barrnap, Infernal, SnoScan.

## Installation

Prerequisites:
- Nextflow (recommended >= 24.04)
- Singularity / Apptainer

Clone the repository:
```bash
git clone https://github.com/dorinemerlat/exogap.git
cd exogap
```

## Usage

### Samplesheet CSV file

Prepare a samplesheet CSV describing genomes to annotate.
## Pipeline inputs & options

- **Samplesheet (required):** a CSV describing genomes. Required columns: `name`, `taxid`, `fasta`. Optional columns supported: `RNASeq-dir` (path to directory with FASTQ files), `SRA` (SRA accessions), `repeats_set` (FASTA library).
- **Main config:** `nextflow.config` (or `example/nextflow_example.config`) controls runtime options. Key params:
	- `input`: path to samplesheet CSV (or pass `--input` on CLI).
	- `module_repeats`, `module_genes`, `module_ncgenes`: enable/disable modules (booleans).
	- `group_consensus_sequences`: (repeats) boolean to share consensus across genomes.
	- `reference_library`: optional repeats library FASTA for repeat annotation.
	- `protein_set`: (genes) path to protein FASTA used for gene annotation.
	- `outdir`: output directory (default `out`).
	- `publish_dir_mode`: how results are published (default `copy`).
	- `max_memory`, `max_cpus`, `max_time`: global resource caps for processes.
	- `env.NXF_SINGULARITY_CACHEDIR`: where Singularity/Apptainer caches images (default `$projectDir/.singularity/`).

### Example runs

Run with your config and samplesheet:
```bash
nextflow run . -c nextflow.config --input samplesheet.csv -profile singularity
```

Use the example config:
```bash
nextflow run . -c example/nextflow_example.config --input example/samplesheet_example.csv -profile singularity
```

To continue after fixing an error, add `-resume`:
```bash
nextflow run . -c nextflow.config --input samplesheet.csv -profile singularity -resume
```

Minimal example samplesheet (samplesheet.csv):
```
name,taxid,fasta
Felis catus,9685,/path/to/felis-catus.fa
```

Column description:

- name (required): display name used in reports (alphanumeric, spaces, '-' and '_' allowed).
- taxid (required): NCBI taxon id (use nearest parent taxid if exact id not available).
- fasta (required): path to fasta assembly file (absolute or relative).
- RNASeq-dir (optional): path to directory containing RNA-Seq FASTQ files for this genome (single-end or paired-end, .fastq or .fastq.gz).
- SRA (optional): comma-separated list of SRA accessions for RNA-Seq data

### Config file

Fill the `nextflow.config` file.

### Launch pipeline

Run the pipeline (example):
```
nextflow run main.nf -profile <singularity|docker|...> --config nextflow.config
```

## Pipeline output
Generated outputs include per-genome annotation files (GFF), fasta sequences and summary reports.

## Credits

Originally written by Dorine Merlat (dorine.merlat@etu.unistra.fr).
Thanks to Arnaud Kress and Odile Lecompte for assistance.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For issues or support, please open an issue on the pipeline GitHub repository.

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  nf-core/test for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->


An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
