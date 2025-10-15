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

Example samplesheet (samplesheet.csv):
```
name,taxid,fasta
Polydesmus complanatus,510027,/path/to/polydesmus-complanatus.fa
```

Column description:

- name (required): display name used in reports (alphanumeric, spaces, '-' and '_' allowed).
- taxid (required): NCBI taxon id (use nearest parent taxid if exact id not available).
- fasta (required): path to fasta assembly file (absolute or relative).

### Config file

Fill the `nextflow.config` file. Here a describe the main options

*Note*: Use absolute paths when targeting cloud storage (S3/GS).

**Input & output**

| Parameter | Type | Default | Description |
|---|---:|---|---|
| input | string | null | Path to samplesheet CSV describing genomes to annotate (required). |
| outdir | string | 'out' | Output directory for pipeline results. |

**Module activation**

| Parameter | Type | Default | Description |
|---|---:|---|---|
| module_repeats | boolean | true | Enable repetitive elements annotation. |
| module_genes | boolean | true | Enable protein-coding gene annotation. |
| module_ncgenes | boolean | true | Enable non-coding gene annotation. |

**Options for repetitive element annotation**

| Parameter | Type | Default | Description |
|---|---:|---|---|
| download_famdb | boolean | false | true = download a FamDB partition from dfam.org; false = use local `famdb_path`. |
| famdb_to_download | string | null | FamDB partition id ('1'..'16') to download when `download_famdb` is true. See [FamDB](https://github.com/dorinemerlat/exogap/README.md#####FamDB) |
| famdb_path | string | null | Local path to FamDB `.h5` file (required when `download_famdb` is false). |
| group_consensus_sequences | boolean | true | Use a single consensus library across all genomes in the run. |

**Resources & limits**

| Parameter | Type | Default | Description |
|---|---:|---|---|
| max_memory | string | '128.GB' | Upper memory limit for jobs (Nextflow format, e.g. "8.GB"). |
| max_cpus | integer | 16 | Upper CPU limit for jobs. |
| max_time | string | '240.h' | Upper walltime limit for jobs (e.g. "2.h", "30.m"). |

#### FamDB

Choose the FamDB partition id (1..16) matching your taxonomic scope when using `download_famdb`.

| Partition | Name | Content |
|---:|---|---|
| 0 | root | root |
| 1 | Brachycera | Brachycera |
| 2 | Archelosauria | Archelosauria |
| 3 | Hymenoptera | Hymenoptera |
| 4 | Otomorpha | Otomorpha |
| 5 | rosids | rosids |
| 6 | Viridiplantae | Saxifragales, asterids, Proteales, Nymphaeales, Amborellales, Caryophyllales, Ranunculales, Mesostigmatophyceae, Chlorokybophyceae, Charophyceae, Lycopodiopsida, Chlorophyta, Liliopsida, Polypodiopsida, Marchantiophyta, Acrogymnospermae, Bryophyta |
| 7 | Mammalia | Mammalia |
| 8 | Noctuoidea | Noctuoidea |
| 9 | Obtectomera | Bombycoidea, Papilionoidea, Pyraloidea, Hesperioidea, Geometroidea, Drepanoidea, Pterophoroidea |
| 10 | Eupercaria | Eupercaria |
| 11 | Ctenosquamata | Ovalentaria, Myctophata, Lampridacea, Carangaria, Holocentrimorphaceae, Batrachoidaria, Anabantaria, Paracanthopterygii, Ophidiaria, Gobiaria, Syngnathiaria, Pelagiaria |
| 12 | Vertebrata | Chondrichthyes, Lepidosauria, Protacanthopterygii, Coelacanthimorpha, Amphibia, Cladistia, Holostei, Cyclostomata, Osteoglossocephala, Stomiati, Dipnomorpha, Elopocephalai, Chondrostei |
| 13 | Coleoptera | Coleoptera |
| 14 | Endopterygota | Gelechioidea, Yponomeutoidea, Incurvarioidea, Tineoidea, Apoditrysia, Nematocera, Strepsiptera, Neuropterida, Siphonaptera, Trichoptera |
| 15 | Protostomia | Nematoda, Chelicerata, Collembola, Polyneoptera, Monocondylia, Palaeoptera, Crustacea, Paraneoptera, Myriapoda, Scalidophora, Spiralia |
| 16 | Riboviria | Diverse groups including many unicellular eukaryotes, viruses and bacteria (see original dfam partition descriptions) |

### Launch pipeline

Run the pipeline (example):
```
nextflow run main.nf -profile <singularity|docker|...> --config nextflow.config
```

## Pipeline output
Generated outputs include per-genome annotation files (GFF), fasta sequences and summary reports and MultiQC reports. See docs/ for detailed output layout (TODO: add full output list).

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
