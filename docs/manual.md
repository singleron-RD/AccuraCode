## Introduction
AccuraCode is a bioinfomatics analysis pipeline developed at Singleron to AccuraCode sequencing data generated with Singleron products. These pipelines take paired-end FASTQ files as input and generate output files which can be used for downstream data analysis.

Each pipeline consists of several steps and they all have two identical pre-processing steps: `barcode` and `cutadapt`. `barcode`step is used for barcode demupltiplexing, correction and read filtering. `cutadapt`step calls [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) for read trimming.

Currently, AccuraCode includes the follwing pipelines:

- `accuracode accura` for Multi samples RNA-seq data generated with AccuraCode kits. It performs preprocessing, genome alignment, feature counting, expression matrix generation.



## [Quick start](quick_start.md)

## [Change log](CHANGELOG.md)

## Pre-processing

- [barcode](tools/barcode.md)
- [cutadapt](tools/cutadapt.md)

## accura
- [mkref](accura/mkref.md)
- [star](accura/star.md)
- [featureCounts](tools/featureCounts.md)
- [count](tools/count.md)
- [multi_accura](accura/multi_accura.md)

