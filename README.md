
# AccuraCode
AccuraCode is a bioinfomatics analysis pipeline developed at Singleron to process AccuraCode sequencing data. This pipeline takes paired-end FASTQ files as input and generates output files which can be used for downstream data analysis.

Detailed docs can be found in [manual](./docs/manual.md).

## Hardware/Software Requirements

- minimum 32GB RAM(to run STAR aligner)
- conda
- git

## Installation

1. Clone repo
```
git clone https://github.com/singleron-RD/AccuraCode.git
```

2. Install conda packages
```
cd AccuraCode
conda create -n accuracode
conda activate accuracode
conda install -y --file conda_pkgs.txt --channel conda-forge --channel bioconda --channel r --channel imperial-college-research-computing
```

3. Install accuracode
```
python setup.py install
```


## Reference genome 

### Homo sapiens

```
mkdir hs_ensembl_99
cd hs_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz

gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.99.gtf.gz

conda activate accuracode
accuracode accura mkref \
 --genome_name Homo_sapiens_ensembl_99 \
 --fasta Homo_sapiens.GRCh38.dna.primary_assembly.fa \
 --gtf Homo_sapiens.GRCh38.99.gtf
```

### Mus musculus

```
mkdir mmu_ensembl_99
cd mmu_ensembl_99

wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-99/gtf/mus_musculus/Mus_musculus.GRCm38.99.gtf.gz

gunzip Mus_musculus.GRCm38.dna.primary_assembly.fa.gz 
gunzip Mus_musculus.GRCm38.99.gtf.gz

conda activate accuracode
accuracode accura mkref \
 --genome_name Mus_musculus_ensembl_99 \
 --fasta Mus_musculus.GRCm38.dna.primary_assembly.fa \
 --gtf Mus_musculus.GRCm38.99.gtf
```

## [Quick start](./docs/quick_start.md)


