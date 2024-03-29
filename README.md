
# AccuraCode
AccuraCode is a bioinfomatics analysis pipeline developed at Singleron to process AccuraCode or AccuraScope sequencing data. This pipeline takes paired-end FASTQ files as input and generates output files which can be used for downstream data analysis.

Detailed docs can be found in [manual](./docs/manual.md).

## NOTICE: Repository No Longer Maintained

Accuracode is no longer actively maintained. We recommend switching to [Celescope](https://github.com/singleron-RD/CeleScope) [bulk_rna assay](https://github.com/singleron-RD/CeleScope/blob/master/doc/assay/multi_bulk_rna.md).

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

Create conda environment and install conda packages. It is recommended to use mamba (which is a faster replacement for Conda):
```
conda install mamba
cd AccuraCode
mamba create -n accuracode -y --file conda_pkgs.txt
```

3. Install python packages and accuracode

Users can specify a mirror source by using the `-i` parameter to speed up the download and installation of packages.
```
conda activate accuracode
pip install -r requirements.txt [ -i https://pypi.tuna.tsinghua.edu.cn/simple ]
python -m pip install .
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


## [New features](./docs/CHANGELOG.md)
