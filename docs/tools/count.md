## Features
- Get read/UMI/gene counts for each well barcode.

- Generate expression matrix.

## Output
- `{sample}_matrix.tsv.gz` The expression matrix of the well barcode, separated by tabs.

- `{sample}_count_detail.txt.gz` 4 columns:
    - barcode
    - gene ID
    - UMI count
    - read_count

- `{sample}_counts.txt` 6 columns:
    - well: Well barcode sequence
    - readcount: read count of each well barcode
    - UMI: UMI count for each well barcode
    - gene: gene count for each well barcode

- `{sample}_downsample.txt` 2 columns：
    - percent: percentage of sampled reads
    - median_geneNum: median gene number per cell




## Arguments
`--genomeDir` Required. Genome directory.

`--chemistry` Required. Default 'accuracode384'.

`--whitelist` Cell barcode whitelist file path, one cell barcode per line.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, accuracode may output addtional file for debugging.

`--bam` Required. BAM file from featureCounts.


