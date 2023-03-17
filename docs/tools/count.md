## Arguments
`--genomeDir` Required. Genome directory.

`--chemistry` Predefined (pattern, barcode whitelist, linker whitelist) combinations. Can be one of:
- `accuracode96` Used for AccuraCode96 libraries.
- `accuracode384` Used for AccuraCode384 libraries.
- `customized` Used for user defined combinations. You need to provide `pattern`, `whitelist` and  `linker` at the
same time.

`--whitelist` Cell barcode whitelist file path, one cell barcode per line.

`--UMI_cutoff` UMI cutoff for output, default 500.

`--skip_umi_correct` For big data to skip umi correct.

`--outdir` Output diretory.

`--assay` Assay name.

`--sample` Sample name.

`--thread` Thread to use.

`--debug` If this argument is used, accuracode may output addtional file for debugging.

`--bam` Required. BAM file from featureCounts.

