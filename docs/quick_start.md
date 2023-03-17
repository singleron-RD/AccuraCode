# Quick start

AccuraCode contains interfaces `multi_{assay}` to generate pipeline scripts for all assays. Assays can be one of:

- accura

Run `multi_{assay} -h` for help.


## Usage Example

- AccuraCode 384

	```
	conda activate accuracode
	multi_accura\
 	--mapfile ./mapfile\
 	--genomeDir /Database/genome/homo\
 	--thread 8\
 	--mod shell
 	```
`--mapfile` Required. Mapfile path.

`--genomeDir` Required. Genome directory after running `accuracode accura mkref`.

`--thread` The recommended setting is 8, and the maximum should not exceed 20.

`--mod` Create `sjm`(simple job manager https://github.com/StanfordBioinformatics/SJM) or `shell` scripts. 

`--skip_umi_correct` For a large data set, users can set it to skip umi corretion to boost speed.

Scripts above will generate a `shell` directory containing `{sample}.sh` files.

You can start your analysis by running:
```
sh ./shell/{sample}.sh
```

- part of AccuraCode 384
	```
	conda activate accuracode
	multi_accura \
 	--mapfile ./mapfile \
 	--genomeDir /Database/genome/homo \
 	--thread 8 \
 	--mod shell \
	--whitelist /path/to/wellbarcode/file
	```

If you use part of the wells, we recommend you set the "--whitelist" parameter and write the well barcodes you used into a wellbarcode file.

The wellbarcoe file format is one barcode per line, eg:
```
AACGGACCT
GACTGGTTG
TCGGTTCGT
GTCTTGCGT
```




## How to write mapfile

Mapfile is a tab-delimited text file with as least three columns. Each line of mapfile represents paired-end fastq files.

1st column: Fastq file prefix.  
2nd column: Fastq file directory path.  
3rd column: Sample name, which is the prefix of all output files.  

### Example

Sample1 has 2 paired-end fastq files located in 2 different directories(fastq_dir1 and fastq_dir2). Sample2 has 1 paired-end fastq file located in fastq_dir1.
```
$cat ./my.mapfile
fastq_prefix1	fastq_dir1	sample1
fastq_prefix2	fastq_dir2	sample1
fastq_prefix3	fastq_dir1	sample2

$ls fastq_dir1
fastq_prefix1_1.fq.gz	fastq_prefix1_2.fq.gz
fastq_prefix3_1.fq.gz	fastq_prefix3_2.fq.gz

$ls fastq_dir2
fastq_prefix2_1.fq.gz	fastq_prefix2_2.fq.gz
```


 
