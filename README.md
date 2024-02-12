# EH5_kmerfilter
k-mer based filtering of ExpansionHunter 5 output.

# Changelog

## v0.12.0
- --save_filtered flag added to save removed calls in a separate output VCF (<output>_removed.vcf)
- --auto flag added to allow auto margin resizing on genotypes below the read length. This allows the separate use of both --rank and --auto in one command. 
- Revise kmer_mul argument to take a float in range 0.0-1.0. Default value = 0.2
- Update margin tightening for genotyped motif lengths < read length

## v0.11.0

- new handling for motifs < read length - if a motif is < the read length the margin size is set to the locus + 20% of the genotyped motif length. This reduces the left-censoring behaviour the filtering was exhibiting at certain loci
- reference-sized STRs should now be correctly reported, unless they occur as the result of a LowDepth call - ie genotypes of ./. in the format string will be filtered, but genotypes of <ref>/<ref> will not.

## v0.10.0

- fixes `--rank` parameter
- adds `auto` option to `--rank`. This is recommended for most uses. `auto` will act as `--rank 1` on loci with only one motif defined, but at motifs where *n* loci are defined, it will act as `--rank n`.

# Installation

## Via pip

Ensure you have Python 3.6 or later installed. You can download Python from the official website: https://www.python.org/

To install Pysam, use pip:

```
pip install pysam
```

If you need to install python libraries to user install directories:

```
pip install pysam --user
```

## Via conda

Alternatively, install Conda, then:

```
conda env create -f environment.yml
source activate eh5_kmer
```

# Usage

```
python kmer_filter.py -b <bam_path> -c <catalog_path> -v <vcf_path> -o <output_path> -m <margin> -r <kmer_mul> --logs
```

## Arguments

* -b, --bam: Path to the EH5 BAM file.
* -c, --catalog: Path to the JSON catalog file used by EH5.
* -v, --vcf: Path to the EH5 VCF file.
* -o, --output: Path for the output files. The script will generate <output>_validated.vcf and optionally <output>_kmers.tsv.gz.
* -m, --margin: Margin for fetching sequence reads from the BAM file (default is 1,000).
* --rank: integer. k-mer rank; will accept a call if it is in the top (rank) of kmers at the locus in BAMs (default = 1).
* --auto: boolean. If present, will act as --rank n at motifs where n loci are defined.
* --keep_lowdepth: boolean. Flag to keep LowDepth calls in VCF
* -r, --kmer_mul: Multiplier for the k-mer size (float in range 0.0 - 1.0, default is 0.2).
* --logs: boolean. If present, will save <output>_kmers.tsv.gz file containing all kmer counts in the BAM.
* --save_filtered: boolean. If present, will save removed calls in a separate file <output>_removed.vcf

## Example

To process a single sample (here imaginatively called 'sample'):

Without the log containing all kmer counts:
```
python kmer_filter.py -b sample.bam -c catalog.json -v sample.vcf -o ./sample_kmer
```

With the log containing all kmer counts:
```
python kmer_filter.py -b sample.bam -c catalog.json -v sample.vcf -o ./sample_kmer --logs
```
