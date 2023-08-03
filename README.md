# EH5_kmerfilter
k-mer based filtering of ExpansionHunter 5 output.

# Installation

## Via pip

Ensure you have Python 3.6 or later installed. You can download Python from the official website: https://www.python.org/

To install Pysam, use pip:

```
pip install pysam
```

## Via conda

Alternatively, install Conda, then:

```
conda create env -f environment.yml
source activate eh5_kmer
```

# Usage

```
python kmer_filter.py -b <bam_path> -c <catalog_path> -v <vcf_path> -o <output_path> -m <margin> -r <kmer_mul>
```

## Arguments

* -b, --bam: Path to the EH5 BAM file.
* -c, --catalog: Path to the JSON catalog file used by EH5.
* -v, --vcf: Path to the EH5 VCF file.
* -o, --output: Path for the output files. The script will generate <output>_validated.vcf and <output>_kmers.tsv.
* -m, --margin: Margin for fetching sequence reads from the BAM file (default is 1000).
* -r, --kmer_mul: Multiplier for the k-mer size (default is 15).

## Example

To process a single sample (here imaginatively called 'sample'):

```
python kmer_filter.py -b sample.bam -c catalog.json -v sample.vcf -o ./sample_kmer
```
