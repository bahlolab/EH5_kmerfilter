import argparse
import json
import re
import pysam
from collections import defaultdict
from math import floor


REVERSE_COMPLEMENTER = str.maketrans("ACGT", "TGCA")
IUPAC_REGEX_DICT = {"R": "[AG]", "Y": "[TC]", "M": "[AC]", "K": "[GT]", "S": "[GC]", "W": "[AT]", "H": "[ACT]",
                    "B": "[GCT]",
                    "V": "[ACG]", "D": "[GAT]", "N": "[ACGT]"}


def lex_min_rotation(motif, reverse_complement=False):
    """Calculate the lexicographically minimal rotation of the supplied motif."""
    mlen = len(motif)
    motif = motif * 2
    cycles = [motif[i:mlen + i] for i in range(mlen)]
    if reverse_complement:
        rcomp = motif.translate(REVERSE_COMPLEMENTER)[::-1]
        rcomp = rcomp * 2
        cycles += [rcomp[i:mlen + i] for i in range(mlen)]
    cycles.sort()
    return cycles[0]

def kmer_count(seq_list, k):
    """Count the occurrences of kmers in a list of sequences."""
    kmer_counts = defaultdict(int)
    for seq in seq_list:
        seqlen = len(seq)
        for i in range(seqlen - k + 1):
            kmer = seq[i:(i + k)]
            kmer_counts[kmer] += 1
    return kmer_counts

def get_motif_regex(iupac_string):
    """Convert an IUPAC string to a regex pattern."""
    return "".join([IUPAC_REGEX_DICT.get(x, x) for x in iupac_string])


def parse_vcf_file(vcf_path):
    """Parse an EH5 VCF file and return a dictionary with relevant info."""
    vcf_catalog = {}
    with open(vcf_path, 'rt') as vcf_file:
        for line in vcf_file:
            # Skip LowDepth calls and headers.
            if line.startswith("#") or "LowDepth" in line:
                continue
            spline = line.strip().split()
            desc_line = spline[7]
            var_id = None
            rep_id = None
            for split in desc_line.split(";"):
                if "VARID=" in split: var_id = split.replace("VARID=", "")
                if "REPID=" in split: rep_id = split.replace("REPID=", "")
            id_string = var_id if var_id else rep_id
            rep_counts = [int(x) for x in spline[9].split(":")[2].split("/")]
            max_allele = max(rep_counts)
            min_allele = min(rep_counts)
            if id_string in vcf_catalog:
                print("ERROR: Overwriting ID string.")
                exit(1)
            vcf_catalog[id_string] = (var_id, rep_id, max_allele, min_allele, spline[3])
    return vcf_catalog

def process_catalog(catalog, vcf_catalog, samfile, motif_dict, regex_dict, kmer_range=1):
    retain_ids = set()
    with open(args.output_path + "_kmers.tsv", 'wt') as kmer_out:
        for val in catalog:
            locus_id_list = val["LocusId"]
            locus_coord = val["ReferenceRegion"]
            if isinstance(locus_coord, list):
                locus_id_list = val["VariantId"]
                locus_struct_list = [x for x in re.findall("\((\w+)\)", val["LocusStructure"])]
            else:
                locus_id_list = [locus_id_list]
                locus_struct_list = [val["LocusStructure"]]
                locus_coord = [locus_coord]
            for motif_id, motif, coord in zip(locus_id_list, locus_struct_list, locus_coord):
                if motif_id not in vcf_catalog:
                    continue
                motif_dict[motif_id] = lex_min_rotation(motif.replace("(", "").replace(")", "").replace("*", ""))
                if motif_dict[motif_id] not in regex_dict:
                    regex_dict[motif_dict[motif_id]] = re.compile(get_motif_regex(motif_dict[motif_id]))
                # This was using a more succint bit of code to get read_sequences, expanded for maintainability. Original was:
                # read_sequences = [read.seq for read in samfile.fetch(*coord.split(":")[0], int(coord.split(":")[1].split("-")[0]) - args.margin, int(coord.split(":")[1].split("-")[1]) + args.margin)]
                read_sequences = list()
                chrom, positions = coord.split(":")
                start, end = [int(x) for x in positions.split("-")]
                for read in samfile.fetch(chrom, start - args.margin, end + args.margin):
                    tag = read.get_tag("XG").split(",")[0]
                    read_sequences.append(read.seq)
                if len(read_sequences) == 0:
                    continue
                genotyped_len = vcf_catalog[motif_id][2]
                read_len = [len(x) for x in read_sequences]
                average_read_len = floor(sum(read_len) / len(read_len))
                kmer_size = max(min(floor(0.2 * genotyped_len) * len(motif_dict[motif_id]),
                                    floor(0.2 * floor(average_read_len / len(motif_dict[motif_id]))) * len(
                                        motif_dict[motif_id])),
                                len(motif_dict[motif_id]))
                kmer_dict = kmer_count(read_sequences, kmer_size)
                lexed_dict = defaultdict(int)
                for kmer in kmer_dict:
                    if kmer.count(kmer[0:len(motif_dict[motif_id])]) * len(motif_dict[motif_id]) == kmer_size:
                        lexed_dict[lex_min_rotation(kmer)[0:len(motif_dict[motif_id])]] += kmer_dict[kmer]
                sort_kmers = [(k, v) for k, v in sorted(lexed_dict.items(), reverse=True, key=lambda item: item[1])]
                freq_kmer_string = " ".join([f'{kmer[0]}:{kmer[1]}' for kmer in sort_kmers[0:1]])
                if any(regex_dict[motif_dict[motif_id]].search(freq_kmer_string) for _ in range(kmer_range)):
                    retain_ids.add(motif_id)
                kmer_out.write("\t".join([motif_id, motif, coord, str(genotyped_len * len(motif)), str(kmer_size),
                                          ";".join([f'{kmer}:{lexed_dict[kmer]}' for kmer in lexed_dict])]))
                kmer_out.write("\n")
    return retain_ids

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", dest="bam_path", help="Path to EH5 BAM output.")
    parser.add_argument("-c", "--catalog", dest="catalog_path", help="Path to EH5 variant catalog")
    parser.add_argument("-v", "--vcf", dest="vcf_path", help="Path to EH5 VCF")
    parser.add_argument("-o", "--output", dest="output_path", help="Prefix (path) for output")
    parser.add_argument("-m", "--margin", dest="margin", default=1000, help="Margin (in nt) around catalog entry to calculate kmers")
    parser.add_argument("-r", "--kmer_mul", dest="kmer_mul", default=15, help="k-mer multiplier; will count repeats of kmer_mul*motif_len")
    parser.add_argument("--rank", dest="rank", default=1, help="k-mer rank; will accept a call if it is in the top (rank) of kmers at the locus in BAMs")
    args = parser.parse_args()
    with open(args.catalog_path, 'rt') as in_file:
        catalog = json.load(in_file)
    vcf_catalog = parse_vcf_file(args.vcf_path)
    samfile = pysam.AlignmentFile(args.bam_path, "rb")
    tagset = set()
    motif_dict = {}
    regex_dict = {}
    retain_ids = process_catalog(catalog, vcf_catalog, samfile, motif_dict, regex_dict)
    with open(args.vcf_path, 'rt') as vcf_file, open(args.output_path + "_validated.vcf", 'wt') as out_file:
        for line in vcf_file:
            if line.startswith("#"):
                out_file.write(line)
                continue
            spline = line.strip().split()
            desc_line = spline[7]
            var_id = None
            rep_id = None
            for split in desc_line.split(";"):
                if "VARID=" in split: var_id = split.replace("VARID=", "")
                if "REPID=" in split: rep_id = split.replace("REPID=", "")
            id_string = var_id if var_id else rep_id
            if id_string in retain_ids:
                out_file.write(line)
