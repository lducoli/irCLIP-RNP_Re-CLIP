import os, sys, collections, re, glob,argparse

############################################################
# Parse arguments.
############################################################

parser = argparse.ArgumentParser(description='Filter out isoforms from a GTF based on RSEM count files.')
parser.add_argument('--gtf', help='Input gtf to filter.')
parser.add_argument('--rsem_dir', help='RSEM counts directory (optional if --rsem_files given)')
parser.add_argument('--rsem_files', help='RSEM counts files (optional if --rsem_dir given)')
parser.add_argument('--out', help='Output GTF filename.')

args = parser.parse_args()
print(args)


"""Input  RSEM_hg38/*.isoforms.results format:
transcript_id	gene_id	length	effective_length	expected_count	TPM	FPKM	IsoPct
ENST00000373020	ENSG00000000003	3768	3500.49	1440.38	12.93	12.03	82.65
ENST00000494424	ENSG00000000003	820	552.49	0.00	0.00	0.00	0.00

TPM: s[5]
"""

#sys.argv = ["", "ref/Homo_sapiens.GRCh38.104.chr.gtf", "RSEM_hg38/", "data/processed/Homo_sapiens.GRCh38.104.chr.filtered.gtf"]
input_gtf = args.gtf
input_counts_dir = args.rsem_dir
input_counts_files = args.rsem_files
output_gtf = args.out

if "," in input_counts_files:
    input_counts_files = input_counts_files.split(',')
print('----', input_counts_files)

############################################################
# Read RSEM counts files.
############################################################

if input_counts_dir is not None:

    genes = {}
    for fname in glob.glob(f"{input_counts_dir}/*.isoforms.results"):
        
        if input_counts_files is not None and fname not in input_counts_files:
            print(f"Skipped using {fname} to prefilter {input_gtf} because it wasn't in the whitelist: {input_counts_files}")
            continue
        print(f'~ {fname}')
        with open(fname) as f:
            next(f)  # Skip header.
            for li in f:
                s = li.split('\t')
                genes.setdefault(s[1], collections.defaultdict(float))
                genes[s[1]][s[0]] += float(s[5])
                
elif input_counts_files is not None:
    genes = {}
    for fname in input_counts_files:
        
        if input_counts_files is not None and fname not in input_counts_files:
            print(f"Skipped using {fname} to prefilter {input_gtf} because it wasn't in the whitelist: {input_counts_files}")
            continue
        print(f'~ {fname}')
        with open(fname) as f:
            next(f)  # Skip header.
            for li in f:
                s = li.split('\t')
                genes.setdefault(s[1], collections.defaultdict(float))
                genes[s[1]][s[0]] += float(s[5])
else:
    raise IOError(f"Error: must set at least one of --rsem_dir or --rsem_files.")

############################################################
# Decide what isoforms to keep.
############################################################

totals = collections.defaultdict(float)
genes_fractions = {}
to_keep = set()
all_isoforms = set()
cutoff = 0.05

for gene, txpt_to_tpm in genes.items():
    totals[gene] = sum(txpt_to_tpm.values())
    if totals[gene] == 0:
        continue
    genes_fractions[gene] = {k:v/totals[gene] for k,v in txpt_to_tpm.items()}
    to_keep |= set([name for name,frac in genes_fractions[gene].items() if frac>=cutoff])
    all_isoforms |= set(txpt_to_tpm.keys())

print(f"Started with {len(all_isoforms)} across {len(genes)} genes. Removed those with " + \
      f"less than {cutoff:.2%} of TPM at that ENSG. Kept {len(to_keep)} " + \
      f"({len(to_keep)/len(all_isoforms):.2%}).")

############################################################
# Filter the gtf based on the to_keep set of OK isoforms.
############################################################

outf = open(output_gtf, 'w')

lines_discarded, lines_kept, lines_discarded_bc_no_enst = 0, 0, 0
with open(input_gtf) as f:
    for li in f:
        if li[0] == '#':
            outf.write(li)
            continue
        enst = re.search('transcript_id "([^"]+)"', li)
        if enst is None:
            lines_discarded_bc_no_enst += 1
            continue
        if enst.group(1) in to_keep:
            outf.write(li)
            lines_kept += 1
        else:
            lines_discarded += 1
            
lines = lines_discarded + lines_kept + lines_discarded_bc_no_enst

print(f"Filtering gtf: kept {lines_kept} ({lines_kept/lines:.2%}) lines out of {lines} lines." + \
      f" Discarded {lines_discarded} as having an ENST below cutoff, and {lines_discarded_bc_no_enst} " + \
      "for having no ENST in the GTF.")

