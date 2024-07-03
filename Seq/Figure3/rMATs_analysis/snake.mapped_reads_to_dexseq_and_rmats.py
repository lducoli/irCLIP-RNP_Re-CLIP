import re, gzip, glob, csv, pandas, collections, functools, subprocess

########################################################################
# Globals.
########################################################################

wildcard_constraints:
    directory=".+\/",
    condition="[^\/]+",
    sample="[^\/]+",
    comparison="[^\/]+",
    
# /oak/stanford/groups/khavari/users/dfporter/seq/genomes/star_gencode_39/
# python scripts/filter_isoforms_using_rsem_results.py --gtf /oak/stanford/groups/khavari/users/dfporter/seq/genomes/gencode_39_primary_files/gencode.v39.annotation.gtf --rsem_dir data/rsem_counts --out data/processed/genome.filtered_by_rsem.gtf --rsem_files data/rsem_counts/Empty_pLEX-anti-FUBP1_shRNA-Batch_3-293T_cells-SRR13169583.isoforms.results

GTF = config['GTF'] #"/oak/stanford/groups/khavari/users/dfporter/seq/genomes/gencode_39_primary_files/gencode.v39.annotation.gtf"
GENOME_FA = config['GENOME_FA']  #"/oak/stanford/groups/khavari/users/dfporter/seq/genomes/gencode_39_primary_files/GRCh38.primary_assembly.genome.fa"
BAM_DIR = config['BAM_DIR']  # = "STAR_genome_hg38/transcriptomes"

GENOMIC_BAM_DIR = "star"
# These are created if they don't exist:
RSEM_REF = config['RSEM_REF']  # "rsem_reference/GRCh38.gencode39.all.chrom.reference"
RSEM_COUNTS_DIR = config['RSEM_COUNTS_DIR']  #"data/rsem_counts"

SAMPLE_TABLE = config['SAMPLE_TABLE'] #sample_table.txt'
comparisons = config['comparisons']
COMPARISONS = ['_vs_'.join(_list) for _list in comparisons]
SAMPLES = [os.path.basename(x).split('_Aligned')[0] for x in glob.glob(BAM_DIR + '/*bam')]

print(config)

########################################################################
# Functions.
########################################################################

def samples_for_comparison(comparison):
    df = pandas.read_csv(SAMPLE_TABLE, sep='\t')
    groups = comparison.split('_vs_')
    return df.loc[[x in groups for x in df.group], 'sample'].values  

def counts_files_for_comparison(comparison):
    return [RSEM_COUNTS_DIR + f"/{samp}.isoforms.results" for samp in samples_for_comparison(comparison)]

def wc_counts_files_for_comparison(wildcards):
    return [RSEM_COUNTS_DIR + f"/{samp}.isoforms.results" for samp in samples_for_comparison(wildcards.comparison)]

def space_delim_counts_files_for_comparison(comparison):
    return " ".join(counts_files_for_comparison(comparison))

def comma_delim_counts_files_for_comparison(comparison):
    return ",".join(counts_files_for_comparison(comparison))

def space_delim_bam_files_for_comparison(wildcards):
    print(f"space_delim_bam_files... arg={wildcards}")
    samps = samples_for_comparison(wildcards.comparison)
    print(f"~~~~~ {samps}")
    return " ".join([f"{GENOMIC_BAM_DIR}/{_sample}/{_sample}_Aligned.out.bam" for _sample in samps])

def wc_bam_files_for_comparison(wildcards):
    print(f'wc_bam_files... arg={wildcards}')
    samps = samples_for_comparison(wildcards.comparison)
    print(f'------ {samps}')
    return [f"{GENOMIC_BAM_DIR}/{_sample}/{_sample}_Aligned.out.bam" for _sample in samps]

def bam_files_for_group(group):
    # Convert all these bam file finders to using a column from the sample sheet.
    df = pandas.read_csv(SAMPLE_TABLE, sep='\t')
    samps = df.loc[[x==group for x in df.group], 'sample'].values
    return [f"{GENOMIC_BAM_DIR}/{_sample}/{_sample}_Aligned.out.bam" for _sample in samps]

def wc_bam_files_for_group(wildcards):
    return bam_files_for_group(wildcards.group)

"""
conda env create -f env/conda.yaml -n splicing

1. Download reference
2. Build reference
3. map_to_rsem_reference
4. (map_to_rsem_reference -> ) count_with_rsem
5. (count_with_rsem -> ) prefilter_isoforms_using_rsem_results [-> makes new gtf]
6. [filtered gtf -> ] format_gtf_for_Rsubread [-> makes new gtf]
7. (map_to_rsem_reference, formatted gtf -> ) count_exons_with_Rsubread
8. (count_with_rsem -> ) differential genes (DESeq2)
9. (count_exons_with_Rsubread -> ) differential exons (DEXSeq)
"""

########################################################################
# Rules.
########################################################################

rule all:
    input:
        #"rsem_reference/GRCh38.gencode39.all.chrom.reference.seq",
        #expand(RSEM_COUNTS_DIR + "/{sample}.isoforms.results", sample=SAMPLES),
        expand("data/processed/{comparison}/subread_fc.out", comparison=COMPARISONS),
        expand("outs/dexseq/dexseq.{comparison}.with_bio_exons.xlsx", comparison=COMPARISONS),
        "outs/dexseq/dexseq.combined_table.xlsx",
        expand("data/processed/rmats_inputs/{comparison}.sbatch", comparison=COMPARISONS),
    shell:
        "echo Worflow complete!"

        
########################################################################
# RMATS.
########################################################################
# samtools view STAR_genome_hg38/Empty_pLEX-anti-FUBP1_shRNA-Batch_3-293T_cells-SRR13169583/Empty_pLEX-anti-FUBP1_shRNA-Batch_3-293T_cells-SRR13169583_Aligned.out.bam | awk '{print length($10)}' | head -1



rule get_read_length:
    input:
        bam_list = wc_bam_files_for_comparison
    output:
        txt = "data/processed/rmats_inputs/bam_file_read_length_{comparison}.txt"
    run:
        cmd = f"samtools view {input.bam_list[0]} | head -n 1 |" + " awk '{print length($10)}' "

        read_len_string = subprocess.getoutput(cmd)
        
        with open(str(output.txt), 'w') as f:
            f.write(read_len_string)

# Functions for the call_rmats rule.
wc_groups_for_comparison = lambda wc: wc.comparison.split('_vs_')
wc_group_1 = lambda wc: wc_groups_for_comparison(wc)[0]
wc_group_2 = lambda wc: wc_groups_for_comparison(wc)[1]
wc_bam_list_1 = lambda wc: bam_files_for_group(wc_group_1(wc))
wc_bam_list_2 = lambda wc: bam_files_for_group(wc_group_2(wc))

wc_bam_file_list_fname_1 = lambda wc: f"data/processed/rmats_inputs/bam_lists/{wc_group_1(wc)}.txt"
wc_bam_file_list_fname_2 = lambda wc: f"data/processed/rmats_inputs/bam_lists/{wc_group_2(wc)}.txt"


rule write_bam_file_lists_for_rmats:
    input:
        bams1 = wc_bam_list_1,
        bams2 = wc_bam_list_2,
    output:
        # This isn't used, it just makes it easier for snakemake to have a pipeline. The individual files are used.
        bam_list_all = "data/processed/rmats_inputs/bam_lists/all_bam_files_{comparison}.txt"
    run:
        bam_list_1 = wc_bam_file_list_fname_1(wildcards)
        bam_list_2 = wc_bam_file_list_fname_2(wildcards)
        with open(str(bam_list_1), 'w') as f:
            f.write(  ",".join(input.bams1)  )
        with open(str(bam_list_2), 'w') as f:
            f.write(  ",".join(input.bams2)  )
        with open(str(output.bam_list_all), 'w') as f:
            f.write(  ",".join(input.bams1) + "," + ",".join(input.bams2) )
            

rule call_rmats:
    input:
        bam_list_all = "data/processed/rmats_inputs/bam_lists/all_bam_files_{comparison}.txt",
        read_length = "data/processed/rmats_inputs/bam_file_read_length_{comparison}.txt",
    output:
        sbatch = "data/processed/rmats_inputs/{comparison}.sbatch",
        rmats_outdir = directory("data/processed/rmats/{comparison}/")
    run:
        read_len = open(str(input.read_length)).readlines()[0]
        bam_list_1 = wc_bam_file_list_fname_1(wildcards)
        bam_list_2 = wc_bam_file_list_fname_2(wildcards)

        with open(str(output.sbatch), 'w') as f:
            #f.write(header)
            g1,g2 = bam_list_1, bam_list_2
            out = str(output.rmats_outdir)#"data/processed/rmats/" #+ "_vs_".join()
            os.makedirs(out, exist_ok=True)
            prefix = f"rmats.py --b1 {g1} --b2 {g2} --gtf {GTF} -t paired  --readLength {read_len} --nthread 8 --od {out} --tmp {out}/tmp"
            f.write(prefix + " --task prep ;\n")
            f.write(prefix + " --task post ;\n")
            f.write(prefix + " --task stat ;\n")
            print(prefix)
            
########################################################################
# Later annotations of DEXSeq results.
########################################################################
"""
The featureCounts formatted gtf - 
"211014_rescue/data/processed/Homo_sapiens.GRCh38.104.chr.filtered.fc.gtf"
- is produced from 211014_rescue/data/processed/Homo_sapiens.GRCh38.104.chr.filtered.gtf
by this shell command produced by snake.mapping_rsubread_rsem_deseq.py 
(rule format_gtf_for_Rsubread):

    python Subread_to_DEXSeq/dexseq_prepare_annotation2.py --aggregate no \
        --featurecountsgtf  data/processed/Homo_sapiens.GRCh38.104.chr.filtered.fc.gtf \
        data/processed/Homo_sapiens.GRCh38.104.chr.filtered.gtf  \
        data/processed/Homo_sapiens.GRCh38.104.chr.filtered.fc.gff

This gtf is then used in
scripts/DEXSeq_of_exons.generalized.R

Specifically, it's used in this function call, as the gtf_fname:
    dxd <- DEXSeqDataSetFromFeatureCounts(
        subread_fname, flattenedfile = gtf_fname, sampleData = sampleData,
        design = ~ sample + (exon + batch) * condition)
        
This "fc" gtf changes the genome annotation from data/processed/Homo_sapiens.GRCh38.104.chr.filtered.gtf
quite substantially. Namely, it splits each of the real exons in the input gtf
into subexons. The output files all reflect these subexons, which do not necessarily
actually exist.

We therefore have to recover the real exons from the subexons.
"""

def reindex(_df):
    _df.index = list(zip(_df['groupID'], _df['featureID'], _df['chrm'], _df['start'], _df['end']))
    
def format_txpts(txpts):
    """Turn the strings output by R into sets of strings."""
    if type(txpts)==type(set()):
        return txpts
    if txpts[:2] == "c(":
        txpts = set([x.strip('"') for x in txpts[2:-1].split(', ')])    
    else:
        txpts = set([txpts])
    return txpts

# For each subexon or exon, it can be defined as: txpt, iv
# To assign a subexon to an exon, we subset to the txpt, the to those iv containing the
# subexon iv.
def read_original_gtf(gtf):
    # enst -> (chrm, start, end, strand, exon_number)
    enst_to_ex = collections.defaultdict(list)
    with open(gtf) as f:
        for li in f:
            if li[0] == '#':
                continue
            #print(li)
            s = li.split('\t')

            if s[2] == 'exon':
                #tx = re.search('transcripts "([^"]+)"', s[-1])
                ex = re.search('exon_number "([^"]+)"', s[-1])
                if ex is None:
                    ex = re.search('exon_number ([^;]+);', s[-1])
                gene = re.search('gene_id "([^"]+)"', s[-1])
                tx = re.search('transcript_id "([^"]+)"', s[-1])
                if tx is not None:
                    _txpt = tx.group(1)
                    #ex = ex.group(1)
                    enst_to_ex[_txpt].append((s[0], int(s[3]), int(s[4]), s[6], ex.group(1)))
    return enst_to_ex

def read_fc_gtf(gtf):
    # enst -> (chrm, start, end, strand, exon_number)
    enst_to_ex = collections.defaultdict(dict)
    with open(gtf) as f:
        for li in f:
            if li[0] == '#':
                continue

            s = li.split('\t')

            if s[2] == 'exon':
                tx = re.search('transcripts "([^"]+)"', s[-1])
                ex = re.search('exon_number "([^"]+)"', s[-1])
                #gene = re.search('gene_id "([^"]+)"', s[-1])
                #tx = re.search('transcript_id "([^"]+)"', s[-1])
                if tx is not None:
                    _txpts = tx.group(1)
                    
                    for txpt in _txpts.split('+'):
                        txpt = txpt.lstrip('c(').rstrip(')')
                        enst_to_ex[txpt][ex.group(1)] = (s[0], int(s[3]), int(s[4]), s[6], ex.group(1))
    return enst_to_ex

def subexon_to_biological_exon_converter(enst_to_ex, enst_to_fc):
    # Create a converter from the subexons to the original, biological exons.
    enst_iv_to_orig_exon = collections.defaultdict(dict)
    for enst, exonsD in enst_to_fc.items():

        for exon_n, sub_ex_iv in exonsD.items():  # For each subexon.
            matches = set()
            for orig_iv in enst_to_ex[enst]:  # For each iv in the original exon.
                if (orig_iv[1] <= sub_ex_iv[1] <= orig_iv[2]) and (
                    orig_iv[1] <= sub_ex_iv[2] <= orig_iv[2]):
                    matches.add(orig_iv)  # This subexon belongs to the original exon orig_iv.
                    #print(f"Match: {sub_ex_iv} is within {orig_iv}")
                    # The subexon is within this exon.
            if len(matches) == 1:
                match = list(matches)[0]
                enst_iv_to_orig_exon[enst][sub_ex_iv] = match
            else:
                raise Exception  # For a single transcript, can't have overlapping exons.
                
    return enst_iv_to_orig_exon

        #########################################################################
        # Load the DEXSeq results txt files into a dict of dataframes.
        # Add biological exons.
        #########################################################################

def convert_subexons_to_biological_exons(df, enst_to_ex, enst_iv_to_orig_exon):
    def convert(txpts, exon_n, chrm, start, end, strand):
        txpts = format_txpts(txpts)
        original_exons = []
        for txpt in txpts:
            
            iv = (chrm, start, end-1, strand, exon_n.lstrip('E'))
            n_exons = len(enst_to_ex[txpt])
            exon_positions = sorted(enst_to_ex[txpt], key=lambda x: x[1])
            this_ex_pos = 1 + exon_positions.index(enst_iv_to_orig_exon[txpt][iv])
            
            fraction = (this_ex_pos)/len(exon_positions)
            fraction = fraction if strand!='-' else 1-fraction
            
            this_ex_pos = this_ex_pos if strand!='-' else 1+n_exons-this_ex_pos
            
            original_exons.append((txpt, n_exons, this_ex_pos, fraction, *enst_iv_to_orig_exon[txpt][iv]))
        #if random.randint(1,1000) == 1:
        #    print(txpts, original_exons, enst_iv_to_orig_exon[txpt][iv])
        return original_exons
    
    df['biological_exons'] = [convert(*tup) for tup in zip(
        df["transcripts"], df['featureID'], df['chrm'], df['start'], df['end'], df['strand'])]

    return df

def split_biological_exon_info_into_columns(df):
    df['# exons'] = [arr[0][1] for arr in df['biological_exons']]
    df['exon #'] = [arr[0][2] for arr in df['biological_exons']]
    df['Pos in gene'] = [arr[0][3] for arr in df['biological_exons']]
    return df

rule recover_biological_exons:
    input:
        original_gtf = GTF,
        fc_gtf = "data/processed/{comparison}/genome.filtered_by_rsem.fc.gtf",
        dex_results = "data/processed/{comparison}/dexseq/DEXSeq_results.annotated.txt", 
    output:
        xlsx = "outs/dexseq/dexseq.{comparison}.with_bio_exons.xlsx"
    run:
        enst_to_ex = read_original_gtf(str(input.original_gtf))
        enst_to_fc = read_fc_gtf(str(input.fc_gtf))
        enst_iv_to_orig_exon = subexon_to_biological_exon_converter(enst_to_ex, enst_to_fc)
        # Fill the biological exons column with:
        # [(txpt, n_exons, this_ex_pos+1, fraction, chrom, start, end, strand, exon_n), ...]
        fname = str(input.dex_results)
        
        name = os.path.basename(fname).split('DEXSeq_results_')[-1].split('.txt')[0]
        df = pandas.read_csv(fname, sep='\t')
        df = convert_subexons_to_biological_exons(df, enst_to_ex, enst_iv_to_orig_exon)
        print(df)
        reindex(df)
        df = split_biological_exon_info_into_columns(df)
        
        df.to_excel(str(output.xlsx), engine='openpyxl')

########################################################################
# DEXSeq.
########################################################################

def space_sep_groups(wildcards):
    return ' '.join(wildcards.comparison.split('_vs_'))

def lookup_tables_from_gtf(gtf):
    exon_pos = {}  # (ENSG, Exon number) to IV.
    transl = collections.defaultdict(set)  # ENSG to gene_id

    # Set these lookup tables from the GTF.
    with open(gtf) as f:
        for li in f:
            if li[0] == '#':
                continue
            ensg = re.search('gene_id "([^"]+)"', li)
            gene_name = re.search('gene_name "([^"]+)"', li)
            if ensg is not None and gene_name is not None:
                transl[ensg.group(1)].add(gene_name.group(1))

            exon_n = re.search('exon_number "([^"]+)"', li)
            if ensg is not None and exon_n is not None:
                s = li.split("\t")
                exon_pos[(ensg.group(1), exon_n.group(1))] = f"{s[0]}:{s[3]}-{s[4]}/{s[6]}"

    transl = {k:list(v)[0] for k,v in transl.items()}

    return transl, exon_pos

rule combine_dexseq_results:
    input:
        dex_results = expand("data/processed/{comparison}/dexseq/DEXSeq_results.annotated.txt", comparison=COMPARISONS)
    output:
        xlsx = "outs/dexseq/dexseq.combined_table.xlsx",
    run:
        writer = pandas.ExcelWriter(str(output.xlsx), engine='openpyxl')
        for fname in input.dex_results:
            df = pandas.read_csv(fname, sep='\t')
            df = df[df['padj']<0.2]
            print("Writing to a sheet: ", re.search('data/processed/(.+)/dexseq', fname).group(1))
            df.to_excel(writer, sheet_name=re.search('data/processed/(.+)/dexseq', fname).group(1))
        writer.save()        
    
rule annotate_dexseq_results:
    input:
        gtf = "data/processed/{comparison}/genome.filtered_by_rsem.fc.gtf",
        dex_results = "data/processed/{comparison}/dexseq/DEXSeq_results.txt", 
    output:
        dex_results = "data/processed/{comparison}/dexseq/DEXSeq_results.annotated.txt", 
    run:
        
        # Look up tables.
        transl, exon_pos = lookup_tables_from_gtf(str(input.gtf))
        
#        for fname in [str(x) for x in input.dex_results]: #glob.glob(f"{top}/outs/R/dexseq/*/DEXSeq_results.txt"):
#            print(fname)
#            """ 
        df = pandas.read_csv(str(input.dex_results), sep='\t')
        print(df.head())
        #df["ENSG"] = [recover_ensg(str(x)) for x in df.groupID]

        df['gene_name'] = [transl.get(x, x) for x in df['groupID']]
        df = df.sort_values(by='padj')

        df['exon_pos'] = [exon_pos.get((a,b.lstrip('E')), "") for a,b in zip(df.groupID, df.featureID)]
        df['iv'] = [x.split('/')[0] for x in df.exon_pos]
        if 'Unnamed: 0' in df.columns:
            del df['Unnamed: 0']
        print(df.head())

        df.to_csv(str(output.dex_results), sep='\t', index=False)

        #"""
        
rule run_dexseq:
    """
    COUNTS = args[1]  # subread_fc.out counts file.
SAMPLE_TABLE = args[2]  # column of sample names match bam files and column of groups.
GTF = args[3]  # Can't just be a GTF. Has to be reformatted.
GROUP_A = args[4]  # Should match the "group" values in SAMPLE_TABLE.
GROUP_B = args[5]

    Sample table is expected to be of this format:
sample  group
Empty_pLEX-anti-PCBP1_shRNA-Batch_1-HCT116_cells-SRR13169561    HCT116_EV_pLEX_and_shPCBP1
Empty_pLEX-anti-PCBP1_shRNA-Batch_1-HCT116_cells-SRR13169562    HCT116_EV_pLEX_and_shPCBP1
Empty_pLEX-control_shRNA-Batch_1-HCT116_cells-SRR13169563       HCT116_EV_pLEX_and_control_shRNA

    """
    input:
        subread_counts = "data/processed/{comparison}/subread_fc.out",
        sample_table = SAMPLE_TABLE,
        gtf = "data/processed/{comparison}/genome.filtered_by_rsem.fc.gtf",
    output:
        subread_counts_simplified_names = 'data/processed/{comparison}/subread_fc.out.simplified_col_names',
        dex_results = "data/processed/{comparison}/dexseq/DEXSeq_results.txt",
    conda:
        "env/bioconda.yaml"
    params:
        groups = space_sep_groups
    shell:
        "Rscript scripts/DEXSeq_of_exons.generalized.R {input.subread_counts} {input.sample_table}" + \
        " {input.gtf} {params.groups}"

########################################################################
# Before DEXSeq.
########################################################################

rule count_exons_with_Rsubread:
    """
    Use the featureCounts() function in Rsubread on bam files to generate exons by samples counts matrices.
    
    Following https://github.com/vivekbhr/Subread_to_DEXSeq:
    
    1. Edit annotation (make a new gtf and gff).
    # Input: data/processed/Homo_sapiens.GRCh38.104.chr.filtered.gtf
python Subread_to_DEXSeq/dexseq_prepare_annotation2.py --aggregate no --featurecountsgtf  data/processed/Homo_sapiens.GRCh38.104.chr.filtered.fc.gtf data/processed/Homo_sapiens.GRCh38.104.chr.filtered.gtf  data/processed/Homo_sapiens.GRCh38.104.chr.filtered.fc.gff

    2. Count reads in exons using the new gtf and bam files:
    sbatch call_feature_counts.sbatch
    Which runs:

#ml load star/2.7.8a
#ml load rsem/1.3.3
ml subread

a=(STAR_genome_hg38/W*/*.sortedByCoord.out.bam)
bams=${a[@]}

# -f use features (exons). -T # CPUs.
featureCounts -f -O -s 2 -p -T 16 -F GTF \
-a data/processed/Homo_sapiens.GRCh38.104.chr.filtered.fc.gtf \
-o data/processed/subread_fc.out $bams

    """
    input:
        bams = wc_bam_files_for_comparison,
        gtf = "data/processed/{comparison}/genome.filtered_by_rsem.fc.gtf"
    output:
        #counts = expand("outs/exon_counts/{sample}", sample=samples),
        counts = "data/processed/{comparison}/subread_fc.out",
    conda:
        "env/bioconda.yaml"
    threads:
        16
    params:
        bam_files = space_delim_bam_files_for_comparison
    shell:
        "featureCounts -f -O -s 2 -p -T 16 -F GTF -a {input.gtf} -o {output.counts} " + \
        "{params.bam_files}"
        
rule format_gtf_for_Rsubread:
    """
    Has two outputs: a featurecounts gtf and a gff.
    python %prog [options] <in.gtf> <out.gff>
    
    python Subread_to_DEXSeq/dexseq_prepare_annotation2.py --aggregate no --featurecountsgtf  data/processed/Homo_sapiens.GRCh38.104.chr.filtered.fc.rescue.gtf data/processed/Homo_sapiens.GRCh38.104.chr.filtered.rescue.gtf data/processed/Homo_sapiens.GRCh38.104.chr.filtered.rescue.gtf
    """
    input:
        gtf = "data/processed/{comparison}/genome.filtered_by_rsem.gtf",
    output:
        gtf = "data/processed/{comparison}/genome.filtered_by_rsem.fc.gtf",
        gff = "data/processed/{comparison}/genome.filtered_by_rsem.fc.gff"
    shell:
        "python Subread_to_DEXSeq/dexseq_prepare_annotation2.py --aggregate no --featurecountsgtf {output.gtf} {input.gtf} {output.gff}"

rule prefilter_isoforms_using_rsem_results:
    """
    https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0862-3
    =>5% cutoff. Removes ~60% of the exons and isoforms.
    
    python scripts/filter_isoforms_using_rsem_results.py ref/Homo_sapiens.GRCh38.104.chr.gtf  rescue_RSEM_hg38/ data/processed/Homo_sapiens.GRCh38.104.chr.filtered.rescue.gtf
    
    python scripts/filter_isoforms_using_rsem_results.py ref/Homo_sapiens.GRCh38.104.chr.gtf 211014_rescue/RSEM_hg38/ 211014_rescue/data/processed/Homo_sapiens.GRCh38.104.chr.filtered.gtf
    
    # From that command got this output:
Started with 180929 across 60605 genes. Removed those with less than 5.00% of TPM at that ENSG. Kept 57920 (32.01%).
Filtering gtf: kept 1013834 (32.23%) lines out of 3145418 lines. Discarded 2070979 as having an ENST below cutoff, and 60605 for having no ENST in the GTF.

    # For the original dataset, get this output:
    Started with 185270 across 60605 genes. Removed those with less than 5.00% of TPM at that ENSG. Kept 62107 (33.52%).
Filtering gtf: kept 1040567 (33.08%) lines out of 3145418 lines. Discarded 2044246 as having an ENST below cutoff, and 60605 for having no ENST in the GTF.
    """

    input:
        counts = wc_counts_files_for_comparison,
        gtf = GTF,
    output:
        gtf = "data/processed/{comparison}/genome.filtered_by_rsem.gtf",
    run:
        shell("python scripts/filter_isoforms_using_rsem_results.py --gtf {input.gtf} --rsem_dir " + RSEM_COUNTS_DIR + "  --out {output.gtf}  --rsem_files " + comma_delim_counts_files_for_comparison(wildcards.comparison))

rule count_with_rsem:
    """
    To make a bunch of input bams in a single directory, run a command like:
    mkdir STAR_genome_hg38/transcriptomes
    for f in STAR_genome_hg38/*/*toTranscriptome.out.bam; do base=$(basename ${f} _Aligned.toTranscriptome.out.bam); out=STAR_genome_hg38/transcriptomes/$(basename ${f}) ; ln -s $(realpath ${f}) $out ; done
    """    
    input:
        rsem_ref = RSEM_REF + '.seq',
        bam = BAM_DIR + "/{sample}_Aligned.toTranscriptome.out.bam"
    output:
        counts = RSEM_COUNTS_DIR + "/{sample}.isoforms.results",
    params:
        counts_prefix = RSEM_COUNTS_DIR + "/{sample}"
    threads:
        8
    conda:
        "env/rsem.yaml",
    shell:
        #"mkdir -p {output.counts_dir} ; " + \
        "rsem-calculate-expression --bam --paired-end -p 8 {input.bam} {RSEM_REF} {params.counts_prefix}"
        
rule build_reference:
    """
        ml load star/2.7.8a
        ml load rsem/1.3.3
    """
    input:
        gtf = GTF, #"ref/Homo_sapiens.GRCh38.104.chr.gtf.gz",
        fa = GENOME_FA, #"ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    output:
        seq = RSEM_REF + ".seq"
    params:
        prefix = RSEM_REF
    conda:
        "env/rsem.yaml",
    shell:
        "rsem-prepare-reference --gtf {input.gtf} {input.fa} {params.prefix}"
        # The --star parameter is necessary. It builds an index for STAR.
        
rule download_reference:
    """
    Don't actually run this rule!
    """
    run:
        shell("mkdir ref/")
        shell("cd ref/")
        shell("wget http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.chr.gtf.gz")
        shell(
    "wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz")
        shell("gunzip *.gz")
        shell("cd ..")