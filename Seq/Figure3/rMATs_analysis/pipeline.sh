######################################################################
# 1. Mapping: mapWithStar.sbatch
######################################################################
ml load star/2.7.8a

for fq1 in raw_data/*/*_1.fq.gz ; 

do base=$(basename ${fq1} _1.fq.gz); fq2=${fq1%_1.fq.gz}_2.fq.gz ; 
mkdir -p star/${base}/ ; 

STAR --runThreadN 12 --quantMode TranscriptomeSAM --genomeDir /oak/stanford/groups/khavari/users/dfporter/seq/genomes/star_gencode_39/ --readFilesIn ${fq1} ${fq2} --readFilesCommand zcat --outFileNamePrefix star/${base}/${base}_ --outSAMtype BAM Unsorted SortedByCoordinate --sjdbGTFfile /oak/stanford/groups/khavari/users/dfporter/seq/genomes/gencode.v39.annotation.gtf ; 

done

######################################################################
# 2. Quantify expression from bams: count_with_rsem.sh
######################################################################

ml load star/2.7.8a
ml load rsem/1.3.3

top=$(realpath . )

mkdir RSEM_hg38/

for bam in star/*/*_Aligned.toTranscriptome.out.bam;
do

bam=$(realpath ${bam})
base=$(basename ${bam} _Aligned.toTranscriptome.out.bam)
out=RSEM_hg38/

mkdir -p ${out}${base}
cd ${out}${base}

# Usage:
#rsem-calculate-expression [options] --alignments [--paired-end] input reference_name sample_name
#p=number of threads.

sbatch --nodes=1 -p apps_khavari --cpus-per-task=8 --qos=normal -t 8:00:00 --mem=20000 --wrap="module load rsem/1.3.3
rsem-calculate-expression --bam --paired-end -p 8 ${bam} /storage/khavari/data/dfporter/ddx21/rna_seq_splicing_other_factors_kd/rsem_reference/GRCh38.gencode39.all.chrom.reference ${base}"

cd ${top}

done

######################################################################
# 3. Make symlinks for the snakemake:
# prepareRsemCountResultsForSnakemake.sh
######################################################################

mkdir -p data/rsem_counts
for f in RSEM_hg38/*/*bam ; do ln -s $(realpath ${f}) data/rsem_counts/$(basename ${f}) ; done
for f in RSEM_hg38/*/*genes.results ; do ln -s $(realpath ${f}) data/rsem_counts/$(basename ${f}) ; done
for f in RSEM_hg38/*/*isoforms.results ; do ln -s $(realpath ${f}) data/rsem_counts/$(basename ${f}) ; done
mkdir star/transcriptomes ; for f in RSEM_hg38/*/*bam ; do ln -s $(realpath ${f}) star/transcriptomes/$(basename ${f}) ; done

######################################################################
# 4. Run the snakemake.
######################################################################

snakemake -s snake.mapped_reads_to_dexseq_and_rmats.py --configfile config.yaml --use-conda -j 1 

######################################################################
# 5. Call the rmats scripts output by the snakemake: callRmats.sh 
######################################################################
for f in data/processed/rmats_inputs/*sbatch;
do

#sbatch -w smsh11dsu-srcf-d15-38  --nodes=1 -p apps_khavari --cpus-per-task=8 --qos=normal -t 8:00:00 --mem=20000 --wrap="
sbatch --nodes=1 -p interactive --cpus-per-task=1 --qos=normal -t 14:00:00 --mem=20000 --wrap="
bash ${f}"

done
