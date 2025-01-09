#!/bin/bash
#SBATCH --partition=employee
#SBATCH --array=1
#SBATCH --cpus-per-task=10
#SBATCH --time=10:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --nodelist=fgcz-c-041
#SBATCH --job-name=SQUANTIQC


# Author: Natalia Zajac

#source /usr/local/ngseq/miniconda3/etc/profile.d/conda.sh
eval "$(/usr/local/ngseq/miniforge3/bin/conda shell.bash hook)"
conda activate gi_SQANTI3.5.2.2
#export PYTHONPATH=$PYTHONPATH:/usr/local/ngseq/src/cDNA_Cupcake/sequence/

k=$SLURM_ARRAY_TASK_ID
name=`sed -n ${k}p < /srv/GT/analysis/zajacn/p28443/SQANTI3_QC/list_of_samples`

REF=/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/sl_unzip_datasets/5c640510-4bee-4ae8-b642-48ed209b79bb/call-unzip_datasets/execution/human_GRCh38_no_alt_analysis_set.fasta
ANN=/misc/sequel2/SMRTLink_v10_dataOutput/cromwell-executions/sl_unzip_datasets/5c640510-4bee-4ae8-b642-48ed209b79bb/call-unzip_datasets/execution/gencode.v39.annotation.sorted.gtf
DIR=/srv/GT/analysis/zajacn/p28443/SQANTI3_QC

python /usr/local/ngseq/src/SQANTI3-5.2.2/sqanti3_qc.py --SR_bam $DIR/${name}/SR_BAM/ -d $DIR/${name} --skipORF --polyA_motif_list /usr/local/ngseq/src/SQANTI3-5.2.2/data/polyA_motifs/mouse_and_human.polyA_motif.txt --cpus 10 --report pdf $DIR/${name}/scisoseq_transcripts.sorted.filtered_lite.gff $ANN $REF

