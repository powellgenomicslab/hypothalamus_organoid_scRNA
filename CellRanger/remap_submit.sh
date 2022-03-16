#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
#          File: main.sge
#                Nextflow script for Cellranger pipeline
#   Modified by: vikgna and Drew neavin
#   Modified on: 2020/02/03
#       Version: 1.1.3
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/


for pool in 1205_GEX_WT_THpos_GEX 1207_GEX_WT_rest_GEX
do

	echo $pool

	Sample_Id=$pool\_0_1_HC3GJDSXY # define sample name
	PIPELINE="/path/to//cellranger.sh"
	OUTPUT_DIR="/path/to/outdir/" # define output path
	LOGS="$OUTPUT_DIR/logs"
	INPUT_DIR=/path/to/fastq/$pool/$Sample_Id # define fastq path
	ref="/path/to/refdata-gex-GRCh38-2020-A_w_tdtomato" # define reference path


	mkdir -p $LOGS


	qsub -S /bin/bash \
		-q short.q \
		-r yes \
		-l mem_requested=8G \
		-l tmp_requested=256G \
		-pe smp 16 \
		-N remapping \
		-cwd \
		-j y \
		-e $LOGS \
		-o $LOGS \
		-V \
		-v pool=$pool,OUTPUT_DIR=$OUTPUT_DIR,INPUT_DIR=$INPUT_DIR,Sample_Id=$Sample_Id,ref=$ref \
		-C '' $PIPELINE

done
