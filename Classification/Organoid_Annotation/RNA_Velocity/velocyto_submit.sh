#!/bin/bash


SCRIPT="/path/to/script/velocyto.sh" ### available on the hypothalamus_organoid_scRNA github (https://github.com/powellgenomicslab/hypothalamus_organoid_scRNA)
TENxDIR="/path/to/tenx_dir/" ### the output form the cellranger.sh script - available on the hypothalamus_organoid_scRNA github (https://github.com/powellgenomicslab/hypothalamus_organoid_scRNA)
OUT="/path/to/velocyto/outdir/"
GTF="/path/to/tenx/genes.gtf"
LOG=$OUT/logs

mkdir -p $LOG

T=32

### Generate a file that can be used downstream to easily call the results
echo "Pool\tName\tDirectory" > $OUT/velocyto_files.tsv


for pool in 1205_GEX_WT_THpos_GEX_0_1_HC3GJDSXY 1207_GEX_WT_rest_GEX_0_1_HC3GJDSXY
do
	echo $pool
    ### Generate a file that can be used downstream to easily call the results
    name=`sed 's/1205_GEX_//g' $pool | sed 's/1207_GEX_//g' | sed 's/_GEX_0_1_HC3GJDSXY//g'`
    echo "$pool\t$name\t$TENxDIR/velocyto/$pool.loom" >> $OUT/velocyto_files.tsv

	qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=4G \
        -l tmp_requested=4G \
		-pe smp $T \
        -N velocyto_$pool \
        -cwd \
        -m e \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v TENxDIR=$TENxDIR,OUT=$OUT,pool=$pool,GTF=$GTF,T=$T \
        -C '' $SCRIPT
done