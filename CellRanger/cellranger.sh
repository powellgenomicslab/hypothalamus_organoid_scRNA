
# copy fastq files to local scratch
echo $TMPDIR
cd  "$TMPDIR"
# copy fastq files into local scratch
mkdir -p "$TMPDIR/fastq"
cp -rp "${INPUT_DIR}" "$TMPDIR/fastq"
# copy reference 
cp -r ${ref} "$TMPDIR/reference"

/directflow/SCCGGroupShare/projects/DrewNeavin/tools/cellranger-6.1.1/cellranger count --id="${Sample_Id}" \
                       --fastqs="$TMPDIR/fastq" \
                       --transcriptome="$TMPDIR/reference" \
                       --sample=${Sample_Id}\
                       --jobmode=local \
                       --localcores=16 \
                       --localmem=256
