bpipe \
    -p caller=$CLINKERDIR/test/caller/bcr_abl1.csv \
    -p col=1,2,3,4 \
    -p genome=19 \
    -p print=true \
    -p competitive=false \
    -p header=true \
    -p align_mem=1000000000 \
    -p genome_mem=360000000000 \
    -p fusions=BCR:ABL1 \
    -p out=$CLINKERDIR/test/results \
    $CLINKERDIR/workflow/clinker.pipe \
    $CLINKERDIR/test/fastq/*.fastq.gz
