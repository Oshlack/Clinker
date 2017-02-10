samtools index results/alignment/Aligned.sortedByCoord.out.bam
bamCoverage -b results/alignment/Aligned.sortedByCoord.out.bam --normalizeUsingRPKM -of bedgraph --binSize 1 -o results/alignment/coverage_rpm.bedgraph
../plotit/fst_plot_prep.sh BCR:ABL1 ../test/results/alignment/BCR_ABL1 BCR_ABL1 ../test/results/annotation ../test/results/alignment ../test/results/reference