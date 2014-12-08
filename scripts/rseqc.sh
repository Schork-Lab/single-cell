
#!/bin/sh

# BAM file that you want to do the analysis on
BAM=$1

# Output directory where you want all the results to go
# Files will be prepended with "$OUT_DIR/rseqc"
OUT_DIR=$2

# BED files of which you want to see the gene body coverage
BED=$3

GENE_BEDS=${@: 4}

# Make plots via PREFIX
PREFIX=$OUT_DIR/rseqc
RSeQCDIR=/home/kbhutani/tools/rna_seq/RSeQC-2.4/build/scripts-2.7
PYTHON=/opt/python/bin/python2.7

# Extract only mapped reads, because of interference with RSeqC
NEWBAM=${BAM##*/}
NEWBAM=${NEWBAM%.bam}
NEWBAM="$OUT_DIR/$NEWBAM.mapped.bam"
samtools view -h -F 4 $BAM | awk '{if($3 !~ "ERCC"){print $0}}' | samtools view -Sb - > $NEWBAM
samtools index $NEWBAM
BAM=$NEWBAM

# Added for this use case because the bed derived from the gtf that was used to align the reads is 
# not working correctly.
BED='/projects/ps-jcvi/tools/references/rnaseq/hg19_RefSeq.bed'
CONTROL='/projects/ps-jcvi/projects/Lasken_Single_Cell/data/nov_24_run/HBPC_RSEM_RNASeq_PV_12012014/BAM_files/genome/3DHCN100pg111314.genome.sorted.bam'

module load scipy/2.7
module load R
echo 'executing: clipping_profile.py' ; $PYTHON $RSeQCDIR/clipping_profile.py --input-file $BAM --out-prefix $PREFIX

for gene_bed in $GENE_BEDS
    do 
        echo "executing: geneBody_coverage.py for $gene_bed"
        s=${gene_bed##*/}
        out_prefix="$PREFIX-${s%.bed}"
        echo $out_prefix
        $PYTHON $RSeQCDIR/geneBody_coverage.py -i $BAM,$CONTROL --refgene $gene_bed --out-prefix $out_prefix
    done
echo 'executing: infer_experiment.py' ; $PYTHON $RSeQCDIR/infer_experiment.py --input-file $BAM --refgene $BED >$PREFIX.infer_experiment.txt
echo 'executing: inner_distance.py' ; $PYTHON $RSeQCDIR/inner_distance.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'executing: junction_annotation.py' ; $PYTHON $RSeQCDIR/junction_annotation.py --input-file $BAM --refgene $BED --out-prefix $PREFIX >$PREFIX.junction_annotation.txt
echo 'executing: junction_saturation.py' ; $PYTHON $RSeQCDIR/junction_saturation.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'executing: read_distribution.py' ; $PYTHON $RSeQCDIR/read_distribution.py --input-file $BAM --refgene $BED > $PREFIX.read_distribution.txt
echo 'executing: read_duplication.py' ; $PYTHON $RSeQCDIR/read_duplication.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: read_GC.py' ; $PYTHON $RSeQCDIR/read_GC.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: read_NVC.py' ; $PYTHON $RSeQCDIR/read_NVC.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: read_quality.py' ; $PYTHON $RSeQCDIR/read_quality.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: RPKM_count.py' ; $PYTHON $RSeQCDIR/RPKM_count.py --input-file $BAM --out-prefix $PREFIX --refgene $BED
echo 'executing: RPKM_saturation.py' ; $PYTHON $RSeQCDIR/RPKM_saturation.py --input-file $BAM --out-prefix $PREFIX --refgene $BED
