RNA_N=rna002_r9_ec_native
RNA_I=rna002_r9_ec_ivt
DNA_R9=dna_r9_dm_ctl
DNA_R10=dna_r10_dm_ctl


if [ ! -d raw ] || [ ! -d ref ] || [ ! -d mm2 ]; then
    echo "Data not found. Please run 'download.sh' first."
    exit 1
fi

#DATASETS=$RNA_N $RNA_N $DNA_r9 $DNA_r9
mkdir -p out

check () {
    name=$1
    cmd=$2 
    out=out/$name.out
    log=out/$name.err
    echo "Running '$cmd 2> $log'"
    `$cmd 2> $log`
    ret=$?
    if [ $ret -ne 0 ]; then
        echo "$name: failed with code $ret"
        echo "$ cat $log"
        cat $log
        #exit $ret
    fi
}

for NAME in $RNA_I $RNA_N; do
    out=${NAME}.align
    check "$out" "uncalled4 align --ref ref/ec_16s.fa --reads raw/${NAME}.fast5 --bam-in mm2/${NAME}.bam --bam-out out/${out}.bam"
done

for NAME in $DNA_R10 $DNA_R9; do
    out=${NAME}.align
    check "$out" "uncalled4 align --ref ref/dm_chr1.fa --reads raw/${NAME}.fast5 --bam-in mm2/${NAME}.bam --bam-out out/${out}.bam"
done

NAME=${RNA_I}
for FMT in pod5 blow5; do
    out="${NAME}.${FMT}.align"
    check $out "uncalled4 align --flowcell FLO-MIN106 --kit SQK-RNA002 --ref ref/ec_16s.fa --reads raw/${NAME}.${FMT} --bam-in mm2/${NAME}.bam --bam-out out/${out}.bam"
done

for NAME in $RNA_I $RNA_N $DNA_R10 $DNA_R9 $RNA_I.blow5 $RNA_I.pod5; do
    check "$NAME.align.index" "samtools index out/$NAME.align.bam"
done

check "rna.refstats" "uncalled4 refstats dtw.current ks,mean --bam-in out/$RNA_I.align.bam out/$RNA_N.align.bam -o out/rna.refstats.tsv"
check "rna.compare" "uncalled4 compare --bam-in out/rna002_r9_ec_ivt.align.bam out/rna002_r9_ec_ivt.pod5.align.bam -o out/rna.comare.tsv"

check "$DNA_R9.train.k6" "uncalled4 train --ref ref/dm_chr1.fa --reads raw/$DNA_R9.fast5 --bam-in mm2/$DNA_R9.bam -k 6 --kmer-shift 3 --out-dir out/$DNA_R9.train.k6 --train-iterations 2 --init-mode moves > out/train.txt"

check "$RNA_I.dotplot" "uncalled4 dotplot --bam-in out/$RNA_I.align.bam -l e679ebd4-2dc5-4fe2-95b0-6ec80df5dc95 -L dtw.current -o out/$RNA_I.dotplot."
check "$RNA_I.trackplot" "uncalled4 trackplot --region ecoli16S:0-1000 --bam-in out/$RNA_I.align.bam --mat dtw.model_diff -o out/$RNA_I.trackplot.html"

echo "Running uncalled4 browser ecoli16S out/$RNA_I.align.bam out/$RNA_N.align.bam"
uncalled4 browser -R ecoli16S --bam-in out/$RNA_I.align.bam out/$RNA_N.align.bam
