#!/bin/sh

./mlat

# fail
./mlat database query test.psl || true

small="-tileSize=6 -stepSize=2 -minScore=0"
./mlat $small data/ref1.fna data/query1.fna test-0.0.psl
./mlat $small -mask=lower -qMask=lower data/ref1.fna data/query1.fna test-0.1.psl
./mlat $small -trimT -trimHardA data/ref1.fna data/query1.fna test-0.2.psl
./mlat $small data/ref1.2bit:seq1:5-85 data/query1.fna test-0.3.psl
./mlat $small data/ref1.2bit:seq1:5-85 data/query1.2bit:read1.2 test-0.4.psl

# fail
./mlat -t=prot -q=prot data/ref1.2bit:seq1:5-85 data/query1.2bit test-0.5.psl || true

./mlat -makeOoc=11.ooc data/hg38_tp53.2bit data/reads_tp53.fa test-1.0.psl
./mlat data/hg38_tp53.2bit data/reads_tp53.fa.gz test-1.1.psl
./mlat -ooc=11.ooc data/hg38_tp53.2bit data/reads_tp53.fa.gz test-1.1.psl

formats=( pslx axt maf sim4 wublast blast blast8 blast9 )
for format in ${formats[@]}; do
	./mlat -out=${format} data/hg38_tp53.2bit data/reads_tp53.fa.gz test-1.1.${format}
done

./mlat -t=dnax -q=dnax data/hg38_tp53.2bit data/reads_tp53.fa.gz test-2.0.psl

./mlat -t=prot -q=prot data/protein_tp53.faa data/peptides_tp53.faa test-3.0.psl

./mlat -t=dnax -q=prot data/hg38_tp53.2bit data/peptides_tp53.faa test-4.0.psl

