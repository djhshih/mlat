#!/bin/sh

mkdir -p out

./mlat

# fail
./mlat database query out/out-fail.psl || true

small="-tileSize=6 -stepSize=2 -minScore=0"
./mlat $small data/ref1.fna data/query1.fna out/out-0.0.psl
./mlat $small -mask=lower -qMask=lower data/ref1.fna data/query1.fna out/out-0.1.psl
./mlat $small -trimT -trimHardA data/ref1.fna data/query1.fna out/out-0.2.psl
./mlat $small data/ref1.2bit:seq1:5-85 data/query1.fna out/out-0.3.psl
./mlat $small data/ref1.2bit:seq1:5-85 data/query1.2bit:read1.2 out/out-0.4.psl

# fail
./mlat -t=prot -q=prot data/ref1.2bit:seq1:5-85 data/query1.2bit out/out-fail.psl || true

./mlat -makeOoc=11.ooc data/hg38_tp53.2bit data/reads_tp53.fa out/out-1.0.psl
./mlat data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.psl
./mlat -ooc=11.ooc data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.psl

formats=( pslx axt maf sim4 wublast blast blast8 blast9 )
for format in ${formats[@]}; do
	./mlat -out=${format} data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.${format}
done

./mlat -t=dnax -q=dnax data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-2.0.psl

./mlat -t=prot -q=prot data/protein_tp53.faa data/peptides_tp53.faa out/out-3.0.psl

./mlat -t=dnax -q=prot data/hg38_tp53.2bit data/peptides_tp53.faa out/out-4.0.psl


echo "Checking outputs..."

cases=( 0.0 0.1 0.2 0.3 0.4 1.0 1.1 2.0 3.0 4.0 )

for x in ${cases[@]}; do
	echo "diff out/out-${x}.psl test/test-${x}.psl"
	diff out/out-${x}.psl test/test-${x}.psl
done

