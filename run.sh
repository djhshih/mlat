#!/bin/sh

mkdir -p out

bin/mlat

# fail
bin/mlat database query out/out-fail.psl || true

small="-tileSize=6 -stepSize=2 -minScore=0"
bin/mlat $small data/ref1.fna data/query1.fna out/out-0.0.psl
bin/mlat $small -mask=lower -qMask=lower data/ref1.fna data/query1.fna out/out-0.1.psl
bin/mlat $small -trimT -trimHardA data/ref1.fna data/query1.fna out/out-0.2.psl
bin/mlat $small data/ref1.2bit:seq1:5-85 data/query1.fna out/out-0.3.psl
bin/mlat $small data/ref1.2bit:seq1:5-85 data/query1.2bit:read1.2 out/out-0.4.psl

# fail
bin/mlat -t=prot -q=prot data/ref1.2bit:seq1:5-85 data/query1.2bit out/out-fail.psl || true

bin/mlat -makeOoc=out/11.ooc data/hg38_tp53.2bit data/reads_tp53.fa out/out-1.0.psl
bin/mlat data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.psl
bin/mlat -ooc=out/11.ooc data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.psl

formats=( pslx axt maf sim4 wublast blast blast8 blast9 )
for format in ${formats[@]}; do
	bin/mlat -out=${format} data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.${format}
done

bin/mlat -t=dnax -q=dnax data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-2.0.psl
../blat-lite/bin/blat -t=dnax -q=dnax data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-2.0b.psl

bin/mlat -t=prot -q=prot data/protein_tp53.faa data/peptides_tp53.faa out/out-3.0.psl

bin/mlat -t=dnax -q=prot data/hg38_tp53.2bit data/peptides_tp53.faa out/out-4.0.psl

bin/mlat-alt data/ref1.fna AGACGGTCGATCGGGATTCGAGGTCGA > out/out-5.0.tsv


echo "Checking outputs..."

cases=( 0.0 0.1 0.2 0.3 0.4 1.0 1.1 2.0 3.0 4.0 5.0 )

for x in ${cases[@]}; do
	echo "diff out/out-${x}.psl test/test-${x}.psl"
	diff out/out-${x}.psl test/test-${x}.psl
done

