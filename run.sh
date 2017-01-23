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

bin/mlat-example data/ref1.fna AGACGGTCGATCGGGATTCGAGGTCGA > out/out-5.0.tsv


echo "Checking outputs..."

cases=( \
	0.0.psl 0.1.psl 0.2.psl 0.3.psl 0.4.psl \
	1.0.psl 1.1.psl \
	1.1.pslx 1.1.axt 1.1.maf 1.1.sim4 1.1.wublast 1.1.blast 1.1.blast8 1.1.blast9 \
	2.0.psl 3.0.psl 4.0.psl 5.0.tsv \
)

for x in ${cases[@]}; do
	echo "diff out/out-${x} test/test-${x}"
	diff out/out-${x} test/test-${x}
done

