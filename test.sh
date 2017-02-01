#!/bin/sh

set -o nounset
set -o pipefail


path=build/bin

mkdir -p out

$path/mlat || true

# fail
$path/mlat database query out/out-fail.psl || true

small="-tileSize=6 -stepSize=2 -minScore=1"
$path/mlat $small data/ref1.fna data/query1.fna out/out-0.0.psl
$path/mlat $small -mask=lower -qMask=lower data/ref1.fna data/query1.fna out/out-0.1.psl
$path/mlat $small -trimT -trimHardA data/ref1.fna data/query1.fna out/out-0.2.psl
$path/mlat $small data/ref1.2bit:seq1:5-85 data/query1.fna out/out-0.3.psl
$path/mlat $small data/ref1.2bit:seq1:5-85 data/query1.2bit:read1.2 out/out-0.4.psl

# fail
$path/mlat -t=prot -q=prot data/ref1.2bit:seq1:5-85 data/query1.2bit out/out-fail.psl || true

$path/mlat -makeOoc=out/11.ooc data/hg38_tp53.2bit data/reads_tp53.fa out/out-1.0.psl
$path/mlat data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.psl
$path/mlat -ooc=out/11.ooc data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.psl

formats=( pslx axt maf sim4 wublast blast blast8 blast9 )
for format in ${formats[@]}; do
	$path/mlat -out=${format} data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-1.1.${format}
done

$path/mlat -t=dnax -q=dnax data/hg38_tp53.2bit data/reads_tp53.fa.gz out/out-2.0.psl

$path/mlat -t=prot -q=prot data/protein_tp53.faa data/peptides_tp53.faa out/out-3.0.psl

$path/mlat -t=dnax -q=prot data/hg38_tp53.2bit data/peptides_tp53.faa out/out-4.0.psl

cd demo
./mlat-demo ../data/ref1.fna AGACGGTCGATCGGGATTCGAGGTCGA > ../out/out-5.0.tsv
./mlat-demo-cpp ../data/ref1.fna AGACGGTCGATCGGGATTCGAGGTCGA > ../out/out-5.1.tsv

./mlat-demo ../data/hg38_tp53.2bit TCATGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGATTCCATCTCAAAAAAAAAAAAAAGGCCTCCCCTGCTT > ../out/out-6.0.tsv
./mlat-demo-cpp ../data/hg38_tp53.2bit TCATGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGATTCCATCTCAAAAAAAAAAAAAAGGCCTCCCCTGCTT > ../out/out-6.1.tsv
cd ..


echo "Checking outputs..."

cases=( \
	0.0.psl 0.1.psl 0.2.psl 0.3.psl 0.4.psl \
	1.0.psl 1.1.psl \
	1.1.pslx 1.1.axt 1.1.maf 1.1.sim4 1.1.wublast 1.1.blast 1.1.blast8 1.1.blast9 \
	2.0.psl 3.0.psl 4.0.psl \
	5.0.tsv 5.1.tsv 6.0.tsv 6.1.tsv \
)

for x in ${cases[@]}; do
	echo "diff out/out-${x} test/test-${x}"
	diff out/out-${x} test/test-${x}
done

