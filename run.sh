#!/bin/sh

./mlat

./mlat database query test.psl || true

./mlat data/ref1.fna data/query1.fna test0.psl
./mlat -mask=lower -qMask=lower data/ref1.fna data/query1.fna test0.psl
./mlat -trimT -trimHardA data/ref1.fna data/query1.fna test0.psl
./mlat data/ref1.2bit:seq1:5-85 data/query1.fna test0.psl

./mlat data/hg38_tp53.2bit data/reads_tp53.fa test1.psl
./mlat data/hg38_tp53.2bit data/reads_tp53.fa.gz test1.psl

./mlat -t=dnax -q=dnax data/hg38_tp53.2bit data/reads_tp53.fa.gz test2.psl

./mlat -t=prot -q=prot data/protein_tp53.faa data/peptides_tp53.faa test3.psl

