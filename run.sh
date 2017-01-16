#!/bin/sh

./mlat

./mlat database query test.psl || true

./mlat data/ref1.fna data/query1.fna test0.psl
./mlat -mask=lower data/ref1.fna data/query1.fna test0.psl
./mlat -trimT -trimHardA data/ref1.fna data/query1.fna test0.psl

./mlat data/hg38_tp53.2bit data/reads_tp53.fa test1.psl
./mlat data/hg38_tp53.2bit data/reads_tp53.fa.gz test1.psl

./mlat -t=prot -q=prot data/protein_tp53.faa data/peptides_tp53.faa test2.psl

