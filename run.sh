#!/bin/sh

./mlat

./mlat database query test.psl || true

./mlat data/hg38_tp53.2bit data/reads_tp53.fa test.psl
./mlat data/hg38_tp53.2bit data/reads_tp53.fa.gz test.psl

./mlat -t=prot -q=prot data/protein_tp53.faa data/peptides_tp53.faa test2.psl

