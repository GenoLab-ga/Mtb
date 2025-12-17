#!/bin/bash



# liste tous les .sra dans les sous-dossiers
for sra_file in raw/SRR*/SRR*.sra; do
    SRR=$(basename "$sra_file" .sra)
    echo "Conversion de $SRR ..."

    # extraction en fastq
    fasterq-dump "$sra_file" --split-files -O raw_data

    # compression
    gzip raw_data/${SRR}_1.fastq
    gzip raw_data/${SRR}_2.fastq
done
