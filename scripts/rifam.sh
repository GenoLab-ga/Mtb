#!/bin/bash

# Liste des identifiants
srr_list=(
SRR13687743 SRR13687744 SRR13687745 SRR13687746 SRR13687747 SRR13687748 SRR13687749 SRR13687750 SRR13687751 # Témoins sans rifampicine

SRR13687752 SRR13687753 SRR13687754 SRR13687755 SRR13687756 SRR13687757 SRR13687758 SRR13687759 SRR13687760 # Traités à 0.008µg/ml de rifampicine

SRR13687761 SRR13687762 SRR13687763 SRR13687764 SRR13687765 SRR13687766 SRR13687767 SRR13687768 SRR13687769 # Traités à 0.002µg/ml de rifampicine
)

# Repertoires de travail
#mkdir -p qc_reports nettoyage mapping reference comptage
input_dir="raw_data"
qc_dir="qc_reports"
clean_dir="nettoyage"
aligned_dir="mapping"
reference_dir="reference"
genome_ref="$reference_dir/h37rv.fna"
gtf_ref="$reference_dir/annotation.gtf"
compte_dir="comptage"

# Traitement 
# ================================================================
# 1) Contrôle qualité des données 
# ================================================================
#for srr_id in "${srr_list[@]}"; do
 #   fastqc "$input_dir/${srr_id}_1.fastq.gz" "$input_dir/${srr_id}_2.fastq.gz" -o "$qc_dir"
#done
  echo " Étape suivante"

# 1.2) Multiqc
#multiqc "$qc_dir"


# ================================================================
# 2) Nettoyage de reads
# ================================================================
#for srr_id in "${srr_list[@]}"; do
 #  fastp -i "$input_dir/${srr_id}_1.fastq.gz" \
  #       -I "$input_dir/${srr_id}_2.fastq.gz" \
   #      -o "$clean_dir/${srr_id}_1_clean.fastq.gz" \
    #     -O "$clean_dir/${srr_id}_2_clean.fastq.gz" 
#done
   
   # 2.2) Contrôle qualité des reads nettoyés
# Cette étape de contrôle qualité n'est pas forcement necessaire, vu que les reads sont déjà de bonne qualité. Mais si vous voulez confirmer, alors faites-le. Vous allez constatez aucun changement.

# mkdir -p qc_clean
# for srr_id in "${srr_list[@]}"; do
 #    fastqc "$clean_dir/${srr_id}_1_clean.fastq.gz" "$clean_dir/${srr_id}_2_clean.fastq.gz" -o qc_clean
#done
   
# ================================================================
# 3) Alignement
# ================================================================
#for srr_id in "${srr_list[@]}"; do
    # indexation
    #bowtie2-build "$genome_ref" "$reference_dir/h37rv_index"

    # mapping
    #bowtie2 -x "$reference_dir/h37rv_index" -1 "$clean_dir/${srr_id}_1_clean.fastq.gz" -2 "$clean_dir/${srr_id}_2_clean.fastq.gz" -S "$aligned_dir/${srr_id}.sam" --very-sensitive -p 8
            
            
   # samtools view -bS "$aligned_dir/${srr_id}.sam" | samtools sort -o "$aligned_dir/${srr_id}.bam"
    #samtools index "$aligned_dir/${srr_id}.bam"
    #rm -f "$aligned_dir/${srr_id}.sam"
            
#done

# ================================================================
# 4) Comptage des reads (featureCounts)
# ================================================================

#for srr_id in "${srr_list[@]}"; do
#   featureCounts -T 8 -p -t exon -g gene_id \
#   -a "$gtf_ref" \
#   -o "$compte_dir/counts.txt" mapping/*.sorted.bam
#done
            
for f in mapping/*.bam.bai.sorted.bam.bai; do
    mv "$f" "${f%.bam.bai.sorted.bam.bai}.sorted.bam.bai"
done
