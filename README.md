# Mtb : ANALYSE RNA-SEQ
Analyse Transcriptomique de Mycobacterium tuberculosis sous Rifampicine

# OBJECTIF DU PROJET

Ce projet vise à analyser la réponse transcriptionnelle de Mycobacterium tuberculosis H37Rv à différentes concentrations de rifampicine au cours du temps. L'objectif est d'identifier les gènes différentiellement exprimés et de comprendre les mécanismes moléculaires de résistance et d'adaptation de la bactérie face à cet antibiotique.

# CONTEXTE BIOLOGIQUE

La rifampicine est un antibiotique de première ligne dans le traitement de la tuberculose. Comprendre la réponse transcriptionnelle de M. tuberculosis à cet antibiotique est crucial pour :

        -Identifier les mécanismes de tolérance et de résistance

        -Découvrir de nouvelles cibles thérapeutiques

        -Améliorer les stratégies de traitement

# DESIGN EXPERIMENTAL

Souche : M. tuberculosis H37Rv

Conditions :

        -Contrôle (sans antibiotique)

        -Rifampicine 0.008 µg/mL

        -Rifampicine 0.02 µg/mL


  +Points temporels : 4h, 24h, 72h post-traitement

  +Réplicats : 3 réplicats biologiques par condition

  +Total : 27 échantillons RNA-seq paired-end (Illumina NextSeq 500)

  +Source des données : GEO (PRJNA701434)



# STRUCTURE DU PROJET


                       SCRIPT

├── scripts/

│   ├── 01_download.sh                 # Téléchargement des données SRA

│   ├── 02_preprocessing.sh            # Pipeline complet (QC, nettoyage, alignement)

│   └── 03_differential_analysis.R     # Analyse différentielle avec DESeq2

├── data/
│   ├── SRR_Acc_List.txt              # Liste des identifiants SRA

│   └── metadata.csv                   # Métadonnées des échantillons

├── reference/
│   ├── h37rv.fna                      # Génome de référence (non inclus, à télécharger)

│   └── annotation.gtf                 # Annotation génique (non inclus, à télécharger)

├── results/
│   ├── DEGs/                          # Gènes différentiellement exprimés (CSV)

│   ├── figures/                       # Visualisations (PCA, Volcano, Heatmaps)

│   └── summary/                       # Résumés et statistiques

└── docs/
    └── workflow.md                    # Documentation détaillée du workflow


# PREREQUIS

    LOGICIELS REQUIS

        # Outils bioinformatiques
- SRA Toolkit (≥ 3.0)
- FastQC (≥ 0.11)
- MultiQC (≥ 1.12)
- fastp (≥ 0.23)
- Bowtie2 (≥ 2.4)
- SAMtools (≥ 1.15)
- Picard (≥ 2.27)
- Subread/featureCounts (≥ 2.0)

# Environnement R
- R (≥ 4.2)
- DESeq2
- ggplot2
- pheatmap
- dplyr
- RColorBrewer

La version des outils peuvent évoluer.


# WORKFLOWS D'ANALYSE

    Étape 1 : Téléchargement des Données
        # Télécharger les données brutes depuis SRA
            bash scripts/01_download.sh

    Étape 2 : Preprocessing
        # QC, nettoyage, alignement, comptage
            bash scripts/02_preprocessing.sh

    Étapes incluses :

        1. Contrôle qualité avec FastQC
        2. Nettoyage et déduplication avec fastp
        3. Alignement sur le génome de référence (Bowtie2)
        4. Conversion et tri des fichiers BAM (SAMtools)
        5. Suppression des duplications PCR (Picard)
        6. Comptage des reads par gène (featureCounts)

    Étape 3 : Analyse Différentielle
        # Analyse avec DESeq2
            Rscript scripts/03_differential_analysis.R

    Comparaisons effectuées :

        Rif 0.008 µg/mL vs Control (à 4h, 24h, 72h)
        Rif 0.02 µg/mL vs Control (à 4h, 24h, 72h)

# RESULTATS

    Les résultats incluent :

    1. Fichiers CSV : Liste des gènes différentiellement exprimés (6 comparaisons)

    2. PCA : Visualisation globale de la variabilité entre échantillons
    3. Volcano plots : Identification des gènes up/down-régulés (6 plots)

    4. Heatmaps : Expression des top DEGs par comparaison (8 heatmaps)

    5. Résumé statistique : Nombre de DEGs par condition


# Références des Données

    Étude originale : Mycobacterium tuberculosis transcriptional response to Rifampicin

    BioProject : PRJNA701434

    GEO Series : GSE166622  


# Notes

Les fichiers de référence (génome et annotation) doivent être téléchargés séparément depuis :

NCBI Reference Sequence: NC_000962.3 (M. tuberculosis H37Rv)
Annotation: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/


Les données brutes (fichiers FASTQ) ne sont pas incluses dans ce dépôt en raison de leur taille. Utilisez 01_download.sh pour les télécharger.

# Auteur
    Keny Karl MOUNGUELE 
Étudiant en Biotechnologie et Bioinformatique 

[Ton email/LinkedIn si tu veux]

# 🙏 Remerciements

Données RNA-seq : [Auteurs de l'étude originale]
Outils bioinformatiques : Communauté open-source