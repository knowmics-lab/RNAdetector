############################### METADATA #######################################
LABEL software="RNAdetector" \
      version="1.0" \
      about.summary="Comprehensive RNA-Seq pipeline for the analysis of genes and ncRNAs" \
      name.maintainer="Alessandro La Ferlita" \
      email.maintainer="alessandrolf90@hotmail.it"

######################## Base images and scripts #######################
# Install Ubuntu
FROM ubuntu:18.04
RUN sudo apt-get update && sudo apt-get install -y \
    curl \
    tar
# Create bash scripts directory
RUN mkdir -p /usr/app/bash_scripts
# Copy bash scripts
COPY /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/trim_galore.bash /usr/app/bash_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/samtools.bash /usr/app/bash_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/tophat.bash /usr/app/bash_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/bwa.bash /usr/app/bash_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/salmon_index_2.sh /usr/app/bash_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/salmon_counting.sh /usr/app/bash_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/ciri.bash /usr/app/bash_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/featurecounts.bash /usr/app/bash_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/bash_scripts/htseqcount.bash /usr/app/bash_scripts

# Create R scripts directory
RUN mkdir -p /usr/app/r_scripts
# Copy R Scripts
COPY /home/alaferlita/Scrivania/RNAdetector/Scripts/R_scripts/deseq.R /usr/app/r_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/R_scripts/edgeR.R /usr/app/r_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/R_scripts/limma.R /usr/app/r_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/R_scripts/normalization.R /usr/app/r_scripts \
     /home/alaferlita/Scrivania/RNAdetector/Scripts/R_scripts/format_count_out.R /usr/app/r_scripts

FROM r-base
FROM python
FROM asdc/apache2-php7
FROM biocontainers/samtools
FROM chrishah/trim_galore
FROM genomicpariscentre/tophat2
FROM biocontainers/bwa
FROM cursecatcher/ciri2
FROM genomicpariscentre/htseq
FROM alexiswl/featurecounts_1.5.3
FROM combinelab/salmon
FROM genomicpariscentre/deseq2
FROM bioneos/edger
FROM biocontainers/bioconductor-limma
# MITHrIL
FROM openjdk
COPY /home/alaferlita/Scrivania/MITHrIL_2/MITHrIL2.jar /usr/app/scripts/
