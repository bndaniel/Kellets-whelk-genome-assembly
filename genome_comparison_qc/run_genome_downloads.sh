#!/bin/bash
# chmod +rwx filename

# conda install -c conda-forge ncbi-datasets-cli

cat ../genome_accessions.txt | while read accession

do 
    echo "Current download: $accession"
    datasets download genome accession $accession --include genome --filename "${accession}_dataset.zip"
done