#!/bin/bash
# chmod +rwx filename

# busco -v ## 5.4.3 in conda environment

# busco -i [SEQUENCE_FILE] -l [LINEAGE/DATASET] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]
# vertebrata_odb10
# metazoa_odb10
# eukaryota_odb10

cat ./paths_to_genomes.txt | while read file

do 
    echo "Starting BUSCO on $file"
    busco -i $file -o "$(basename "$file")_busco_output" -m genome -l metazoa_odb10 --cpu 40
done