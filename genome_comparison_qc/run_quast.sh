#!/bin/bash
# chmod +rwx filename

# quast version 5.2.0

cat ./paths_to_genomes.txt | while read file

do 
    echo "Starting QUAST on $file"
    quast.py -o "$(basename "$file")_quast_output" -m 0 -l "$(basename "$file")" $file
done