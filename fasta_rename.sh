#!/bin/bash

#SCRIPT TO MARGE FASTA SEQEUNCES FROM IRMA FOLDER AND THEN RENAME GENE SEQUENCES.

folder="$1" #Folder CONTAINING IRMA ASSEMBLIES

# Check if the provided path is a directory
if [ ! -d "$folder" ]; then
    echo "Error: $folder is not a valid directory"
    exit 1
fi

# Get the full path of the folder
full_path=$(readlink -f "$folder")
OUTPUT_DIR=$full_path/merged_sequences
mkdir -p $OUTPUT_DIR

for directory in $full_path/*; do
 	name=$(basename $directory)
	if [ $name != 'merged_sequences' ]; then #LOOP THROUGH FOLDERS 	EXCEPT merged_sequences
 		echo Processsing $name ...................
 		cd $directory
	 	cat *.fasta > $name'_merged.fasta'
	 	sed "s/>/>${name}_/" $name'_merged.fasta' >  $OUTPUT_DIR/$name'_merged.fasta'

	 	echo deleting some files
	 	rm $name'_merged.fasta'
	 	echo Complete Processsing $name ...................
	 	cd ../
	fi

done
