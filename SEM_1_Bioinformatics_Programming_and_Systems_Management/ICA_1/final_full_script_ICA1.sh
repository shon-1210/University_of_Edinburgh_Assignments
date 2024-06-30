#!/bin/bash

###Entire script for the ICA

###Initialising variables that will be used

###these are variables for model 3###
export FQ_FILES=()
export filename=""
export REFERENCE=""
export GENE_FILES=()

echo "please run this script from the homespace if there are any errors"

mkdir folder_to_run_ICA

cd folder_to_run_ICA

###Code for task 1

cp -r /localdisk/data/BPSM/ICA1/ .       ###dont forget the . at the end ### -r is to copy all files and subdirectories recursively
cd ICA1/
mkdir fastqc_after_unzip
cd fastq
gunzip *.gz                              #### unzipping the .gz files
mv *.fq ../fastqc_after_unzip           ##moving to new folder to do fastqc
cd ../fastqc_after_unzip
fastqc *.fq
unzip '*.zip'                           ###unzipping the zip files


### Code for task 2


echo " The quality of the raw sequence data can be seen looking at the summary file for each sample, example path is (homespace)/folder_to_run_ICA/ICA1/fastqc_after_unzip/Tco-106_1_fastqc/summary.txt"




### Code for task 3
cd ~/folder_to_run_ICA/ICA1
mkdir bowtie_task
cd Tcongo_genome/
gunzip TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta.gz
cd ../bowtie_task/
cp -r ../Tcongo_genome/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta .        ##### dont forget the . at the end
cp -r ../fastqc_after_unzip/*.fq .                               ####dont forget the . at the end, we need the .fq files to do bowtie




cd ../fastqc_after_unzip/

# Create an empty array
FQ_FILES=()

for file in *.fq;
do
file_without_extension="${file%%_*}"
FQ_FILES+=("$file_without_extension")
done


#####The %% operator in ${file%%_*} removes the longest matching suffix pattern. It tries to match '_*' (underscore followed by any string) from the back of the string and removes it.

export  FQ_FILES

#if you want to view the array files and confirm the array works


filename="fq_files_for_alignment.txt"
for element in "${FQ_FILES[@]}"
do
        echo "$element" >> "$filename"
done

#cat fq_files_for_alignment.txt

#### code to do the actual alignment####
### make sure it is a local alignment ###



### make the index ##

cd ../bowtie_task/
bowtie2-build TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta  bowtie_database_index #this step is to create the index



# Path to the reference genome index prefix
REFERENCE="${HOME}/folder_to_run_ICA/ICA1/bowtie_task/bowtie_database_index" #IS THE SAME AS BELOW


# Array of your gene sequence files (without extensions)
GENE_FILES=("${FQ_FILES[@]}") # and so on

# Loop through the array of genes
for GENE in "${GENE_FILES[@]}"
do
  # Make a new directory for the current gene
  mkdir -p "${GENE}_output"

  # Navigate into the directory (safe practice for file generation)
  pushd "${GENE}_output"

  # As we have paired-end reads, align with something like:
  bowtie2 -x "${REFERENCE}" -1 "../${GENE}_1.fq" -2 "../${GENE}_2.fq" -S "${GENE}.sam"

  # Convert SAM to BAM, sort and index it
  samtools view -S -b "${GENE}.sam" | samtools sort -o "${GENE}.bam"
  samtools index "${GENE}.bam"

  # Remove the SAM file if you no longer need it
  rm "${GENE}.sam"

  # Navigate back to the parent directory
  popd
done






######This is the code for task 4########



####THIS IS THE CODE FOR MODEL 4 WHICH USES BEDTOOLS


cd ~/folder_to_run_ICA/ICA1
mkdir bedtools_multiple_full_1

cd bedtools_multiple_full_1
mkdir counts_data

cp ../*.bed .      #bed file from ICA1 folder
cp -r ../bowtie_task/*_output .    ###Dont forget the . at the end


for file in *output
        {
                cd ${file}
                cp -r *.bam ../counts_data
                cd ..
        }

cd counts_data
for file in *bam
        {
                bedtools coverage -counts -a ../TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${file} > ${file}.txt

        }




###we need to group according to the groups as seen in more Tco2.fqfiles


#### now we have to group the samples before proceeding###


####code for task 5 full thing###
###arrive to the correct folder###

cd ~/folder_to_run_ICA/ICA1/

mkdir sorting_groups
cd sorting_groups
cp ~/folder_to_run_ICA/ICA1/fastq/Tco2.fqfiles .     ###dont forget the . at the end





# Arrays for the different variables

clones=("WT" "Clone1" "Clone2")
times=("0" "24" "48")
states=("Induced" "Uninduced")



### this is to make the final group table ###
### needed for referencing if we got the right result ###

# Loop through each combination of variables
for clone in "${clones[@]}"; do
    for time in "${times[@]}"; do
        for state in "${states[@]}"; do
           # Extracting data from Tco2.fqfiles based on the combination of variables
            result=$(cat Tco2.fqfiles | cut -f 2,4,5 | grep "$clone" | grep "$time" | grep "$state")
            echo -e  "$clone\t$time\t$state" >> sorted_output.txt
        done
    done
done



awk -F'\t' ' !($2 == "0" && $3 == "Induced") ' sorted_output.txt > refined.txt    ###this is remove groups with 0h and Induced parameter

cat Tco2.fqfiles | cut -f 1,2,4,5 > all_the_fqfiles.txt      ##putting all the names of Tco2.fqfiles into a txt file for easy use



###code to make groups


mkdir groups_of_samples

# Define the arrays with the different parameters
clones=("WT" "Clone1" "Clone2")
times=("0" "24" "48")
states=("Induced" "Uninduced")

# Looping through all combinations of clone, time, and state
for clone in "${clones[@]}"; do
    for time in "${times[@]}"; do
        for state in "${states[@]}"; do
            # Constructing the filename for the current group
            group_filename="${clone}_${time}h_${state}_group.txt"

            # Use awk to filter the original file based on the current parameters
            awk -v clone="$clone" -v time="$time" -v state="$state" \
            '$2 == clone && $3 == time && $4 == state' "all_the_fqfiles.txt" > "$group_filename"
        done
    done
done

rm *0h_Induced_group.txt   ### to remove the txt files with 0h and induced as parameters

mv *group.txt groups_of_samples/

cd groups_of_samples

wc -l *group.txt > number_of_files_per_group.txt  ###to see number of samples that are in each group




##### This is to calculate mean and generate files with the mean value as the last column ####
cd ~/folder_to_run_ICA/ICA1/sorting_groups/groups_of_samples/

# Path to the directory containing the Tco-*.txt files
DATA_DIR="${HOME}/folder_to_run_ICA/ICA1/bedtools_multiple_full_1/counts_data"

# Create a new file for output
OUTPUT_FILE="all_samples_with_means_output.txt"
echo -e "SampleName\tSampleType\tTime\tTreatment\tMean" > "$OUTPUT_FILE"

# Reading all_the_fqfiles.txt line by line
while IFS=$'\t' read -r SampleName SampleType Time Treatment; do
  if [ "$SampleName" != "SampleName" ]; then  # Skip the header

    # Remove 'Tco' from the sample name and construct the file name
    stripped_sample_name="${SampleName#Tco}"
    FILENAME="$DATA_DIR/Tco-${stripped_sample_name}.bam.txt"
    #Checking if the file exists
    if [[ -f "$FILENAME" ]]; then
      # Calculate the mean using awk. This assumes that the values are in the last column.
      mean=$(awk '{ total += $NF; count++ } END { print total/count }' "$FILENAME")

      # Append the data to the output file
      echo -e "$SampleName\t$SampleType\t$Time\t$Treatment\t$mean" >> "$OUTPUT_FILE"
    else
      echo "File not found: $FILENAME"
    fi
  fi
done < "${HOME}/folder_to_run_ICA/ICA1/sorting_groups/all_the_fqfiles.txt"

# The results will be in all_samples_with_means_output.txt #



### now we need to calculate the average for the groups ##



# Name of the input and output file
INPUTFILE="${HOME}/folder_to_run_ICA/ICA1/sorting_groups/groups_of_samples/all_samples_with_means_output.txt"
OUTPUTFILE="${HOME}/folder_to_run_ICA/ICA1/sorting_groups/groups_of_samples/grouped_means.txt"

# First we will print the header to the output file. This is so that we obtain the output with header line first
echo -e "SampleType\tTime\tTreatment\tMean" > $OUTPUTFILE

# Using awk to process the data file, we calculate the mean for each group, and add the data below our header line in the output file.
awk 'BEGIN {
    FS=OFS="\t"  # Setting field separator to tab
}
NR > 1 {  # Here we are Ignoring the header line (first line) of the input file
    group = $2 FS $3 FS $4  # Grouping by SampleType, Time, Treatment
    total[group] += $5  # Summing the Mean values
    count[group]++  # Counting the number of entries per group
}
END {
    for (group in total) {
        # Calculating the mean and printing the group along with its calculated mean
        print group, total[group] / count[group]
    }
}' $INPUTFILE | sort -k1,1 -k2,2n -k3,3 -s >> $OUTPUTFILE  # Sorting data and appending to output file, below the header


### in the above awk script in the sort line, n is for sorting numerically and s is for maintaining the order of equal key records in the output as was in the input ###

echo "The means of the groups are stored in: $OUTPUTFILE"

echo "path of the file is ~/folder_to_run_ICA/ICA1/sorting_groups/groups_of_samples/grouped_means.txt"

echo "The gene descriptions were not added"



