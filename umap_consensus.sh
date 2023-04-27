#!/bin/bash

echo "Please type in the path to the UMAP TSV file: "
read input_file
echo "Please type in the path to the FASTQ file (with the file): "
read input_fastq
mkdir "$input_file/split_OTU"

# Count the total number of bin_id
total_bin_id_count=$(awk 'NR>1{count[$5]++} END{for (bin_id in count) sum += count[bin_id]; print sum}' "$input_file"/umap_hdbscan_plot.tsv)

# Read the input file line by line
while IFS=$'\t' read -r read length D1 D2 bin_id; do
  # Create a separate .txt file for each bin_id in the output directory
  output_file="$input_file/split_OTU/$bin_id.txt"
  
  # Calculate the percentage of bin_id count
  bin_id_count=$(awk -v bin_id="$bin_id" 'NR>1{if ($5 == bin_id) count++} END{print count}' "$input_file"/umap_hdbscan_plot.tsv)
  percentage=$(echo "scale=2; ($bin_id_count * 100) / $total_bin_id_count" | bc)

  # Check if the percentage is more than or equal to 3.00
  if (( $(echo "$percentage >= 3.00" | bc -l) )); then
    # Extract read name from the line
    read_name=$(echo "$read" | sed 's/>//')
    
    # Search for the corresponding sequence in the input FASTQ file
    sequence=$(grep -A1 "$read_name" "$input_fastq" | tail -n1)
    
    # Append "@" to the front of the read name
    read_name=">${read_name}"
    
    # Append the read name and sequence to the respective bin_id .txt file in the output directory if the condition is true
    echo -e "$read_name\n$sequence" >> "$output_file.fasta"
  fi
done < "$input_file"/umap_hdbscan_plot.tsv

# Loop through each file in the directory
for file in "$input_file"/split_OTU/*; do
  # Get the file name without extension
  file_name=$(basename "$file" | cut -d. -f1)
  
  # Create a directory with the file name
  mkdir "$input_file"/split_OTU/"$file_name"

  # Move the file to the corresponding directory
  mv "$file" "$input_file"/split_OTU/"$file_name" 

  # Correct sequence with Canu
  canu -correct -p corrected_reads -d "$input_file"/split_OTU/"$file_name"/correct -nanopore-raw "$input_file"/split_OTU/"$file_name"/*.fasta genomeSize=1000 stopOnLowCoverage=0 minInputCoverage=0 minReadLength=500 minOverlapLength=200

  # Unzip corrected read
  gunzip -r "$input_file"/split_OTU/*

  # Create a directory with the file name for spliting read
  mkdir "$input_file"/split_OTU/"$file_name"/split_correct
  mkdir "$input_file"/split_OTU/"$file_name"/read_list

  # Split read
  split -l 2 "$input_file"/split_OTU/"$file_name"/correct/corrected_reads.correctedReads.fasta "$input_file"/split_OTU/"$file_name"/split_correct/split
  find "$input_file"/split_OTU/"$file_name"/split_correct/split* > "$input_file"/split_OTU/"$file_name"/read_list/read_list.txt

  #Draft selection
  mkdir "$input_file"/split_OTU/"$file_name"/draft
  fastANI --ql "$input_file"/split_OTU/"$file_name"/read_list/read_list.txt --rl "$input_file"/split_OTU/"$file_name"/read_list/read_list.txt -o "$input_file"/split_OTU/"$file_name"/read_list/fastani_output.txt -t 48 -k 16 --fragLen 160
  chmod +x "$input_file"/split_OTU/"$file_name"/read_list/fastani_output.txt
  awk 'NR>1{name[$1] = $1; arr[$1] += $3; count[$1] += 1}  END{for (a in arr) {print arr[a] / count[a], name[a] }}' "$input_file"/split_OTU/"$file_name"/read_list/fastani_output.txt > "$input_file"/split_OTU/"$file_name"/draft/selected.fasta
  sort -rg "$input_file"/split_OTU/"$file_name"/draft/selected.fasta | cut -d " " -f2 | head -n1 > "$input_file"/split_OTU/"$file_name"/draft/draft_dir.txt
  
  # Load the selected sequence to generate draft sequence
  file_dir="$input_file"/split_OTU/"$file_name"/draft/draft_dir.txt
    while IFS= read -r line; do
      # Process the current line (e.g., display or store it)
      cat "$line" > "$input_file"/split_OTU/"$file_name"/draft/draft.fasta
    done < "$file_dir"

  # Medaka_consensus pipeline
  medaka_consensus -i "$input_file"/split_OTU/"$file_name"/correct/corrected_reads.correctedReads.fasta -d "$input_file"/split_OTU/"$file_name"/draft/draft.fasta -o "$input_file"/split_OTU/"$file_name"/consensus -m r941_min_sup_g507
done 
