#Importing the necessary libraries:

import argparse
import os
from Bio import SeqIO
from Bio.Seq import Seq
from matplotlib import pyplot as plt
from collections import defaultdict, Counter
from pyfaidx import Fasta
import logging

#Setting up Logger:

logging.basicConfig(filename = "error_log.log", level = logging.INFO,
                    format = '%(asctime)s - %(levelname)s - %(message)s')

#Setting up the Command Line Parser:

def parse_args():
    
    """This Function Helps to Interact with the Script Through the Command Line."""
    
    parser = argparse.ArgumentParser(description = "Processing the VCF, GFF and Genome Fasta Files" )
    parser.add_argument("--vcf", required = True, help = "Provide VCF Filepath")
    parser.add_argument("--gff", required = True, help = "Provide GFF Filepath")
    parser.add_argument("--fasta", required = True, help = "Provide Genome Fasta Filepath")
    return parser.parse_args()

#Function to Handle the GFF file:

def GFF_file_handler(path_to_gff):
    
    """Handles the GFF file and extracts the coding sequence (CDS) regions."""
    
    gff_content = []
    with open(path_to_gff) as file:
        for line in file:
            if not line.startswith("#"):
                sections = line.strip().split('\t')
                
                #Extracting Coding Sequence Regions Only
                
                if sections[2] == "CDS":
                    
                    #Storing only information of relevance:
                    
                    gff_content.append({
                        "sequence_id": sections[0],
                        "start": int(sections[3]),
                        "end": int(sections[4]),
                        "strand": sections[6],
                        "attributes": sections[8]
                    })
    return gff_content

#Function to Handle the VCF file:

def VCF_file_handler(path_to_vcf, gff_content, fasta):
    
    """This function handles a VCF file to categorize the variants into: Non-coding, Synonymous and Non-Synonymous Variants."""
    
    #Initiating Empty List to Store Results and Tracker to Count Variants of Poor Confidence:
    
    results = []
    poor_quality_count = 0 
    
    #Prepping the File to Write Output Too:
    
    with open(path_to_vcf) as vcf, open("Output_Variants.tsv", 'w') as output_tsv:
        output_tsv.write("Chrom\tPos\tRef\tAlt\tType\tTranscript\tProtein Location\tRef AA\tAlt AA\n")

        for line in vcf:
            if line.startswith("#"):  #Searching for header row
                continue
            
            sections = line.strip().split("\t")
            
            # Handling Missing Fields in VCF data:
            
            if len(sections) < 6:
                logging.error(f" VCF file supposed to have atleast 6 columns, but only has {len(sections)}): {line.strip()}")
                continue 
            
            # Ensuring the VCF Fields are Loaded Correctly:
            
            chrom = sections[0]
            position = int(sections[1])
            original = sections[3] if len(sections) > 3 else "N/A"
            alternative = sections[4] if len(sections) > 4 else "N/A"
            quality = float(sections[5]) if len(sections) > 5 else 0.0  # If missing, value is accomodated

            # Ensuring Variants are of Good Quality:
            
            if quality <= 20:
                poor_quality_count += 1
                continue

            # Creating Dictionary for Variants:
            
            variant = {"Chromosome": chrom, "Position": position, "Original": original, "Alternate": alternative}
            coding_region = coding_region_checker(variant, gff_content)

            if coding_region:
                original_amino, ugly_amino = amino_seq_checker(variant, coding_region, fasta)
                
                #Ensuring Sequences are Read in the Correct Direction using Strand for Reference:
                
                if coding_region["strand"] == "+":
                    protein_location = (variant["Position"] - coding_region["start"]) // 3 + 1
                elif coding_region["strand"] == "-":
                    protein_location = (coding_region["end"] - variant["Position"]) // 3 + 1
                else:
                    protein_location = "N/A"
                
                #Extracting the Required Info from Annotations in GFF File:
                
                transcript_id = coding_region["attributes"].split(";")[0].split("=")[-1]
                
                #Categorizing the Variants into Respective Sets:
                
                variant_type = "Non-Synonymous" if original_amino != ugly_amino else "Synonymous"
            else:
                transcript_id, protein_location, original_amino, ugly_amino = "NA", "NA", "NA", "NA"
                variant_type = "Non-coding"

            # Store the result and write to output file
            
            results.append(variant_type)
            output_tsv.write(f"{chrom}\t{position}\t{original}\t{alternative}\t{variant_type}\t{transcript_id}\t{protein_location}\t{original_amino}\t{ugly_amino}\n")

    return results, poor_quality_count

#Function to Verify Coding Region:

def coding_region_checker(variant, gff_content):
    
    """This function verifies whether a variant exists within a coding region."""
    
    for section in gff_content:
        
        #Checking if Variant is within a CDS region using boolean logic:
        
        if section["sequence_id"] == variant["Chromosome"] and section["start"] <= variant["Position"] <= section["end"]:
            return section
    return None

#Function to Determine the Change in Amino Acids due to Variants:

def amino_seq_checker(variant, coding_region, fasta):
    
    """This Function verifies and notes the amino acid change caused by a genetic variant."""
    
    try:
        reference_sequence = fasta[coding_region["sequence_id"]][coding_region["start"]-1:coding_region["end"]].seq 
        translated_sequence = Seq(reference_sequence).translate()

        #Calculating the Position of the Mutation in Protein Sequence:
        
        protein_position = (variant["Position"] - coding_region["start"]) // 3
        
        #Ensuring that the Position of the Protein is Always Within Range:
        
        protein_position = min(protein_position, len(translated_sequence) - 1) 
        
        #Updating Logger if Protein Position is not Within Valid Range:
        
        if not (0 <= protein_position < len(translated_sequence)):
            logging.error(f"Invalid position: {protein_position}")
            return "N/A", "N/A"

        #Obtaining the Original Amino Acid:
        
        original_amino = translated_sequence[protein_position]

        #Casting the Sequence from String to a List:
        
        mutated_sequence = list(reference_sequence)
        
        #Obtaining the position of the Mutation:
        
        mutation_index = variant["Position"] - coding_region["start"]

        #Checking if the Mutation is within Bounds and Logs an Error Otherwise:
        
        if mutation_index < 0 or mutation_index >= len(mutated_sequence):
            logging.error(f"Mutation index {mutation_index} out of range of sequence length {len(mutated_sequence)}")
            return "N/A", "N/A"
        
        #Accounting for Deletions, if not then we Replace with the Alternate:
        
        if variant["Alternate"] != "-":
            mutated_sequence[mutation_index] = variant["Alternate"]
        
        #Converting the Sequence with the Mutation into a String and Translating:
        
        mutated_sequence = "".join(mutated_sequence)
        mutated_translation = Seq(mutated_sequence).translate()

        #Checking if Protein Position is within Valid Range otherwise Logs an Error:
        
        if protein_position < 0 or protein_position >= len(mutated_translation):
            logging.error(f"Protein position {protein_position} out of range in mutated sequence length {len(mutated_translation)}.")
            return "N/A", "N/A"
        
        #Obtaining the Amino Acid from the Mutated Position:
        
        ugly_amino = mutated_translation[protein_position]

        return original_amino, ugly_amino
    
    #Stores Null for any Error and updates Logs:
    
    except Exception as e:
        logging.error(f"Error in computing amino acid")
        return "N/A", "N/A"
    
#Configuring Final Output Format and Logger Logs:

def final_output_generator(results, poor_quality_count):
    
    """This Function ensures that all the output, i.e., Summary Statistics, log results and a plot to visualise variant distribution, is formatted and shown to us clearly."""
    
    total = Counter(results)
    
    #Prepping Trackers to Log Number of Variants in Each Category:
    
    non_coding_count = total.get("Non-coding", 0)
    synonymous_count = total.get("Synonymous", 0)
    nonsynonymous_count = total.get("Non-Synonymous", 0)
    
    #Visualising the Classification Distribution:
    
    plt.bar(total.keys(), total.values())
    plt.xlabel("Variant_Type")
    plt.ylabel("Amount")
    plt.title("Proportion_of_Variant_Types")
    plt.savefig("Variant_Numbers.png")
    
    #Updating Log Messages to Show all Desirable Info:
    
    logging.info(f"Low Quality Variants (Quality <= 20): {poor_quality_count}")
    logging.info(f"Variants that have been categorized: {len(results)} Variants")
    logging.info(f" - Non-coding variants: {non_coding_count}")
    logging.info(f" - Synonymous variants: {synonymous_count}")
    logging.info(f" - Non-Synonymous variants: {nonsynonymous_count}")
    logging.info("Results saved: View Output_Variants.tsv and Variant_Numbers.png")

#Main function:

if __name__ == "__main__":
    args = parse_args()
    
    #Accounting for Missing Files:

    missing_files = [file for file in [args.vcf, args.gff, args.fasta] if not os.path.exists(file)]
    if missing_files:
        logging.error(f"Missing files: \n{missing_files}")
        exit(1)

    try:
        
        #Putting it all Together and calling Functions in Order:
        
        gff_content = GFF_file_handler(args.gff)
        fasta = Fasta(args.fasta)
        
        results, poor_quality_count = VCF_file_handler(args.vcf, gff_content, fasta)
        final_output_generator(results, poor_quality_count)
        
    except Exception as e:
        logging.error(f"An error occured: {e}")
        