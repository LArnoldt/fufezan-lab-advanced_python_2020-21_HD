import csv
import argparse
from collections import Counter

def import_fasta_file(fasta_file_string):
    """Imports the fasta file.

    Args:
        fasta_file_string: directory of the fasta file with the sequences to be analysed

    Returns:
        fasta_file: Loaded fasta file
    """

    fasta_file = open(fasta_file_string, "r")
    return fasta_file

def aggregate_rows_in_fasta_file(fasta_file):
    """Aggregates all rows in the fasta file.

    Args:
        fasta_file: Loaded fasta file

    Returns:
        aggregated_fastas: aggregated fasta sequences
    """

    aggregated_fastas = ""

    for row in fasta_file:
        if row[0] != ">":
            row = row[:-1]
            aggregated_fastas += row

    return aggregated_fastas

def count_amino_acid_prevalences_in_aggregated_fastas(aggregated_fastas):
    """Function outputs the calculated prevalences of the amino acids and saves them to the CSV

    Args:
        aggregated_fastas: aggregated fasta sequences

    Returns:
        amino_acids_dict: dictionary (Counter type) of amino acids with prevalences
    """

    amino_acids_dict = Counter(aggregated_fastas)
    return amino_acids_dict

def output_amino_acid_prevalences_and_save_to_csv(amino_acids_dict):
    """Function outputs the calculated prevalences of the amino acids and saves them to the CSV

    Args:
        amino_acids_dict: dictionary (Counter type) of amino acids with prevalences
    """

<<<<<<< HEAD
    with open('./exercises/e02/mouse_aa_distribution.csv', 'w', newline='') as csvfile:
        aa_distribution_csv = csv.writer(csvfile, delimiter=' ')
        aa_distribution_csv.writerow(["aa, count"])
        for amino_acid in amino_acids_dict.keys():
            print(str(amino_acid) + ": " + str(amino_acids_dict[amino_acid]))
            aa_distribution_csv.writerow([str(amino_acid) + ", " + str(amino_acids_dict[amino_acid])])
=======
    with open('./exercises/e02/aa_distribution.csv', 'w', newline='') as csvfile:
        human_aa_distribution_csv = csv.writer(csvfile, delimiter=' ', quoting=csv.QUOTE_MINIMAL)
        for amino_acid in amino_acids_dict.keys():
            print(str(amino_acid) + ": " + str(amino_acids_dict[amino_acid]))
            human_aa_distribution_csv.writerow([str(amino_acid) + ", " + str(amino_acids_dict[amino_acid])])
>>>>>>> Exercise 2 finalized.

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", help="FASTA file with the sequences to be analyzed. ", type=str)
    args = parser.parse_args()
    fasta_file_string = args.fasta_file

    fasta_file = import_fasta_file(fasta_file_string)

    aggregated_fastas = aggregate_rows_in_fasta_file(fasta_file)

    amino_acids_dict = count_amino_acid_prevalences_in_aggregated_fastas(aggregated_fastas)

    output_amino_acid_prevalences_and_save_to_csv(amino_acids_dict)

    fasta_file.close()