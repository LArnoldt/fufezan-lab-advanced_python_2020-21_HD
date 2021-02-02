import argparse
import matplotlib.pyplot as plt
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

def plot_amino_acid_prevalences(amino_acids_dict):
    """Function plots a histogram for the prevalences of the amino acids in the dictionary.

    Args:
        amino_acids_dict: dictionary with prevalences for single amino acids
    """

    prevalences_amino_acids = list(amino_acids_dict.values())
    amino_acids = list(amino_acids_dict.keys())

    plt.bar(amino_acids, prevalences_amino_acids, color="green")
    plt.title("Prevalences of Amino Acids for Uniprot sequences of a species.")
    plt.ylabel("Prevalences")
    plt.xlabel("Amino Acids")
    plt.savefig("./exercises/e02/prevalence_distribution_amino_acids.pdf")
    plt.xlabel("Amino Acids")
    plt.savefig("./exercises/e02/prevalence_distribution_amino_acids.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", help="FASTA file with the sequences to be analyzed. ", type=str)
    args = parser.parse_args()
    fasta_file_string = args.fasta_file

    fasta_file = import_fasta_file(fasta_file_string)

    aggregated_fastas = aggregate_rows_in_fasta_file(fasta_file)

    amino_acids_dict = count_amino_acid_prevalences_in_aggregated_fastas(aggregated_fastas)

    plot_amino_acid_prevalences(amino_acids_dict)

    fasta_file.close()
