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
        longest_fasta: length of the longest protein sequence
        shortest_fasta: length of the shortest protein sequence
    """

    aggregated_fastas = ""
    longest_fasta = 0
    shortest_fasta = float('inf')
    counter = 0
    aggregated_fastas_protein = ""

    for row in fasta_file:
        if row[0] != ">":
            row = row[:-1]
            aggregated_fastas += row
            counter += 1
            aggregated_fastas_protein += row
        if row[0] == ">":
            counter = 0
            if len(aggregated_fastas_protein) > longest_fasta:
                longest_fasta = len(aggregated_fastas_protein)
            if len(aggregated_fastas_protein) < shortest_fasta and len(aggregated_fastas_protein) > 0:
                shortest_fasta = len(aggregated_fastas_protein)
            aggregated_fastas_protein = ""

    return aggregated_fastas, longest_fasta, shortest_fasta

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

    with open('./exercises/e02/mouse_aa_distribution.csv', 'w', newline='') as csvfile:
        aa_distribution_csv = csv.writer(csvfile, delimiter=' ')
        aa_distribution_csv.writerow(["aa, count"])
        for amino_acid in amino_acids_dict.keys():
            print(str(amino_acid) + ": " + str(amino_acids_dict[amino_acid]))
            aa_distribution_csv.writerow([str(amino_acid) + ", " + str(amino_acids_dict[amino_acid])])


def import_amino_acid_properties_and_create_pi_dict(directory_amino_acid_properties):
    '''Function imports a table with amino acid properties and extracts the hydropathy scores for certain bases.

    Args:
        directory_amino_acid_properties: directory of table with amino acid properties, e.g. hydropathy
    Returns:
        pi_dict: dict of pi scores for certain amino acids
    '''

    pi_dict = {}

    with open(directory_amino_acid_properties, newline='') as csvfile:
        amino_acid_properties = csv.reader(csvfile, delimiter=',', quotechar='|')
        for pos, row in enumerate(amino_acid_properties):
            if pos > 0:
                pi_dict[row[2]] = float(row[10])

    return pi_dict

def calculation_pi(amino_acids_dict, pi_dict):
    '''Function calculates the pI of a sequence.

    Args:
        amino_acids_dict: dictionary (Counter type) of amino acids with prevalences
        pi_dict: dict of pi scores for certain amino acids
    Returns:
        pi: pi score for a sequence
    '''

    pi = 0

    for amino_acid, count in amino_acids_dict.items():
        try:
            pi += pi_dict[amino_acid] * count
        except:
            print("The amino acid " + amino_acid + " was not given in the file amino_acid_properties.csv. Therefore it is ignored for the futher calculation.")

    pi = pi / sum(amino_acids_dict.values())

    return pi

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", help="FASTA file with the sequences to be analyzed. ", type=str)
    parser.add_argument("--directory_amino_acid_properties", help="Directory of CSV with amino acid properties",
                        type=str, required=False)
    args = parser.parse_args()
    fasta_file_string = args.fasta_file
    directory_amino_acid_properties = args.directory_amino_acid_properties

    fasta_file = import_fasta_file(fasta_file_string)

    aggregated_fastas, longest_fasta, shortest_fasta = aggregate_rows_in_fasta_file(fasta_file)

    amino_acids_dict = count_amino_acid_prevalences_in_aggregated_fastas(aggregated_fastas)

    output_amino_acid_prevalences_and_save_to_csv(amino_acids_dict)

    pi_dict = import_amino_acid_properties_and_create_pi_dict(directory_amino_acid_properties)
    pi = calculation_pi(amino_acids_dict, pi_dict)

    print("The total pI of the sequence is " + str(pi) + ".")
    print("The longest FASTA in the inputted FASTA file has a length of " + str(longest_fasta) + ".")
    print("The shortest FASTA in the inputted FASTA file has a length of " + str(shortest_fasta) + ".")

    fasta_file.close()