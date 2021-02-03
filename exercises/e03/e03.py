import argparse
import requests
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "browser"
from pathlib import Path
import csv
from collections import deque
import pandas as pd

def uniprot_download_single_entry(acc_uniprot):
   """Function downloads the corresponding .fasta file for a uniprot accession number and saves it.

   Args:
       acc_uniprot: uniprot accession number
   """

   url_entry_uniprot = "https://www.uniprot.org/uniprot/" + acc_uniprot + ".fasta"
   request_entry_uniprot = requests.get(url_entry_uniprot, allow_redirects=True)
   open("./exercises/e03/" + acc_uniprot + ".fasta", 'wb').write(request_entry_uniprot.content)

def extract_fasta(fasta_file_string):
    '''Function extracs FASTA string from fasta file.

    Args:
        fasta_file_string: Directory of FASTA file to be loaded.
    Returns:
        fasta: string with bases
    '''

    fasta = ""

    with open(fasta_file_string, "r") as fasta_file:
        fasta_file_lines = fasta_file.readlines()
        for row in fasta_file_lines:
            if row[0] != ">":
                row = row[:-1]
                fasta += row

    return fasta

def import_amino_acid_properties_and_create_hydropathy_dict(directory_amino_acid_properties):
    '''Function imports a table with amino acid properties and extracts the hydropathy scores for certain bases.

    Args:
        directory_amino_acid_properties: directory of table with amino acid properties, e.g. hydropathy
    Returns:
        hydropathy_dict: dict of hydropathy scores for certain bases
    '''
    
    hydropathy_dict = {}

    with open(directory_amino_acid_properties, newline='') as csvfile:
        amino_acid_properties = csv.reader(csvfile, delimiter=',', quotechar='|')
        for pos, row in enumerate(amino_acid_properties):
            if pos > 0:
                hydropathy_dict[row[2]] = float(row[11])

    return hydropathy_dict

def create_hydropathy_value_list_for_fasta(fasta, hydropathy_dict, length_sliding_window):
    '''Function generates a list of hydropathy scores for the sequence making use of a sliding window.

    Args:
        fasta: fasta sequence to be converted into hydropathy scores
        hydropathy_dict: dict of hydropathy scores for certain bases
        length_sliding_window: length of the sliding window used in the determination of the hydropathy
    Returns:
        hydropathy_sequence_list: hydropathy values for sequence positions
    '''

    hydropathy_sequence_list = []

    '''
    #Version for Task B)
    for base in fasta:
        hydropathy_sequence_list.append(hydropathy_dict[base])

    return hydropathy_sequence_list
    '''

    #Version for Task D)

    window = deque([], maxlen=length_sliding_window)

    for pos, aa in enumerate(fasta):
        window.append(aa)
        average_hydropathy_score = 0
        for base in list(window):
            average_hydropathy_score += hydropathy_dict[base]
        average_hydropathy_score = average_hydropathy_score / len(window)
        hydropathy_sequence_list.append(average_hydropathy_score)

    return hydropathy_sequence_list

def plot_hydropathy_list(hydropathy_sequence_list, acc_uniprot, length_sliding_window):
    ''' Function generates a pandas table with the sequence positions and the hydropathy_values and plots them in a bar chart.

    Args:
        hydropathy_sequence_list: hydropathy values for sequence positions
        acc_uniprot: uniprot accession number
        length_sliding_window: length of the sliding window used in the determination of the hydropathy
    '''

    sequence_positions = []
    [sequence_positions.append(sequence_position) for sequence_position in range(0, len(hydropathy_sequence_list))]

    sequence_position_hydropathy_tuple = list(zip(sequence_positions, hydropathy_sequence_list))
    sequence_position_hydropathy_df = pd.DataFrame(sequence_position_hydropathy_tuple, columns=['Sequence Position', 'Hydropathy'])

    hydropathy_plot = px.bar(sequence_position_hydropathy_df, x = 'Sequence Position', y = 'Hydropathy', title="Hydropathy List for Uniprot Acc: " + acc_uniprot + " with a sliding window of " + length_sliding_window)
    hydropathy_plot.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", help="FASTA file with the sequences to be analyzed.", type=str, default="NaN", required=False)
    parser.add_argument("--acc_uniprot", help="Uniprot accession number.", type=str, default="NaN", required=False)
    parser.add_argument("--length_sliding_window", help="Length of the sliding window.", type=int, required=True)
    parser.add_argument("--directory_amino_acid_properties", help="Directory of CSV with amino acid properties", type=str, required=True)
    args = parser.parse_args()

    fasta_file_string = args.fasta_file
    acc_uniprot = args.acc_uniprot
    length_sliding_window = args.length_sliding_window
    directory_amino_acid_properties = args.directory_amino_acid_properties

    breaking = False

    if acc_uniprot != "NaN":
        uniprot_download_single_entry(acc_uniprot)
        fasta = extract_fasta("./exercises/e03/" + acc_uniprot + ".fasta")
        acc_uniprot = Path(acc_uniprot)
    elif fasta_file_string != "NaN":
        fasta = extract_fasta(fasta_file_string)
        acc_uniprot = Path(fasta_file_string)
    else:
        print("No FASTA file, neither a Uniprot accession number was added. Break!")
        breaking = True

    if breaking is not True:

        hydropathy_dict = import_amino_acid_properties_and_create_hydropathy_dict(directory_amino_acid_properties)
        hydropathy_sequence_list = create_hydropathy_value_list_for_fasta(fasta, hydropathy_dict, length_sliding_window)
        plot_hydropathy_list(hydropathy_sequence_list, acc_uniprot.stem, length_sliding_window)