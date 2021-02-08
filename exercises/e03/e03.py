import argparse
import requests
import plotly.express as px
import plotly.io as pio
pio.renderers.default = "browser"
import plotly.graph_objects as go
from pathlib import Path
import csv
from collections import deque
import pandas as pd
import xmlschema

def uniprot_download_single_entry(acc_uniprot):
   """Function downloads the corresponding .fasta file for a uniprot accession number and saves it.

   Args:
       acc_uniprot: uniprot accession number
   """

   url_entry_uniprot = "https://www.uniprot.org/uniprot/" + acc_uniprot + ".fasta"
   request_entry_uniprot = requests.get(url_entry_uniprot, allow_redirects=True)
   open("./exercises/e03/" + acc_uniprot + ".fasta", 'wb').write(request_entry_uniprot.content)

def extract_fasta(fasta_file_string):
    '''Function extracts FASTA string from fasta file.

    Args:
        fasta_file_string: Directory of FASTA file to be loaded.
    Returns:
        fasta: string with amino acids
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

    hydropathy_plot = px.bar(sequence_position_hydropathy_df, x = 'Sequence Position', y = 'Hydropathy', title="Hydropathy List for Uniprot Acc: " + acc_uniprot + " with a sliding window of " + str(length_sliding_window))
    hydropathy_plot.show()

def plot_hydropathy_list_and_transmembrane_region(hydropathy_sequence_list, acc_uniprot, length_sliding_window, domain_sites_list):
    ''' Function generates a pandas table with the sequence positions, the hydropathy_values and region sites and plots them in a bar chart.

    Args:
        hydropathy_sequence_list: hydropathy values for sequence positions
        acc_uniprot: uniprot accession number
        length_sliding_window: length of the sliding window used in the determination of the hydropathy
        domain_sites_list: list with topological domains and transmembrane region sites of the corresponding protein
    '''

    sequence_positions = []
    [sequence_positions.append(sequence_position) for sequence_position in range(0, len(hydropathy_sequence_list))]

    sequence_position_hydropathy_domains_tuple = list(zip(sequence_positions, hydropathy_sequence_list, domain_sites_list))
    sequence_position_hydropathy_domains_df = pd.DataFrame(sequence_position_hydropathy_domains_tuple, columns=['Sequence Position', 'Hydropathy', 'Domain'])

    data = [
        go.Bar(
            x=sequence_position_hydropathy_domains_df['Sequence Position'],
            y=sequence_position_hydropathy_domains_df['Hydropathy'],
            showlegend=False,
            marker_color=sequence_position_hydropathy_domains_df['Domain'].map(
                {
                    "transmembrane_region_sites": "blue",
                    "topological_domain_site": "green",
                }
            )
        )
    ]

    data.append(go.Scatter(x=[None], y=[None], mode='markers',
                           marker=dict(size=10, color='blue'),
                           legendgroup='Transmembrane region site', showlegend=True, name='Transmembrane region site'))
    data.append(go.Scatter(x=[None], y=[None], mode='markers',
                           marker=dict(size=10, color='green'),
                           legendgroup='Topological domain site', showlegend=True, name='Topological domain site'))

    fig = go.Figure(data=data)
    fig.update_layout(template="plotly_dark", title="Hydropathy List for Uniprot Acc: " + acc_uniprot + " with a sliding window of " + str(length_sliding_window))
    fig.update_layout(showlegend=True)
    fig.layout.xaxis.title = "Sequence Positions"
    fig.layout.yaxis.title = "Hydropathy of the protein residues (Kyte-Doolittle Method)"
    fig.layout.legend.title = "Region Sites"
    fig.show()

def extract_xml_and_extract_topological_domains(directory_xml_file):
    '''Function imports the XML file of a protein and extracts all topological domains.

    Args:
        directory_xml_file: directory of the XML file for a corresponsing protein
    Returns:
        domain_sites_list: list with topological domains and transmembrane region sites of the corresponding protein
    '''

    schema = xmlschema.XMLSchema('https://www.uniprot.org/docs/uniprot.xsd')
    entry_dict = schema.to_dict(directory_xml_file)
    entry_dict.keys()
    content = entry_dict['entry'][0]

    domain_sites_list = []
    end_old = 0

    for element in content['feature']:
        if element['@type'] == 'topological domain':
            begin, end = element['location']['begin']['@position'], element['location']['end']['@position']
            for amino_acid_residue_position in range(end_old, begin-1):
                domain_sites_list.append("transmembrane_region_sites")
            for amino_acid_residue_position in range(begin-1, end):
                domain_sites_list.append("topological_domain_site")
                end_old = end

    return domain_sites_list

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta_file", help="FASTA file with the sequences to be analyzed.", type=str, default="NaN", required=False)
    parser.add_argument("--acc_uniprot", help="Uniprot accession number.", type=str, default="NaN", required=False)
    parser.add_argument("--length_sliding_window", help="Length of the sliding window.", type=int, required=True)
    parser.add_argument("--directory_amino_acid_properties", help="Directory of CSV with amino acid properties", type=str, required=True)
    parser.add_argument("--directory_xml_file", help="Directory of CSV with amino acid properties",type=str, required=True)
    args = parser.parse_args()

    fasta_file_string = args.fasta_file
    acc_uniprot = args.acc_uniprot
    length_sliding_window = args.length_sliding_window
    directory_amino_acid_properties = args.directory_amino_acid_properties
    directory_xml_file = args.directory_xml_file

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
        domain_sites_list = extract_xml_and_extract_topological_domains(directory_xml_file)
        plot_hydropathy_list_and_transmembrane_region(hydropathy_sequence_list, acc_uniprot.stem, length_sliding_window, domain_sites_list)
        #Without the domain sites:
        #plot_hydropathy_list(hydropathy_sequence_list, acc_uniprot.stem, length_sliding_window)