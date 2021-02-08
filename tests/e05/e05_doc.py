import pytest
import sys, os

from exercises.e05.e05_doc import Protein

def test_get_sequence_1():

    protein = Protein("P32249", 5, "./data/amino_acid_properties.csv", "./exercises/e03/P32249.xml")
    fasta = protein.get_sequence()

    assert fasta == "MDIQMANNFTPPSATPQGNDCDLYAHHSTARIVMPLHYSLVFIIGLVGNLLALVVIVQNRKKINSTTLYSTNLVISDILFTTALPTRIAYYAMGFDWRIGDALCRITALVFYINTYAGVNFMTCLSIDRFIAVVHPLRYNKIKRIEHAKGVCIFVWILVFAQTLPLLINPMSKQEAERITCMEYPNFEETKSLPWILLGACFIGYVLPLIIILICYSQICCKLFRTAKQNPLTEKSGVNKKALNTIILIIVVFVLCFTPYHVAIIQHMIKKLRFSNFLECSQRHSFQISLHFTVCLMNFNCCMDPFIYFFACKGYKRKVMRMLKRQVSVSISSAVKSAPEENSREMTETQMMIHSKSSNGK"

def test_get_sequence_2():

    protein = Protein("", 5, "./data/amino_acid_properties.csv", "./exercises/e03/P32249.xml")
    fasta = protein.get_sequence()

    assert fasta == ""

def test_get_sequence_3():

    protein = Protein("DS", 5, "./data/amino_acid_properties.csv", "./exercises/e03/P32249.xml")
    fasta = protein.get_sequence()

    assert fasta == ""

def test_import_amino_acid_properties_and_create_hydropathy_dict_1():

    protein = Protein("P32249", 5, "./data/amino_acid_properties.csv", "./exercises/e03/P32249.xml")
    hydropathy_dict = protein.import_amino_acid_properties_and_create_hydropathy_dict()

    assert hydropathy_dict == {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}


def test_import_amino_acid_properties_and_create_hydropathy_dict_2():

    protein = Protein("P32249", 5, "", "./exercises/e03/P32249.xml")
    hydropathy_dict = protein.import_amino_acid_properties_and_create_hydropathy_dict()

    assert hydropathy_dict == {}

def test_create_hydropathy_value_list_for_fasta_1():

    protein = Protein("P32249", 5, "./data/amino_acid_properties.csv", "./exercises/e03/P32249.xml")
    hydropathy_dict = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2,
                        'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
                        'Y': -1.3, 'V': 4.2}
    hydropathy_sequence_list = protein.create_hydropathy_value_list_for_fasta("MDIQMANNFTPPSATPQGNDCDLYAHHSTARIVMPLHYSLVFIIGLVGNLLALVVIVQNRKKINSTTLYSTNLVISDILFTTALPTRIAYYAMGFDWRIGDALCRITALVFYINTYAGVNFMTCLSIDRFIAVVHPLRYNKIKRIEHAKGVCIFVWILVFAQTLPLLINPMSKQEAERITCMEYPNFEETKSLPWILLGACFIGYVLPLIIILICYSQICCKLFRTAKQNPLTEKSGVNKKALNTIILIIVVFVLCFTPYHVAIIQHMIKKLRFSNFLECSQRHSFQISLHFTVCLMNFNCCMDPFIYFFACKGYKRKVMRMLKRQVSVSISSAVKSAPEENSREMTETQMMIHSKSSNGK", hydropathy_dict)
    hydropathy_sequence_list_solution = [1.9,
 -3.5,
 4.5,
 -3.5,
 1.9,
 1.8,
 -3.5,
 -3.5,
 2.8,
 -0.7,
 -1.6,
 -1.6,
 -0.8,
 1.8,
 -0.7,
 -1.6,
 -3.5,
 -0.4,
 -3.5,
 -3.5,
 2.5,
 -3.5,
 3.8,
 -1.3,
 1.8,
 -3.2,
 -3.2,
 -0.8,
 -0.7,
 1.8,
 -4.5,
 4.5,
 4.2,
 1.9,
 -1.6,
 3.8,
 -3.2,
 -1.3,
 -0.8,
 3.8,
 4.2,
 2.8,
 4.5,
 4.5,
 -0.4,
 3.8,
 4.2,
 -0.4,
 -3.5,
 3.8,
 3.8,
 1.8,
 3.8,
 4.2,
 4.2,
 4.5,
 4.2,
 -3.5,
 -3.5,
 -4.5,
 -3.9,
 -3.9,
 4.5,
 -3.5,
 -0.8,
 -0.7,
 -0.7,
 3.8,
 -1.3,
 -0.8,
 -0.7,
 -3.5,
 3.8,
 4.2,
 4.5,
 -0.8,
 -3.5,
 4.5,
 3.8,
 2.8,
 -0.7,
 -0.7,
 1.8,
 3.8,
 -1.6,
 -0.7,
 -4.5,
 4.5,
 1.8,
 -1.3,
 -1.3,
 1.8,
 1.9,
 -0.4,
 2.8,
 -3.5,
 -0.9,
 -4.5,
 4.5,
 -0.4,
 -3.5,
 1.8,
 3.8,
 2.5,
 -4.5,
 4.5,
 -0.7,
 1.8,
 3.8,
 4.2,
 2.8,
 -1.3,
 4.5,
 -3.5,
 -0.7,
 -1.3,
 1.8,
 -0.4,
 4.2,
 -3.5,
 2.8,
 1.9,
 -0.7,
 2.5,
 3.8,
 -0.8,
 4.5,
 -3.5,
 -4.5,
 2.8,
 4.5,
 1.8,
 4.2,
 4.2,
 -3.2,
 -1.6,
 3.8,
 -4.5,
 -1.3,
 -3.5,
 -3.9,
 4.5,
 -3.9,
 -4.5,
 4.5,
 -3.5,
 -3.2,
 1.8,
 -3.9,
 -0.4,
 4.2,
 2.5,
 4.5,
 2.8,
 4.2,
 -0.9,
 4.5,
 3.8,
 4.2,
 2.8,
 1.8,
 -3.5,
 -0.7,
 3.8,
 -1.6,
 3.8,
 3.8,
 4.5,
 -3.5,
 -1.6,
 1.9,
 -0.8,
 -3.9,
 -3.5,
 -3.5,
 1.8,
 -3.5,
 -4.5,
 4.5,
 -0.7,
 2.5,
 1.9,
 -3.5,
 -1.3,
 -1.6,
 -3.5,
 2.8,
 -3.5,
 -3.5,
 -0.7,
 -3.9,
 -0.8,
 3.8,
 -1.6,
 -0.9,
 4.5,
 3.8,
 3.8,
 -0.4,
 1.8,
 2.5,
 2.8,
 4.5,
 -0.4,
 -1.3,
 4.2,
 3.8,
 -1.6,
 3.8,
 4.5,
 4.5,
 4.5,
 3.8,
 4.5,
 2.5,
 -1.3,
 -0.8,
 -3.5,
 4.5,
 2.5,
 2.5,
 -3.9,
 3.8,
 2.8,
 -4.5,
 -0.7,
 1.8,
 -3.9,
 -3.5,
 -3.5,
 -1.6,
 3.8,
 -0.7,
 -3.5,
 -3.9,
 -0.8,
 -0.4,
 4.2,
 -3.5,
 -3.9,
 -3.9,
 1.8,
 3.8,
 -3.5,
 -0.7,
 4.5,
 4.5,
 3.8,
 4.5,
 4.5,
 4.2,
 4.2,
 2.8,
 4.2,
 3.8,
 2.5,
 2.8,
 -0.7,
 -1.6,
 -1.3,
 -3.2,
 4.2,
 1.8,
 4.5,
 4.5,
 -3.5,
 -3.2,
 1.9,
 4.5,
 -3.9,
 -3.9,
 3.8,
 -4.5,
 2.8,
 -0.8,
 -3.5,
 2.8,
 3.8,
 -3.5,
 2.5,
 -0.8,
 -3.5,
 -4.5,
 -3.2,
 -0.8,
 2.8,
 -3.5,
 4.5,
 -0.8,
 3.8,
 -3.2,
 2.8,
 -0.7,
 4.2,
 2.5,
 3.8,
 1.9,
 -3.5,
 2.8,
 -3.5,
 2.5,
 2.5,
 1.9,
 -3.5,
 -1.6,
 2.8,
 4.5,
 -1.3,
 2.8,
 2.8,
 1.8,
 2.5,
 -3.9,
 -0.4,
 -1.3,
 -3.9,
 -4.5,
 -3.9,
 4.2,
 1.9,
 -4.5,
 1.9,
 3.8,
 -3.9,
 -4.5,
 -3.5,
 4.2,
 -0.8,
 4.2,
 -0.8,
 4.5,
 -0.8,
 -0.8,
 1.8,
 4.2,
 -3.9,
 -0.8,
 1.8,
 -1.6,
 -3.5,
 -3.5,
 -3.5,
 -0.8,
 -4.5,
 -3.5,
 1.9,
 -0.7,
 -3.5,
 -0.7,
 -3.5,
 1.9,
 1.9,
 4.5,
 -3.2,
 -0.8,
 -3.9,
 -0.8,
 -0.8,
 -3.5,
 -0.4,
 -3.9]

    assert hydropathy_sequence_list == hydropathy_sequence_list_solution

def test_create_hydropathy_value_list_for_fasta_2():

    protein = Protein("P32249", 5, "./data/amino_acid_properties.csv", "./exercises/e03/P32249.xml")
    hydropathy_dict = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'E': -3.5, 'Q': -3.5, 'G': -0.4, 'H': -3.2,
                        'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9,
                        'Y': -1.3, 'V': 4.2}
    hydropathy_sequence_list = protein.create_hydropathy_value_list_for_fasta("", hydropathy_dict)

    assert hydropathy_sequence_list == []

def test_create_hydropathy_value_list_for_fasta_3():

    protein = Protein("P32249", 5, "./data/amino_acid_properties.csv", "./exercises/e03/P32249.xml")
    hydropathy_sequence_list = protein.create_hydropathy_value_list_for_fasta("MDIQMANNFTPPSATPQGNDCDLYAHHSTARIVMPLHYSLVFIIGLVGNLLALVVIVQNRKKINSTTLYSTNLVISDILFTTALPTRIAYYAMGFDWRIGDALCRITALVFYINTYAGVNFMTCLSIDRFIAVVHPLRYNKIKRIEHAKGVCIFVWILVFAQTLPLLINPMSKQEAERITCMEYPNFEETKSLPWILLGACFIGYVLPLIIILICYSQICCKLFRTAKQNPLTEKSGVNKKALNTIILIIVVFVLCFTPYHVAIIQHMIKKLRFSNFLECSQRHSFQISLHFTVCLMNFNCCMDPFIYFFACKGYKRKVMRMLKRQVSVSISSAVKSAPEENSREMTETQMMIHSKSSNGK", {})

    assert hydropathy_sequence_list == []

def test_extract_xml_and_extract_topological_domains_1():

    protein = Protein("P32249", 5, "./data/amino_acid_properties.csv", "./exercises/e03/P32249.xml")
    domain_sites_list = protein.extract_xml_and_extract_topological_domains()

    assert domain_sites_list == []

def test_extract_xml_and_extract_topological_domains_2():

    protein = Protein("P32249", 5, "./data/amino_acid_properties.csv", "")
    domain_sites_list = protein.extract_xml_and_extract_topological_domains()
    domain_sites_list_solution = ['topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'transmembrane_region_sites',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site',
 'topological_domain_site']

    assert domain_sites_list == domain_sites_list_solution

