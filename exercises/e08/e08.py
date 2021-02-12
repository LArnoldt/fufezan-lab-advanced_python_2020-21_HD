import argparse
import subprocess
from os import path
import pandas as pd
import numpy as np

def download_msp_file(url_msp_file, string_msp_file):
    print("The download and extraction of the MSP file started. Downloading may take some time.")
    subprocess.call([f"bash download_extract_msp_file.sh {url_msp_file} {string_msp_file}"], shell=True)
    print("The download and extraction of the MSP file finished.")

def msp_to_df(
        string_msp_file,
        max_seq_len=30,
        min_ce=36,
        max_ce=40,
        mz_min=135,
        mz_max=1400,
    ):
        """
        Function to read spectrum data from .msp file and convert to dataframe.
        Args:
            string_msp_file (str): path to .msp file
            max_seq_len (int): maximum acceptable sequence length
            min_ce (int): minimum collision energy of spectra to be included in df
            max_ce (int): maximum collision energy of spectra to be included in df
            mz_min (int): lower boundary for m/z to be included in df
            mz_max (int): lower boundary for m/z to be included in df

        Returns:
            df (pd.DataFrame or np.array): spectrum information within defined parameters [n_spectra, n_features]
            seqs (pd.DataFrame or np.array): sequences
        """

        df_mz = pd.DataFrame(np.zeros, index=range(0, 0), columns=range(mz_min, mz_max + 1))
        df_sequences = pd.DataFrame(np.nan, index=range(0, 0), columns=["Sequences"])

        count = 0
        go = False
        not_needed = 0
        mz_old = 0

        msp_file = open(string_msp_file, 'r')
        for line in msp_file:
            if not_needed != 0:
                not_needed -= 1
                continue
            elif line == "\n":
                go = False
            else:
                if line[0:4:] == "Name":
                    line = line[:-1]
                    line = line.split(": ")[1]
                    sequence, energy = line.split("/")
                    energy = float(energy.split("_")[-1][:-2])
                    if len(sequence) < max_seq_len and max_ce > energy > min_ce:
                        df_sequences.loc[count] = sequence
                        count += 1
                        not_needed = 3
                        go = True
                        mz_old = 0
                elif go == True:
                    mz = int(round(float(line.split("\t")[0]),0))
                    intensity = float(line.split("\t")[1])
                    if mz >= mz_min and mz <=mz_max:
                        if mz_old == mz:
                            intensity = max(intensity,  df_mz.loc[count, mz])
                        df_mz.loc[count, mz] = intensity
                    mz_old = mz

        df_mz.replace(np.nan, 0)

        return df_mz, df_sequences

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--url_msp_file", help="URL of the MSP file.", type=str, required=True)
    args = parser.parse_args()

    url_msp_file = args.url_msp_file  # url_msp_file = "ftp://chemdata.nist.gov/download/peptide_library/libraries/cptaclib/2015/cptac2_mouse_hcd_selected.msp.tar.gz"

    firstpos = url_msp_file.rfind("/")
    lastpos = len(url_msp_file)
    string_msp_file = "./exercises/e08/" + url_msp_file[firstpos + 1:lastpos]

    if not path.exists(string_msp_file):
        download_msp_file(url_msp_file, string_msp_file)

    string_msp_file = string_msp_file[:-7]