import pandas as pd
import argparse
import plotly.io as pio
pio.renderers.default = "browser"
import plotly.graph_objects as go




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory_table", help="Directory of the imported table.", type=str, required=True)
    parser.add_argument("--selected_columns", help="Columns in table to be kept in pandas dataframe. Please type in as 'Column1,Column2,...'", type=str, required=True)
    args = parser.parse_args()

    directory_table = args.directory_table #directory_table = "./data/arabica_data_cleaned.csv"

    dataframe = import_table_with_pandas(directory_table)