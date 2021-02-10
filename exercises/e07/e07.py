import pandas as pd
import argparse
import plotly.io as pio
pio.renderers.default = "browser"
import plotly.graph_objects as go

import ssl
import json
import urllib

def import_table_with_pandas_via_urllib(url):
    ssl._create_default_https_context = ssl._create_unverified_context

    json_unformated = urllib.request.urlopen(url).read().decode("utf-8")
    json = json.loads(json_unformated)
    df = pd.DataFrame(json['records'])

    return df

def rename_notification_rate_per_100000_population_14_days_to_14d_incidence(df):

    df.rename(
        columns={
            "notification_rate_per_100000_population_14-days": "14d-incidence",
        },
        inplace=True
    )

    return df

def create_datetime_objects(df):

    df['date_reported'] = pd.to_datetime(df['dateRep'])

    #TREATMENT OPTION
    #df['date_reported'].dt #.day

def create_column(df, column_name):

    column_name = "deltaTime_since_start_of_recording"

    df[column_name] = 10433 * [None]


def plot_histogram_for_pandas_column(dataframe, column):
    '''Function plots a histogram for a single column in a pandas dataframe.

    Args:
        dataframe: pandas dataframe
        column: column in dataframe with the data supposed to be plotted
    '''

    data = [
        go.Histogram(
            x=dataframe[column],
            )
    ]

    fig = go.Figure(data=data)
    fig.update_layout(template="plotly_dark", title="Histogram for the imported table for the column '" + column +"'")
    fig.layout.xaxis.title = 'Processing Methods'
    fig.layout.yaxis.title = 'Number of arabica bones, which we treated with this processing method'
    fig.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", help="URL of the imported table.", type=str, required=True) #url = "https://opendata.ecdc.europa.eu/covid19/casedistribution/json/"
    parser.add_argument("--selected_columns", help="Columns in table to be kept in pandas dataframe. Please type in as 'Column1,Column2,...'", type=str, required=True)
    args = parser.parse_args()

    directory_table = args.directory_table #directory_table = "./data/arabica_data_cleaned.csv"

    ####dataframe = import_table_with_pandas(directory_table)

    for column in df.columns:
        plot_histogram_for_pandas_column(df, column)
