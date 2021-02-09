import pandas as pd
import argparse
import plotly.io as pio
pio.renderers.default = "browser"
import plotly.graph_objects as go

def import_table_with_pandas(directory_table):
    '''Function to import a table as a pandas dataframe.

    Args:
        directory_table: directory of a table
    Returns:
        dataframe: imported pandas dataframe
    '''

    dataframe = pd.read_csv(directory_table)

    return dataframe

def select_columns_in_pandas_table(dataframe, selected_columns):
    '''Function to select certain columns in a pandas dataframe and aggregate them in a new table, so all other columns are deleted.

    Args:
        dataframe: imported pandas dataframe
        selected_columns: selected columns to be kept in the dataframe
    Returns:
        dataframe: pandas dataframe with selected columns
    '''

    dataframe_selected_columns = pd.DataFrame([])
    for column in selected_columns:
        dataframe_selected_columns = dataframe_selected_columns.append(dataframe[column])
    dataframe_selected_columns = dataframe_selected_columns.T

    return dataframe_selected_columns

def rename_columns_in_pandas_table_with_spaces_instead_of_dots(dataframe):
    '''Function renames the columns if a pandas dataframe, so they do include a " " instead of ".".

    Args:
        dataframe: imported pandas dataframe
    Returns:
        dataframe: pandas dataframe with renamed columns
    '''

    new_dataframe_columns = []

    for element in dataframe.columns:
        new_dataframe_columns.append(element.replace(".", " "))

    dataframe.columns = new_dataframe_columns

    return dataframe

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

def identification_tasks(dataframe):
    '''Function fulfills several identification requirements defined in the tasks.

    Args:
        dataframe: pandas dataframe
    '''

    frequencies_countries = dataframe['Country.of.Origin'].value_counts().to_dict()

    countries_less_10_entries = []
    countries_less_30_entries = []

    for country, frequency in frequencies_countries.items():
        if frequency < 30:
            countries_less_30_entries.append(country)
            if frequency < 10:
                countries_less_10_entries.append(country)

    print("List of countries with less than 10 entries: " + str(countries_less_10_entries))
    print("List of countries with less than 30 entries: " + str(countries_less_30_entries))

    frequencies_producer = dataframe['Producer'].value_counts().to_dict()
    producer_most_entries = list(frequencies_producer.keys())[0]

    print("The producer with most entries ist: " + producer_most_entries)

    frequencies_processing_method = dataframe['Processing Method'].value_counts().to_dict()

    most_common_processing_method = list(frequencies_processing_method.keys())[0]
    least_common_processing_method = list(frequencies_processing_method.keys())[-1]

    print("Most common processing method: " + most_common_processing_method)
    print("Least common processing method: " + least_common_processing_method)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory_table", help="Directory of the imported table.", type=str, required=True)
    parser.add_argument("--selected_columns", help="Columns in table to be kept in pandas dataframe. Please type in as 'Column1,Column2,...'", type=str, required=True)
    args = parser.parse_args()

    directory_table = args.directory_table #directory_table = "./data/arabica_data_cleaned.csv"
    selected_columns = args.selected_columns #selected_columns = ["Country.of.Origin", "Producer", "Processing.Method"]
    selected_columns = selected_columns.split(",")

    dataframe = import_table_with_pandas(directory_table)
    dataframe = select_columns_in_pandas_table(dataframe, selected_columns)
    dataframe = rename_columns_in_pandas_table_with_spaces_instead_of_dots(dataframe)

    for column in dataframe.columns:
        plot_histogram_for_pandas_column(dataframe, column)

    identification_tasks(dataframe)
