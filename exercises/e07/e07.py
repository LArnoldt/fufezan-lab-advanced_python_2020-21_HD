import pandas as pd
import argparse
import plotly.io as pio
pio.renderers.default = "browser"
import plotly.graph_objects as go
import ssl
import json
import urllib
import numpy as np
from datetime import datetime, date

def import_table_with_pandas_via_urllib(url):
    ssl._create_default_https_context = ssl._create_unverified_context

    json_unformated = urllib.request.urlopen(url).read().decode("utf-8")
    json_file = json.loads(json_unformated)
    df = pd.DataFrame(json_file['records'])

    return df

def rename_notification_rate_per_100000_population_14_days_to_14d_incidence_and_preprocess_data(df):

    df.rename(
        columns={
            "notification_rate_per_100000_population_14-days": "14d-incidence",
        },
        inplace=True
    )

    df['14d-incidence'] = pd.to_numeric(df['14d-incidence'])

    return df

def create_datetime_objects(df):

    df['date_reported'] = pd.to_datetime(df['dateRep'], format='%d/%m/%Y', errors='raise')
    df['year'] = df['date_reported'].dt.year.values

    return df

def create_deltaTime(df):

    column_name = "deltaTime_since_start_of_recording"

    today = date.today()
    today = pd.to_datetime(today)

    df[column_name] = df.shape[0] * [None]
    for line in range(0,df.shape[0]):
        df.loc[line, column_name] = df.loc[line, 'date_reported'] - today

    df['cases_weekly'].mask(df['cases_weekly'] < 0, np.nan, inplace=True)
    df['deaths_weekly'].mask(df['deaths_weekly'] < 0, np.nan, inplace=True)
    df['14d-incidence'].mask(df['14d-incidence'] < 0, np.nan, inplace=True)
    df = df[df['cases_weekly'].notna()]
    df = df[df['deaths_weekly'].notna()]
    df = df[df['14d-incidence'].notna()]

    return df

def plot_histogram_for_pandas_column(df, column):
    '''Function plots a histogram for a single column in a pandas dataframe.

    Args:
        dataframe: pandas dataframe
        column: column in dataframe with the data supposed to be plotted
    '''

    data = [
        go.Histogram(
            x=df[column],
            )
    ]

    fig = go.Figure(data=data)
    fig.update_layout(template="plotly_dark", title="Histogram for the imported table for the column '" + column +"'")
    fig.show()

def top_values(df, col='diff_14d-incidence'):
    return df.sort_values(col, ascending=False).head(2)

def bottom_values(df, col='diff_14d-incidence'):
    return df.sort_values(col, ascending=True).head(2)

def visualization_countries_grouped_by_countries_for_14dincidence(df):

    grp = df.groupby(["countriesAndTerritories"])
    for country in grp["countriesAndTerritories"].unique():
        grp_country = grp.get_group(country[0]).sort_values('deltaTime_since_start_of_recording')
        for pos, incidence in enumerate(grp_country['14d-incidence']):
            if pos == 0:
                df.loc[grp_country.index[pos], 'diff_14d-incidence'] = float('NaN')
            else:
                df.loc[grp_country.index[pos], 'diff_14d-incidence'] = (grp_country.iloc[pos]['14d-incidence'] - grp_country.iloc[pos - 1]['14d-incidence'])

    continent_countries_top = df.groupby(['continentExp']).apply(top_values)

    data = [
        go.Bar(
            x = continent_countries_top["countriesAndTerritories"],
            y = continent_countries_top["diff_14d-incidence"],
        )
    ]

    fig = go.Figure(data=data)
    fig.update_layout(template="plotly_dark", title="Countries with the steepest increase of the 14d-incidence.", xaxis=dict(title="Countries"), yaxis=dict(title="Change of 14 days incidence"))
    fig.show()

    continent_countries_bottom = df.groupby(['continentExp']).apply(bottom_values)

    data = [
        go.Bar(
            x = continent_countries_bottom["countriesAndTerritories"],
            y = continent_countries_bottom["diff_14d-incidence"],
        )
    ]

    fig = go.Figure(data=data)
    fig.update_layout(template="plotly_dark", title="Countries with the steepest decrease of the 14d-incidence.", xaxis=dict(title="Countries"), yaxis=dict(title="Change of 14 days incidence"))
    fig.show()

    largest_diff_14d_incidence = df.nlargest(1, ['diff_14d-incidence'])
    country_largest_diff_14d_incidence = largest_diff_14d_incidence['countriesAndTerritories']
    print("The country with the largest fluctuation of the 14d-incidence is " + country_largest_diff_14d_incidence.unique()[0] + " with an incidence of " + str(float(largest_diff_14d_incidence['diff_14d-incidence'])))

    df['diff_14d-incidence'] = abs(df['diff_14d-incidence'])

    smallest_diff_14d_incidence = df.nsmallest(1, ['diff_14d-incidence'])
    country_smallest_diff_14d_incidence = smallest_diff_14d_incidence['countriesAndTerritories']
    print("The country with the smallest fluctuation of the 14d-incidence is " + country_smallest_diff_14d_incidence.unique()[0] + " with an incidence of " + str(float(smallest_diff_14d_incidence['diff_14d-incidence'])))

def line_plot_14dincidence(df):

    df_europe = df.groupby(["continentExp"]).get_group('Europe')
    df_europe = df_europe.groupby(["countriesAndTerritories"])

    fig = go.Figure()

    for country_territory in df_europe["countriesAndTerritories"].unique():
        country_territory_europe = df_europe.get_group(country_territory[0]).sort_values('deltaTime_since_start_of_recording')
        fig.add_trace(go.Scatter(x=country_territory_europe['date_reported'],
                             y=country_territory_europe['14d-incidence'],
                             mode='lines+markers',
                             name=country_territory[0]))

    fig.update_layout(template="plotly_dark", title="14 days incidence over time for European countries", xaxis=dict(title="Time (Months.)"), yaxis=dict(title="14 days incidence"))
    fig.show()

def smoothed_line_plot_14d_incidence(df):

    smoothed_df_europe = df.groupby(["continentExp"]).get_group('Europe').copy()
    smoothed_df_europe_14d_incidence = smoothed_df_europe['14d-incidence'].rolling(window=13, center = True).mean()
    smoothed_df_europe['smoothed_14d-incidence'] = smoothed_df_europe_14d_incidence
    smoothed_df_europe = smoothed_df_europe.groupby(["countriesAndTerritories"])

    fig = go.Figure()

    for country_territory in smoothed_df_europe["countriesAndTerritories"].unique():
        smoothed_country_territory_europe = smoothed_df_europe.get_group(country_territory[0]).sort_values(
            'deltaTime_since_start_of_recording')
        fig.add_trace(go.Scatter(x=smoothed_country_territory_europe['date_reported'],
                                 y=smoothed_country_territory_europe['smoothed_14d-incidence'],
                                 mode='lines+markers',
                                 name=country_territory[0]))

    fig.update_layout(template="plotly_dark", title="Smoothed 14 days incidence over time for European countries",
                      xaxis=dict(title="Time (Months.)"), yaxis=dict(title="14 days incidence"))
    fig.show()


def radial_axis_plot(df):
    mask_countries = (df['countriesAndTerritories'] == 'Germany') | (df['countriesAndTerritories'] == 'Italy') | (df['countriesAndTerritories'] == 'Sweden') | (df['countriesAndTerritories'] == 'Greece')
    df_countries = df[mask_countries].copy()

    df_countries = df_countries[['date_reported', 'popData2019', 'countriesAndTerritories', 'deaths_weekly']].groupby(['countriesAndTerritories'])

    data_radial_plot = []

    for country, df in df_countries:
        day_in_year = pd.to_datetime(df['date_reported']) - pd.to_datetime(2020, format='%Y')

        data_radial_plot.append(
            go.Scatterpolar(r=(df['deaths_weekly'] * 100000) / df['popData2019'],
                            theta=day_in_year.dt.days * 360 / 365,
                            name=country,))

    layout = {
        'template': 'plotly_dark',
        'title': {'text': 'Death rate for Covid-19 for the countries Germany, Italy, Sweden and Greece'},
        'polar': {
            'angularaxis': {
                'tickmode': 'array',
                'tickvals': [0, 72, 144, 216, 288],
                'ticktext': ['Day 0', 'Day 73', 'Day 146', 'Day 219', 'Day 292']
            },
            'radialaxis': {'dtick': 2, }
        }
    }

    fig = go.Figure(data=data_radial_plot, layout=layout)
    fig.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--url", help="URL of the imported table.", type=str, required=True) #url = "https://opendata.ecdc.europa.eu/covid19/casedistribution/json/"
    args = parser.parse_args()

    url = args.url

    df = import_table_with_pandas_via_urllib(url)

    df = rename_notification_rate_per_100000_population_14_days_to_14d_incidence_and_preprocess_data(df)
    df = create_datetime_objects(df)
    df = create_deltaTime(df)

    for column in df.columns:
        plot_histogram_for_pandas_column(df, column)

    visualization_countries_grouped_by_countries_for_14dincidence(df)
    line_plot_14dincidence(df)
    smoothed_line_plot_14d_incidence(df)
    radial_axis_plot(df)