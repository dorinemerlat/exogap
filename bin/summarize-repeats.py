#!/usr/bin/env python3

import argparse
import pandas as pd
import plotly.express as px
import os
import kaleido

def split_te_name(df):
    # Split 'te_name' into 'class' and 'family' columns
    df[['class', 'family']] = df['te_name'].str.split('/', expand=True)

    # If there is no '/', set 'family' column to None
    df['family'].fillna(value=pd.NA, inplace=True)

    # Drop the original 'te_name' column
    df.drop(columns=['te_name'], inplace=True)

    return df

def combine_csv_files_in_directory():
    # Get the current directory
    current_directory = os.getcwd()

    # Get a list of all CSV files in the current directory
    csv_files = [f for f in os.listdir(current_directory) if f.endswith('.out.reformated')]

    # Create an empty DataFrame to store the combined data
    combined_data = pd.DataFrame()

    # Iterate through the list of CSV files, read each file, and concatenate its content
    for csv_file in csv_files:
        file_path = os.path.join(current_directory, csv_file)
        df = pd.read_csv(file_path)
        combined_data = pd.concat([combined_data, df], ignore_index=True)

    return combined_data

def create_genome_length_bar_plot(dataframe):
    fig = px.bar(dataframe, x='genome_name', y='genome_length',
                title='Genome Lengths',
                labels={'genome_name': 'Genome Name', 'genome_length': 'Genome Length'},
                orientation='h')
    fig.write_image("/home/merlat/sync/test.png")

# function to change genome_length from bp to Mb
def genome_length_to_Mb(df): 
    df['genome_length'] = df['genome_length'].apply(lambda x: x / 1000000)
    return df

def clean (df) :
    df = split_te_name(df)
    df = genome_length_to_Mb(df)

    return df

def overview_repeats(directory, classification, newick):
    df = combine_csv_files_in_directory()
    df = split_te_name(df)
    
df['genome_length'] = df['genome_length'].apply(lambda x: x / 1000000)

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Summarize repeats')
    parser.add_argument('-c', '--classification', help='Classification file')
    parser.add_argument('-n', '--newick', help='Newick file')
    args = parser.parse_args()

    # Collect all data in a single DataFrame
    df = combine_csv_files_in_directory()

###### DROP THIS SECTION ######
    # # save combined_df to csv
    # combined_df.to_csv('combined_df.csv', index=False)

    # # open combined_df.csv
    # df = pd.read_csv('combined_df.csv')
###############################

    # Reformat data
    df = split_te_name(df)

    # Genome lengths
    df['genome_length'] = df['genome_length'].apply(lambda x: x / 1000000)
    genome_lengths = df[['genome_id', 'genome_name', 'genome_length']].drop_duplicates()


    print(genome_lengths)
    create_genome_length_bar_plot(genome_lengths)
    print("finished")

