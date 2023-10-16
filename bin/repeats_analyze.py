# Split 'te_name' into 'class' and 'family' columns (family = None if there is only class indicated)
def split_te_name(df):
    df[['class', 'family']] = df['te_name'].str.split('/', expand=True)
    df['family'].fillna(value=pd.NA, inplace=True)
    df.drop(columns=['te_name'], inplace=True)

    return df

# function to change genome_length from bp to Mb
def genome_length_to_Mb(df): 
    df['genome_length'] = df['genome_length'].apply(lambda x: x / 1000000)
    return df

# 
def get_repeats_df(df):
    df = split_te_name(df)
    df = genome_length_to_Mb(df)
    return df
