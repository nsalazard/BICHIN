import pandas as pd

# Load the CSV data into a pandas DataFrame
df = pd.read_csv('config.csv')
TMAX=df['TMAX'][0]
# Print the DataFrame to see its contents
print(df)
print(TMAX)