#!/usr/local/bin/python3
import panedr

# Read the EDR file
path = 'ener.edr'
df = panedr.edr_to_df(path)

# The `verbose` optional parameter can be set to True to display the
# progress on stdout
df = panedr.edr_to_df(path, verbose=True)

# Get the average pressure after the first 10 ns
pressure_avg = df[u'Pressure'][df[u'Time'] > 1000].mean()
print(df[u'Pressure'][df[u'Time']])
print(df.keys())
