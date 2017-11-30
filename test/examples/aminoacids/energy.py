#!/usr/local/bin/python3
import panedr
import numpy

# Read the EDR file
path = 'ener.edr'
df = panedr.edr_to_df(path)

# The `verbose` optional parameter can be set to True to display the
# progress on stdout
#df = panedr.edr_to_df(path, verbose=True)
df = panedr.edr_to_df(path, verbose=False)

# Get the average pressure after the first 10 ns
pressure_avg = df[u'Temperature'][df[u'Time']].mean()
print(pressure_avg)
nframes=len(df[u'Temperature'][df[u'Time']])
print(nframes)
print(df.keys())
#for i in df[u'T-Argon'][df[u'Time']]:
#    print(i)
if 'Time' in df.keys():
    print("Yes")

for dfkey in df.keys():
    if 'Time' in dfkey:
        time = numpy.asarray(df[dfkey])
    if 'LJ' in dfkey:
        lj = numpy.asarray(df[dfkey])
    if 'Potential' in dfkey:
        pot = numpy.asarray(df[dfkey])
print(time)
print(lj)
print(pot)
vir = numpy.zeros([nframes,3,3])
for step in range(nframes):
    vir[step] = numpy.asarray([
        df[u'Vir-XX'][time[step]], df[u'Vir-XY'][time[step]], df[u'Vir-XZ'][time[step]],
        df[u'Vir-YX'][time[step]], df[u'Vir-YY'][time[step]], df[u'Vir-YZ'][time[step]],
        df[u'Vir-ZX'][time[step]], df[u'Vir-ZY'][time[step]], df[u'Vir-ZZ'][time[step]]
        ]).reshape((3, 3))
print(vir)
print(vir.shape)
print(df[u'Time'])





