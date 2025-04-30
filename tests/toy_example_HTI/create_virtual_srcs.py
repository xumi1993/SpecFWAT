import numpy as np
import pandas as pd


def to_forcesolution(lat, lon):
    solution = f'''
FORCE  001
time shift:     0.0000
f0:             0.5
latorUTM:       {lat}
longorUTM:      {lon}
depth:          0.0000
source time function:            0
factor force source:             1.d15
component dir vect source E:     0.d0
component dir vect source N:     0.d0
component dir vect source Z_UP:  1.d0
'''
    return solution

def create_virtual_srcs():
    stas = pd.read_csv('DATA/STATIONS', header=None, sep=r'\s+')
    sources = []
    for i in range(stas.shape[0]-1):
        sources.append([
            stas.iloc[i, 0],  # station name
            stas.iloc[i, 2],  # latitude
            stas.iloc[i, 3],  # longitude
            stas.iloc[i, 4],  # elevation
            stas.iloc[i, 5],  # depth
            1.0
        ])
        stations = []
        for j in range(i + 1, stas.shape[0]):
            stations.append([
                stas.iloc[j, 0],
                stas.iloc[j, 1],  # station name
                stas.iloc[j, 2],  # latitude
                stas.iloc[j, 3],  # longitude
                stas.iloc[j, 4],  # elevation
                stas.iloc[j, 5],  # depth
            ])
        stations = pd.DataFrame(stations)
        stations.to_csv(
            f'src_rec/STATIONS_{stas.iloc[i, 0]}',
            header=False,
            index=False,
            sep=' ',
        )
        # create the force solution file
        with open(f'src_rec/FORCESOLUTION_{stas.iloc[i, 0]}', 'w') as f:
            f.write(to_forcesolution(stas.iloc[i, 2], stas.iloc[i, 3]))
    sources = pd.DataFrame(sources)
    sources.to_csv(
        'src_rec/sources_noise.dat',
        header=False,
        index=False,
        sep=' ',
    )

if __name__ == '__main__':
    create_virtual_srcs()
