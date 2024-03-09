#!/usr/bin/python3

import glob
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.transforms import blended_transform_factory
from obspy import read, Stream
from obspy.geodetics import gps2dist_azimuth

st = Stream()
for sacf in glob.glob('processedSeismograms/Event_2014_02_26_21_13_40/*HHZ.sac'):
    st += read(sacf)

for tr in st:
    tr.stats.distance = tr.stats.sac.dist
    tr.data=np.require(tr.data, dtype=np.float32)  ### solve the error of data type 

#st.filter('bandpass', freqmin=0.02, freqmax=0.5)

fig = plt.figure()
st.plot(type='section', plot_dx=20e3,recordlength=300,
        linewidth=1.5, grid_linewidth=0.25, show=False, fig=fig)
plt.show()
