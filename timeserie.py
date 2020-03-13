'''
Purpose: Create a time serie from a folder with HMI fits files.

Author: Willinton Caicedo
        @ Universidad Nacional de Colombia - OAN - Mar 2020

Inputs:     * Folder Path
            * Flare peak time

Output: PDF file containing a time series plot.

History: Mar-2020: Created - Willy.
         Mar-2020: Modified - Milo. Increased efficiency.
'''

# User inputs:
tflare = "2013-11-08T04:26"  # time of flare peak (keep format)
path = 'HOME/20131108T0426/' # folder path


# Packages
from sunpy.map import Map
from astropy.table import Table
from scipy import ndimage
import sunpy.timeseries as ts
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import astropy.units as u
import datetime
import glob

# Read fits files and make maps
file_list = []
file_list.append(glob.glob(path+'*.i*.fits'))
file_list=sorted(file_list)
Mfiles=Map(file_list)

for mf in Mfiles: # changing 'CRDER1' and 'CRDER2' keywords to avoid annoying warnings
    mf.meta['CRDER1'] = 0.0
    mf.meta['CRDER2'] = 0.0

# Máscara:
imask = 10
diff = Mfiles[imask+1].data-Mfiles[imask].data
Mdiff = Map(np.nan_to_num(np.abs(diff)),Mfiles[imask].meta)
Mdiffrot = Mdiff.rotate(angle=Mdiff.meta['crota2'] * u.deg)
mask = Mdiffrot.data > 10
Gdiff = ndimage.gaussian_filter(Mdiffrot.data,40)
Gmask = Gdiff > Gdiff.max()*0.5
labels, n = ndimage.label(Gdiff*Gmask)
flaremask = (labels==3)

# Creating arrays containing time-serie data
Int_Inc = []
tiempos = []
for i in range(20):
    diff = Mfiles[i+1].data-Mfiles[i].data
    Mdiff = Map(np.nan_to_num(np.abs(diff)),Mfiles[i].meta)
    Mdiffrot = Mdiff.rotate(angle=Mdiff.meta['crota2'] * u.deg)
    Int_Inc.append((Mdiffrot.data*flaremask).sum())
    tiempos.append(datetime.datetime.strptime(Mdiffrot.meta["date-obs"],'%Y-%m-%dT%H:%M:%S.%f'))

tbl_meta = {'t_key':'t_value'}
table = Table([tiempos, Int_Inc/np.max(Int_Inc)], names=['time', 'Inclination'], meta=Mfiles[i].meta)
table.add_index('time')
ts_table = ts.TimeSeries(table)

# PLOT
fig, ax = plt.subplots(figsize=(10,4))
ts_table.plot(marker='o',linewidth=3)
ax.axvline(tflare, color="gray", linestyle="--")
ax.text(tflare, np.min(Int_Inc/np.max(Int_Inc)), 'flare peak', fontsize=14,color='gray',rotation=90, rotation_mode='anchor')
ax.tick_params(axis='both',labelsize=14)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
ax.set_xlabel('Time [hour:min]',fontsize=16)
ax.set_ylabel('Normalized B-Inclination diff ',fontsize=14)
ax.set_title(Mfiles[imask].meta['date-obs'],fontsize=24)
plt.legend('',frameon=False) # fix bug with legend
fig.savefig('20131108T0422.pdf',dpi=150,bbox_inches='tight')
plt.show()
