# -*- coding: cp1251 -*-

def init_matplotlib(name):
    import matplotlib.font_manager as font_manager

    name.rcParams.update(plt.rcParamsDefault)
    name.style.use('classic')
    prop = font_manager.FontProperties(fname='c:\\windows\\fonts\\times.ttf')
    name.rcParams['font.family'] = prop.get_name()
    name.rcParams['lines.linewidth'] = 0.5
    name.rcParams["figure.figsize"] = (140.0/25.4,150.0/25.4) #mm to inch
    name.rcParams["figure.dpi"] =500
    name.rcParams["savefig.dpi"] = plt.rcParams["figure.dpi"]
    name.rcParams["savefig.format"] = 'jpg'
    name.rcParams["font.size"] = 12
    name.rcParams["legend.fontsize"] = plt.rcParams["font.size"]
    name.rcParams["legend.numpoints"] = 1
    name.rcParams["axes.titlesize"] = plt.rcParams["font.size"]
    name.rcParams["axes.labelsize"] = plt.rcParams["font.size"]
    name.rcParams["xtick.labelsize"] = plt.rcParams["font.size"]
    name.rcParams["ytick.labelsize"] = plt.rcParams["font.size"]
    name.rcParams["legend.columnspacing"] = 1.5
    name.rcParams["markers.fillstyle"] = 'none'
    name.rcParams["lines.markersize"] = 3.0
    name.rcParams["legend.handlelength"] = 2.0
    name.rcParams["lines.markeredgewidth"] = 1.0
    name.rcParams["axes.formatter.limits"] = (-7,7)
    name.rcParams["axes.formatter.use_locale"] = True

def set_param_matplotlib(name,rcParamsName,rcParamsVal):
    """
    Изменение базовых настроек графиков
    """
    name.rcParams[rcParamsName] = rcParamsVal

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
#from labellines import *

init_matplotlib(plt)
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.hsv(np.linspace(0,1,10)))

data1 = pd.read_csv("./data/s2d_film_point1.dat",sep='\s+')
x1, y1 = np.array(data1["Time"]), np.array(data1["H"])
xp1 = np.arange(np.min(x1),np.max(x1),0.001)
yp1 = np.interp(xp1,x1,y1)
data2 = pd.read_csv("./data/s2d_film_point2.dat",sep='\s+')
x2, y2 = np.array(data2["Time"]), np.array(data2["H"])
xp2 = np.arange(np.min(x2),np.max(x2),0.001)
yp2 = np.interp(xp2,x2,y2)
data3 = pd.read_csv("./data/s2d_film_point3.dat",sep='\s+')
x3, y3 = np.array(data3["Time"]), np.array(data3["H"])
xp3 = np.arange(np.min(x3),np.max(x3),0.001)
yp3 = np.interp(xp3,x3,y3)

def freq_plot(sig, max_freq, fs):
    if(max_freq>fs/2):
        error('Max freq must be less than Nyquist frequency')
    X = np.abs(np.fft.fft(sig))
    N = np.size(sig)
     
    f = np.arange(0.0,fs+fs/N, fs/N, dtype=float)
    num_bins = np.count_nonzero(f> max_freq)
    X = X[0:num_bins]
    f = f[0:num_bins]
    return f, X/(N)*2

patch2fileout = './/figout//test_'

init_matplotlib(plt)
plt.rcParams["figure.figsize"] = (160.0/25.4,60.0/25.4)

fig, ax = plt.subplots(nrows=1,ncols = 1)

nameX = '$f \mathrm{,\ Гц}$'
nameY = '$(\\delta - \\overline{\\delta})/\\overline{\\delta} \\cdot 100 \%$'

fs = 1000.0
max_freq = 90.0
xmin, xmax, dx = 0.0, max_freq, 5.0
ymin, ymax, dy = 0.0, 25.0, 5.0

ax1 = ax

x0, y0 = xp1[xp1>=0.5], yp1[xp1>=0.5]
X,Y = freq_plot((y0-np.average(y0))/np.average(y0)*100, max_freq, fs)
ax1.stem(X, Y)

ax1.set_xlim(xmin = xmin, xmax=xmax)
ax1.set_ylim(ymin = 0.0, ymax=ymax)
ax1.set_xlabel(nameX)
ax1.set_ylabel(nameY)
ax1.set_xticks(np.arange(xmin,xmax+dx,dx))
#ax1.set_yticks(np.arange(ymin,ymax+dy,dy))

plt.tight_layout()
plt.savefig(patch2fileout)
plt.show()