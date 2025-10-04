# -*- coding: cp1251 -*-

def init_matplotlib(name):
    import matplotlib.font_manager as font_manager

    name.rcParams.update(plt.rcParamsDefault)
    name.style.use('classic')
    prop = font_manager.FontProperties(fname='c:\\windows\\fonts\\times.ttf')
    name.rcParams['font.family'] = prop.get_name()
    name.rcParams['lines.linewidth'] = 0.5
    name.rcParams["figure.figsize"] = (140.0/25.4,70.0/25.4) #mm to inch
    name.rcParams["figure.dpi"] =150
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
from scipy.signal import correlate
#from labellines import *

init_matplotlib(plt)

filename = "./data/s2d_film_time_statistic.dat"
with open(filename, "r") as f:
    first_line = f.readline().strip()   # читаем первую строку
    dxmesh = float(first_line)

data = pd.read_csv("./data/s2d_film_time_statistic.dat",sep=',',skiprows=1)
last_column_name = data.columns[-1]

nmax = int(last_column_name[1:])
lmax = nmax*dxmesh


print(f"dxmesh = {dxmesh:.4g}")
print(f"NMax = {nmax:.4g}")
print(f"LMax = {lmax:.4g}")

def X2N(x,dx,nx):
    return int(x/dx*nx-0.5)
def N2X(n,dx):
    return n*dx
tshift = 0.5
#выбираем точки
n1 = 1000
xp1 = N2X(n1,dxmesh)
n2 = 1020
xp2 = N2X(n2,dxmesh)
print(f"delta21 = {xp2-xp1:.4g}")
name1 = "N" + str(n1)
name2 = "N" + str(n2)
tval = np.array(data["Time"])
delta1 = np.array(data[name1])
delta2 = np.array(data[name2])


patch2fileout = './/figout//timesignal_'

init_matplotlib(plt)

fig, ax = plt.subplots(nrows=1,ncols = 1)

nameX = '$\mathrm{Time,\ ms}$'
nameY = '$\\delta \mathrm{,\ mm}$'

ax1 = ax

ax1.plot(tval*1.0e3,delta1*1.0e3,'-b',label = "-1")
ax1.plot(tval*1.0e3,delta2*1.0e3,'--g',label = "-2")

ax1.set_xlim(xmin = 0.0)
ax1.set_ylim(ymin = 0.0)
ax1.set_xlabel(nameX)
ax1.set_ylabel(nameY)
#ax1.set_xticks(np.arange(xmin,xmax+dx,dx))
#ax1.set_yticks(np.arange(ymin,ymax+dy,dy))

plt.tight_layout()
plt.savefig(patch2fileout)
plt.show()
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
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
    
def dominant_wavelengths(thickness, dx, top_n=1, max_lambda=None):
    fs = 1.0/dx
    f, X = freq_plot(thickness, max_freq=fs/2, fs=fs)

    f = np.asarray(f)
    X = np.asarray(X)
    nonzero = f > 0
    f = f[nonzero]
    X = X[nonzero]

    lambdas = 1.0 / f

    if max_lambda is not None:
        mask = lambdas <= max_lambda
        lambdas = lambdas[mask]
        X = X[mask]

    idx_sorted = np.argsort(X)[::-1]
    top_idx = idx_sorted[:top_n]
    dominant_lambdas = lambdas[top_idx]
    dominant_amps = X[top_idx]

    return dominant_lambdas, dominant_amps, lambdas, X
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
patch2fileout = './/figout//freq_'

init_matplotlib(plt)

fig, ax = plt.subplots(nrows=1,ncols = 1)

nameX = '$f \mathrm{,\ Гц}$'
nameY = '$(\\delta - \\overline{\\delta})/\\overline{\\delta} \\cdot 100 \%$'

fs = 1000.0
max_freq = 100.0
xmin, xmax, dx = 0.0, max_freq, 5.0
ymin, ymax, dy = 0.0, 25.0, 5.0

ax1 = ax

x1, y1 = tval[tval>=tshift], delta1[tval>=tshift]
X,Y = freq_plot((y1-np.average(y1))/np.average(y1)*100, max_freq, fs)
ax1.stem(X, Y)

ax1.set_xlim(xmin = xmin, xmax=xmax)
ax1.set_ylim(ymin = 0.0)
ax1.set_xlabel(nameX)
ax1.set_ylabel(nameY)
ax1.set_xticks(np.arange(xmin,xmax+dx,dx))
#ax1.set_yticks(np.arange(ymin,ymax+dy,dy))

plt.tight_layout()
plt.savefig(patch2fileout)
plt.show()
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------    
timelook = 1.0
idx = (data["Time"] - timelook).abs().idxmin()
deltaxdistr = np.array(data.loc[idx])[1:]
xval = np.linspace(dxmesh/2, dxmesh/2 + (nmax-1)*dxmesh, nmax)

lambdadominant = dominant_wavelengths(deltaxdistr, dxmesh, top_n=1, max_lambda=None)[0][0]

patch2fileout = './/figout//xsignal_'

init_matplotlib(plt)

fig, ax = plt.subplots(nrows=1,ncols = 1)

nameX = '$\mathrm{Distance,\ mm}$'
nameY = '$\\delta \mathrm{,\ mm}$'

ax1 = ax

ax1.plot(xval*1.0e3,deltaxdistr*1.0e3,'-b',label = "-1")

ax1.set_xlim(xmin = 0.0)
ax1.set_ylim(ymin = 0.0)
ax1.set_xlabel(nameX)
ax1.set_ylabel(nameY)
ax1.set_title(f"Основная длина волны: {lambdadominant*1.0e3:.3g} mm")
#ax1.set_xticks(np.arange(xmin,xmax+dx,dx))
#ax1.set_yticks(np.arange(ymin,ymax+dy,dy))

plt.tight_layout()
plt.savefig(patch2fileout)
plt.show()

#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
y2 = delta2[tval>=tshift]
corr = correlate(y1, y2, mode='full')
lags = np.arange(-len(y1)+1, len(y2))

lag_max = lags[np.argmax(corr)]
time_shift = lag_max / fs

print(f"Сдвиг между сигналами: {lag_max} отсчётов, или {time_shift:.6f} s")
uw = abs(xp2-xp1)*1.0e3/abs(time_shift)
print(f"Фазовая скорость: {uw:.3g} mm/s")
patch2fileout = './/figout//corr_'

fig, ax = plt.subplots(nrows=1,ncols = 1)

ax1 = ax

ax1.plot(lags/fs, corr)
ax1.axvline(time_shift, color='r', linestyle='--', label=f"Сдвиг = {time_shift:.6f} c")
ax1.set_xlabel("Лаг (сек)")
ax1.set_ylabel("Корреляция")
ax1.legend()
ax1.set_title("Взаимная корреляция")

plt.tight_layout()
plt.savefig(patch2fileout)
plt.show()
