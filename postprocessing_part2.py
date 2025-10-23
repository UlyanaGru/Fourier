import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.signal import correlate
from numpy.fft import fft2, fftshift, fftfreq

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

nx_max = 1500
ny_max = 30
dt_step = 0.001
x_frac_start = 0.3 
x_frac_end = 0.8
t_start = 0.7

filename = f"./data/s2d_film_time_statistic_{nx_max}_{ny_max}.dat"
with open(filename, "r") as f:
    first_line = f.readline().strip()
    dx_step = float(first_line)

tf_ind = int(t_start/dt_step)
xf_ind = int(x_frac_start*nx_max)
xl_ind = int(x_frac_end*nx_max)
data = np.loadtxt(filename, skiprows=2, delimiter=',')

with open(filename, 'r') as f:
    first_line = f.readline().strip() #1st strip
    #for structured mesh read width CV dx
    dxmesh = float(first_line)
last_column_name = data.columns[-1]
nmax = int(last_column_name[1:])
lmax = nmax*dxmesh
print(f"dxmesh = {dxmesh:.4g}")
print(f"NMax = {nmax:.4g}")
print(f"LMax = {lmax:.4g}")

def N2X(n,dx):
    return n*dx
#choose points
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

data = data[:, 1:]*1.0e3 
data = data[tf_ind:, xf_ind:xl_ind]
data_detrend = data - np.mean(data, axis=1, keepdims=True)
nt, nx = data_detrend.shape

# --- 2D FFT ---
win_t = np.hanning(nt)[:, None]
win_x = np.hanning(nx)[None, :]
data_win = data_detrend * win_t * win_x
spec2 = fftshift(fft2(data_win))
power = np.abs(spec2)**2 
freqs = fftshift(fftfreq(nt, d=1.0))     # шт / t-ед.
wavenums = fftshift(fftfreq(nx, d=1.0))

nameX = r'$\mathrm{Время,\ мс}$'
nameY = r'$\delta \mathrm{,\ мм}$'
fig, ax = plt.subplots(nrows=1,ncols = 1)
ax1 = ax
ax1.plot(tval*1.0e3,delta1*1.0e3,'-b',label = f"{xp1*1e3:.2g} мм")
ax1.set_xlim(xmin = 0.0)
ax1.set_ylim(ymin = 0.0)
ax1.set_xlabel(nameX)
ax1.set_ylabel(nameY)
ax1.legend(loc='lower right')
plt.tight_layout()
plt.show()