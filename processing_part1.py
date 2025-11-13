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
from numpy.fft import fft2, fftshift, fftfreq

nx_max = 1500
ny_max = 60
v = 5
dt_step = 0.001
x_frac_start = 0.3 
x_frac_end = 0.8
t_start = 0.7

filename = f"./data/s2d_film_time_statistic_{nx_max}_{ny_max}_{v}.dat"
with open(filename, "r") as f:
    first_line = f.readline().strip()
    dx_step = float(first_line)

tf_ind = int(t_start/dt_step)
xf_ind = int(x_frac_start*nx_max)
xl_ind = int(x_frac_end*nx_max)
data = np.loadtxt(filename, skiprows=2, delimiter=',', usecols=range(1,1500))
data = data[:, 1:]*1.0e3 

fig, ax = plt.subplots(nrows=1,ncols = 1)
ax1 = ax
xmesh = np.linspace(0, 114.770, data.shape[1]) 
tmesh = np.linspace(0, data.shape[0]*0.001, data.shape[0])
im1 = ax1.pcolor(xmesh,tmesh,data, cmap='Greys') 
plt.colorbar(im1, label='$\\delta \mathrm{,\ мм}$')
ax1.set_xlabel("Длина, мм")
ax1.set_ylabel("Время, с")
plt.tight_layout()
plt.show()

data = data[tf_ind:, xf_ind:xl_ind]
data_detrend = data - np.mean(data, axis=1, keepdims=True)
nt, nx = data_detrend.shape

fig, ax = plt.subplots(nrows=1,ncols = 1)
ax1 = ax
im1 = ax1.pcolor(xmesh[xf_ind:xl_ind],tmesh[tf_ind:],data_detrend, cmap='Greys') 
plt.colorbar(im1, label='$\\delta - \\overline{\\delta}\mathrm{,\ мм}$')
ax1.set_xlabel("Длина, мм")
ax1.set_ylabel("Время, с")
plt.tight_layout()
plt.show()

# --- 2D FFT ---
win_t = np.hanning(nt)[:, None]
win_x = np.hanning(nx)[None, :]
data_win = data_detrend * win_t * win_x

spec2 = fftshift(fft2(data_win))
power = np.abs(spec2)**2 
freqs = fftshift(fftfreq(nt, d=1.0))     # шт / t-ед.
wavenums = fftshift(fftfreq(nx, d=1.0))  # шт / x-ед

# --- Поиск пика ---
power_flat = power.copy()
peak_idx = np.unravel_index(np.argmax(power_flat), power_flat.shape)
f_peak = freqs[peak_idx[0]]
k_peak = wavenums[peak_idx[1]]

# --- Фазовая скорость ---
if k_peak != 0:
    u_phase = f_peak / k_peak   # в единицах (пространственных точек / временной шаг)
else:
    u_phase = np.nan

#print(f"Пик: f = {f_peak:.5f}, k = {k_peak:.5f}, u = {u_phase:.5f} точек / шаг по времени")

lambda_phys = 1 / (k_peak / dx_step)       # длина волны, м
freq_phys = f_peak / dt_step                     # частота, Гц
u_phys = freq_phys * lambda_phys            # фазовая скорость, м/с
print(f"Длина волны: {np.abs(lambda_phys*100.0):.4f} см")
print(f"Частота: {np.abs(freq_phys):.4f} Гц")
print(f"Фазовая скорость: {np.abs(u_phys*100.0):.4f} см/с")
print(f"Отклонение фазвой скорости от экспериментальной: {np.abs((u_phys*100.0 - 21.7)/21.7):.4f} %")

fig, ax = plt.subplots(nrows=1,ncols = 1)
ax1 = ax
im1 = ax1.imshow(np.log10(power + 1e-12), cmap='jet', aspect='auto', origin='lower', 
                    extent=(wavenums[0], wavenums[-1], freqs[0], freqs[-1]))
plt.colorbar(im1, label=r'$\mathrm{Мощность (log10)}$')
ax1.set_xlabel("Волновое число k, шт / x-ед.")
ax1.set_ylabel("Частота f, шт / t-ед.")
plt.tight_layout()
plt.show()

#Попытка 2D-анализа для сигнала в каждой отдельной точке с выводом на график
#список с значениям получившихся длин волн в точках от 0 до 1500
lenght_waves = list()
#список с значениям получившихся частот в точках от 0 до 1500
freq_waves = list()
#список с значениям получившихся фазовых скоростей в точках от 0 до 1500
uw_waves = list()

for x_id in range(nx_max):
    end_x = min(nx_max, x_id+1)
    start_x = max(0,x_id)
    data_detrend = data - np.mean(data, axis=1, keepdims=True)
    nt, nx = data_detrend.shape
    #тут 2D FFT
    win_t = np.hanning(nt)[:, None]
    win_x = np.hanning(nx)[None, :]
    data_win = data_detrend * win_t * win_x
    spec2 = fftshift(fft2(data_win))
    power = np.abs(spec2)**2 
    freqs = fftshift(fftfreq(nt, d=1.0))     # шт / t-ед.
    wavenums = fftshift(fftfreq(nx, d=1.0))  # шт / x-ед
    power_flat = power.copy()
    peak_idx = np.unravel_index(np.argmax(power_flat), power_flat.shape)
    f_peak = freqs[peak_idx[0]]
    k_peak = wavenums[peak_idx[1]]
    if k_peak != 0:
        u_phase = f_peak / k_peak
    else:
        u_phase = np.nan
    lambda_phys = 1 / (k_peak / dx_step)       # длина волны, м
    freq_phys = f_peak / dt_step                     # частота, Гц
    u_phys = freq_phys * lambda_phys            # фазовая скорость, м/с
    #обработка
    lenght_waves.append(np.abs(lambda_phys*100) if not np.isnan(lambda_phys) else 0)
    freq_waves.append(np.abs(freq_phys) if not np.isnan(freq_phys) else 0)
    uw_waves.append(np.abs(u_phys*100) if not np.isnan(u_phys) else 0)

lenght_waves = np.array(lenght_waves)
freq_waves = np.array(freq_waves)
uw_waves = np.array(uw_waves)

ymin0, ymax0 = 0.815, 0.825
ymin1, ymax1 = 27.0, 27.1
ymin2, ymax2 = 21., 23.

# fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
# axes[0].plot(range(nx_max), lenght_waves, 'b-', linewidth=1)
# axes[0].set_ylabel('Длина волны, см')
# axes[0].set_ylim(ymin0, ymax0)
# axes[0].grid(True)
# axes[1].plot(range(nx_max), freq_waves, 'r-', linewidth=1)
# axes[1].set_ylabel('Частота, Гц')
# axes[1].set_ylim(ymin1, ymax1)
# axes[1].grid(True)
# axes[2].plot(range(nx_max), uw_waves, 'g-', linewidth=1)
# axes[2].set_ylabel('Фазовая скорость, см/с')
# axes[2].set_xlabel('x-позиция')
# axes[2].set_ylim(ymin2, ymax2)
# axes[2].grid(True)
# plt.tight_layout()
# plt.show()