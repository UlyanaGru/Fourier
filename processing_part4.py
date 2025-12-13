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

#----Поиск максимальной амлитуды в каждой точке
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import pandas as pd

nx_max = 1500
ny_max = 60
v = 6
dt_step = 0.001
x_frac_start = 0.3 
x_frac_end = 1.0
t_start = 0.7

filename = f"./data_variousG_6-12/s2d_film_time_statistic_v{v}.dat"
with open(filename, "r") as f:
    first_line = f.readline().strip()
    dx_step = float(first_line)
tf_ind = int(t_start/dt_step)
xf_ind = int(x_frac_start*nx_max)
xl_ind = int(x_frac_end*nx_max)
data = np.loadtxt(filename, skiprows=2, delimiter=',')
""" data = np.genfromtxt(
    filename, 
    skip_header=2,
    delimiter=',', 
    usecols=range(1314),
    missing_values='')
data_time = data[:, 0]
data = data[:, 1:]*1.0e3
index = np.argmax(data_time[data_time>0.2]) """

#-Для последнего момента времени

time_slice = -1  #последний момент времени
data_slice = data[time_slice, :]  #массив по пространству
peaks_ind = []
peaks_values = []
#находим пики по пространственной координате
peaks, properties = find_peaks(data_slice)
peaks_ind.append(peaks)
peaks_values = data_slice[peaks]
x_coords = np.arange(0, dx_step * len(data_slice), dx_step)
#print(f"Размер x_coords: {x_coords.shape}")
#print(f"Размер data_slice: {data_slice.shape}")
plt.plot(x_coords, data_slice, 'r-', label='Толщина пленки, мм')
plt.plot(x_coords[peaks], peaks_values, 'xb', label='Пики')
plt.xlabel('Координата по пластине')
plt.ylabel('Толщина пленки')
plt.legend()
plt.grid(True, alpha=0.3)
plt.title(f'Пики толщины пленки в последний момент времени')
plt.show()

#-Для всех времен
data = np.loadtxt(filename, skiprows=2, delimiter=',')
""" data = np.genfromtxt(
    filename, 
    skip_header=2,
    delimiter=',', 
    usecols=range(1314),
    missing_values='') """
data = data[:, 1:]*1.0e3 
data_amlitude = np.copy(data[tf_ind:, :])
max_val = np.max(data_amlitude, axis=0, keepdims=False)
min_val = np.min(data_amlitude, axis=0, keepdims=False)
amlitude = (max_val - min_val)/(max_val + min_val)
x_coords = np.arange(0,amlitude.shape[0],1)
plt.plot(x_coords, amlitude, 'b-')
plt.xlabel('Длина, КО', fontsize=18)
plt.ylabel(r'$\alpha = \frac{\alpha_{\max} - \alpha_{\min}}{\alpha_{\max} + \alpha_{\min}}$', fontsize=18)
plt.grid(True, alpha=0.3)
plt.title(f'Амплитуда волн, б/м', fontsize=20)
plt.show()

print(dx_step)