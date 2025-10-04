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
from scipy.signal import correlate
#from labellines import *

init_matplotlib(plt)
plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.hsv(np.linspace(0,1,10)))

data1 = pd.read_csv("./data/s2d_film_point1.dat",sep='\s+')
x1, y1 = np.array(data1["Time"]), np.array(data1["H"])
xp1 = np.arange(0.45000,np.max(x1),0.001)
yp1 = np.interp(xp1,x1,y1)
data2 = pd.read_csv("./data/s2d_film_point2.dat",sep='\s+')
x2, y2 = np.array(data2["Time"]), np.array(data2["H"])
xp2 = np.arange(0.45000,np.max(x2),0.001)
yp2 = np.interp(xp2,x2,y2)
data3 = pd.read_csv("./data/s2d_film_point3.dat",sep='\s+')
x3, y3 = np.array(data3["Time"]), np.array(data3["H"])
xp3 = np.arange(0.45000,np.max(x3),0.001)
yp3 = np.interp(xp3,x3,y3)

plt.plot(x1+0.014,y1,'-b')
plt.plot(x3,y3,'-g')
plt.show()
exit()
fs = 1000
x0 = yp2[xp2>=0.5]
y0 = yp3[xp3>=0.5]

x0 = x0 - np.mean(x0)
y0 = y0 - np.mean(y0)

plt.plot(x0,'-b')
plt.plot(y0,'--g')
#plt.plot(np.roll(y0,10))
plt.show()

corr = correlate(y0, x0, mode='full')
lags = np.arange(-len(x0)+1, len(y0)) 

lag_max = lags[np.argmax(corr)]
time_shift = lag_max / fs

print(f"Сдвиг между сигналами: {lag_max} отсчётов, или {time_shift:.6f} секунд")

# график корреляционной функции
plt.figure(figsize=(8,4))
plt.plot(lags/fs, corr)
plt.axvline(time_shift, color='r', linestyle='--', label=f"Сдвиг = {time_shift:.6f} c")
plt.xlabel("Лаг (сек)")
plt.ylabel("Корреляция")
plt.legend()
plt.title("Взаимная корреляция")
plt.grid()
plt.show()
