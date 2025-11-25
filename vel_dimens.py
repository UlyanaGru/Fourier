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

filename = f"./data_vel/v2_udimnesless.csv"

data = pd.read_csv(filename, sep=';', usecols=[0,1,2,4])
Ugy = data.iloc[:,0]
Vof = data.iloc[:,1]
Udim = data.iloc[:,2]
CoordY = data.iloc[:,3]

init_matplotlib(plt)

plt.plot(CoordY, Ugy, 'r-', label='Профиль скорости')
plt.plot(CoordY, Udim, 'b-', label='Безразмерный профиль')
plt.xlabel('Координата Y по пластине')
plt.ylabel('Скорость')
plt.legend(loc='lower left')
plt.grid(True, alpha=0.3)
plt.title(f'Сравнение и безразмерного профиля скорости')
plt.show()