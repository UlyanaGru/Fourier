import matplotlib.pyplot as plt
import pandas as pd

HTCCalc_mean = pd.read_csv("HTCalc_mean.csv",sep=';',header=None)
HTCCalc_mean.columns = ["coordinate","HTC"]
HTCCalc_mean['coordinate'] = HTCCalc_mean['coordinate'].str.replace(',', '.').astype(float)
HTCCalc_mean['HTC'] = HTCCalc_mean['HTC'].str.replace(',', '.').astype(float)

HTCCalc_small = pd.read_csv("HTCalc_small.csv",sep=';',header=None)
HTCCalc_small.columns = ["coordinate","HTC"]
HTCCalc_small['coordinate'] = HTCCalc_small['coordinate'].str.replace(',', '.').astype(float)
HTCCalc_small['HTC'] = HTCCalc_small['HTC'].str.replace(',', '.').astype(float)

HTCexp = pd.read_csv("HTCexp.csv",sep=';',header=None)
HTCexp.columns = ["coordinate","HTC"]
HTCexp['coordinate'] = HTCexp['coordinate'].str.replace(',', '.').astype(float)
HTCexp['HTC'] = HTCexp['HTC'].str.replace(',', '.').astype(float)


fig, ax = plt.subplots()
htcmean, = ax.plot(HTCCalc_mean['coordinate'],HTCCalc_mean['HTC'], label='КТО средн.',linestyle='',markersize=7,
                   marker='D',markerfacecolor='white',markeredgecolor='black')
htcsmall, = ax.plot(HTCCalc_small['coordinate'],HTCCalc_small['HTC'], label ='КТО мал.',linestyle='',markersize=7,
                    marker='o',markerfacecolor='white',markeredgecolor='black')
htcexp, = ax.plot(HTCexp['coordinate'],HTCexp['HTC'], label = 'КТО эксп.',linestyle='',markersize=7,marker='s',
                  markerfacecolor='red',markeredgecolor='black')
ax.legend(handles=[htcexp,htcsmall,htcmean],fontsize=14,loc='upper left')
#ax.set_title("Зависимость коэффициента теплоотдачи от координаты",fontsize=18,pad=20)
ax.grid(True,alpha=0.4)
ax.set_xlim([0,1])
ax.set_ylim([0,4000])
ax.set_ylabel(r"Коэффициент теплоотдачи $\frac{Вт}{м^{2} \cdot К}$",fontsize=16)
ax.set_xlabel("Координата x",fontsize=16)
ax.annotate("R245fa\n" + r"G = 100 $\frac{кг}{м^{2} \cdot с}$",xy=(0.7,3500),fontsize=14)
ax.tick_params(labelsize=14)
plt.show()