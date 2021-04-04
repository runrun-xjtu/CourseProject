import math
import numpy as np
import matplotlib.pyplot as plt

T_plot = [[]]
p_plot = [[]]
Tc = 647.14
pc = 22.064
h_range = [0,500,1000,2000,3000,4000,5000]
p_range = []
T_range = []

"""水 p-T 图"""
for T in np.arange(353.15,647.14,0.5):
    p = math.log(T/Tc) * (7.60794067 + 10.1932439*pow((1-T/Tc),1.89) + 21.1083545*pow((1-T/Tc),5.67))
    p = math.exp(p)
    p = pc * p

    T_plot[0].append(T)
    p_plot[0].append(p)

plt.plot(T_plot[0], p_plot[0], linewidth=2, color='#007500', label='log1.5(x)')
plt.show()

"""不同海拔的大气压和饱和温度计算"""
print('h(m)     p(MPa)    T(K)')
for h in h_range:
    p = math.exp(5.25885*math.log(288.15-0.0065*h)-18.2573)
    p_range.append(p)

for p in p_range:
    p = p / 1000000
    for i in range(0, 600, 1):
        if p >= p_plot[0][i] and p <= p_plot[0][i + 1]:
            T_range.append(T_plot[0][i])
            break
for i in range(7):
    print('%4d      %.3f     %.2f' %(h_range[i],p_range[i]/1000000,T_range[i]))

