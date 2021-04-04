import math
import numpy as np
import matplotlib.pyplot as plt

R=8.31451
"""用字典、列表记录各个制冷剂的热物性参数"""

"""
dic_refri = {0:'R290',1:'R600a',2:'R1234yf',3:'R1234ze',4:'R290|R600a 50%|50%'}
list_M = [44.096e-3,58.122e-3,114.04e-3,114.04e-3]
list_Tc = [351.255,407.81,367.85,382.52]
list_pc = [5.78,3.629,3.3822,3.6363]
list_w = [0.277,0.184,0.276,0.313]

"""

dic_refri = {0:'R290',1:'R600a',2:'R1234yf',3:'R1234ze'}
list_M = [44.096e-3,58.122e-3,114.04e-3,114.04e-3]
list_Tc = [369.89,407.81,367.85,382.52]
list_pc = [4.2512,3.629,3.3822,3.6363]
list_w = [0.1521,0.184,0.276,0.313]

nprange = [[],[],[],[]]
x_plot = [[]]
y_plot = [[]]

def NewtonIteration(A,B,Z):
    """牛顿迭代法"""
    for i in range(0,100,1):
        Z0=Z
        f = Z * Z * Z - (1 - B) * Z * Z + (A - 3 * B * B - 2 * B) * Z - (A * B - B * B - B * B * B)
        df = 3 * Z * Z - 2 * (1 - B) * Z + (A - 3 * B * B - 2 * B)
        Z = Z - f / df
        if (abs(Z-Z0)<1e-7):  break
    return Z

def calc_ab(T,p,Tc,pc,M,w):
    """单质制冷剂的a、b计算"""
    k = 0.37464 + 1.54226 * w - 0.26992 * w * w
    Tr = T/Tc
    Alpha = pow((1+k*(1-pow(Tr,0.5))),2)
    a = 0.45727 * Alpha * R * R * Tc * Tc / pc
    b = 0.07780*R*Tc/pc
    return a,b,M


"""开始主程序"""
while True:
    i = int(input('请输入选择的制冷剂序号 0:R290,1:R600a,2:R1234yf,3:R1234ze\n'))
    if i <= 3 and i >= 0:
        print(f'您选择了{dic_refri[int(i)]}')
        break
    else:
        print('您输入有误，请输入i的范围在0-3之间')

print('p-T数据列表：')
print('P（MPa）     T（K）  ')
nprange[0] = np.arange(0.5,1.4,0.1) ; nprange[1] = np.arange(0.4,1.2,0.05)
nprange[2] = np.arange(0.8,1.45,0.05) ;nprange[3] = np.arange(0.6,1.4,0.1)
Trange = [270,300,300,300]

for p in nprange[i]:
    p=p*1000000
    T = Trange[i]
    while T<500 :
        ab_result = calc_ab( T, p, list_Tc[i], list_pc[i] * 1000000, list_M[i], list_w[i])
        a = ab_result[0] ; b = ab_result[1] ; M = ab_result[2]
        A = a * p / (R * R * T * T)
        B = b * p / (R * T)
        Z1 = NewtonIteration(A,B,1.1)
        Z2 = NewtonIteration(A,B,0.001)
        #print(Z1,Z2)
        faiv =math.exp(math.exp((Z1 - 1) - math.log(Z1 - B) - A / (2 * 1.414 * B) * math.log((Z1 + 2.414 * B) / (Z1 - 0.414 * B))))
        fail =math.exp(math.exp((Z2 - 1) - math.log(Z2 - B) - A / (2 * 1.414 * B) * math.log((Z2 + 2.414 * B) / (Z2- 0.414 * B))))
        if abs(faiv - fail) < 1e-4 :  break
        T += 0.001
    p = p/1000000
    print('%.2f    %.3f' %(p,T))


    x_plot[0].append(T)
    y_plot[0].append(p)
plt.plot(x_plot[0], y_plot[0], linewidth=2, color='#007500', label='   ')
plt.show()
