import math
import numpy as np
import matplotlib.pyplot as plt

R=8.31451
"""用字典、列表记录各个制冷剂的热物性参数"""
dic_refri = {0:'R290',1:'R600a'}
list_M = [44.096e-3,58.122e-3]
list_Tc = [369.89,407.81]
list_pc = [4.2512,3.629]
list_w = [0.1521,0.184]
x_plot = [[],[]]
x2_plot = [[],[]]
y_plot = [[],[]]
Tm = [-60+273.15,273.15]

k12 = 0.01 ; T1 = 200 ; T2 = 260

def NewtonIteration(A,B,Z):
    """牛顿迭代法"""
    for i in range(0,100,1):
        Z0=Z
        f = Z * Z * Z - (1 - B) * Z * Z + (A - 3 * B * B - 2 * B) * Z - (A * B - B * B - B * B * B)
        df = 3 * Z * Z - 2 * (1 - B) * Z + (A - 3 * B * B - 2 * B)
        Z = Z - f / df
        if (abs(Z-Z0)<1e-5)and(f<1e-5): break
    return Z

def calc_ab(T,p,Tc,pc,M,w):
    """单质制冷剂的a、b计算"""
    k = 0.37464 + 1.54226 * w - 0.26992 * w * w
    Tr = T/Tc
    Alpha = pow((1+k*(1-pow(Tr,0.5))),2)
    a = 0.45727 * Alpha * R * R * Tc * Tc / pc
    b = 0.07780*R*Tc/pc
    return a,b,M

def calc_abmix(T,p,Tc1,pc1,M1,w1,Tc2,pc2,M2,w2,Z,y1,y2):
    """混合制冷剂的a、b计算"""
    k12 = 0.01
    ab1 = calc_ab(T,p,Tc1,pc1,M1,w1)
    ab2 = calc_ab(T,p,Tc2,pc2,M2,w2)
    a1 = ab1[0] ; b1 =  ab1[1]
    a2 = ab2[0] ; b2 =  ab2[1]
    x1 = 1 / (1 + M1/M2) ; x2 = 1 / (1 + M2/M1)
    a = x1 * x1 * a1 + x2 * x2 * a2 + 2 * x1 * x2 * (1 - k12) * math.sqrt(a1 * a2)
    b = x1 * b1 + x2 * b2
    M = x1 * M1 + x2 * M2
    A = a * p / (R * R * T * T)
    B = b * p / (R * T)
    Z = NewtonIteration(A, B, Z)
    mf1 = math.exp(b1/b*(Z-1) - math.log(Z-B) - A/(2*1.1412*B)*(2*(y1*a1+y2*(1-k12)*math.sqrt(a1*a2))/a - b1/b)*math.log((Z+2.414*B)/(Z-0.414*B)))
    mf2 = math.exp(b2/b*(Z-1) - math.log(Z-B) - A/(2*1.1412*B)*(2*(y2*a2+y1*(1-k12)*math.sqrt(a1*a2))/a - b2/b)*math.log((Z+2.414*B)/(Z-0.414*B)))
    return mf1,mf2

"""开始主程序"""
for p in [1.013e5,10*1.013e5]:
    if p == 1.013e5 : m = 0
    else : m =1
    print(f'p={p/1.013e5} atm')
    print(f'x   y   T')

    y1 = 0
    while y1<=1:
        y2 = 1 - y1 ; x1 = 0.1
        x2 = 1 - x1 ; x = 0
        T = Tm[m]
        while abs(x - 1) >= 0.01:
            T = T + 0.1
            result1 = calc_abmix(T,p,list_Tc[0],list_pc[0]*1000000,list_M[0],list_w[0],list_Tc[1],list_pc[1]*1000000,list_M[1],list_w[1],0.001,x1,x2)
            result2 = calc_abmix(T,p,list_Tc[0],list_pc[0]*1000000,list_M[0],list_w[0],list_Tc[1],list_pc[1]*1000000,list_M[1],list_w[1],1.001,y1,y2)
            mfl1 = result1[0] ; mfl2 = result1[1]
            mfg1 = result2[0] ; mfg2 = result2[1]

            k1 = mfg1 / mfl1
            k2 = mfg2 / mfl2
            x1 = k1 * y1 / (k1 * y1 + k2 * y2)
            x2 = k2 * y2 / (k1 * y1 + k2 * y2)
            x0 = x
            x = k1 * y1 + k2 * y2
            while abs(x - x0) > 0.001:
                result1 = calc_abmix(T,p,list_Tc[0],list_pc[0]*1000000,list_M[0],list_w[0],list_Tc[1],list_pc[1]*1000000,list_M[1],list_w[1],0.001,x1,x2)
                mfl1 = result1[0] ; mfl2 = result1[1]
                k1 = mfg1 / mfl1
                k2 = mfg2 / mfl2
                x1 = k1 * y1 / (k1 * y1 + k2 * y2)
                x2 = k2 * y2 / (k1 * y1 + k2 * y2)
                x0 = x
                x = k1 * y1 + k2 * y2
        print(f'{x1}     {y1}     {T}')

        x_plot[m].append(x1)
        x2_plot[m].append(y1)
        y_plot[m].append(T)
        y1+=0.01
    plt.plot(x_plot[m], y_plot[m], linewidth=2, color='#007500', label='   ')
    plt.plot(x2_plot[m], y_plot[m], linewidth=2, color='#9F35FF', label='   ')
    plt.show()
