import math
import numpy as np
import matplotlib.pyplot as plt

R=8.31451
"""用字典、列表记录各个制冷剂的热物性参数"""
dic_refri = {0:'R290',1:'R600a',2:'R1234yf',3:'R1234ze',4:'R290|R600a 50%|50%'}
list_M = [44.096e-3,58.122e-3,114.04e-3,114.04e-3]
list_Tc = [369.89,407.81,367.85,382.52]
list_pc = [4.2512,3.629,3.3822,3.6363]
list_w = [0.1521,0.184,0.276,0.313]
x_plot = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
x2_plot = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
y_plot = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
j = 0

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

def calc_abmix(T,p,Tc1,pc1,M1,w1,Tc2,pc2,M2,w2):
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
    return a,b,M

"""开始主程序"""

while True:
    Option = int(input('请选择功能：根据p,T求v（1）或pvT做图（2）\n'))
    #Option = 2
    if Option == 1 :
        while True:
            i = int(input('请输入选择的制冷剂序号 0:R290,1:R600a,2:R1234yf,3:R1234ze,4:R290|R600a\n'))
            if i <= 4 and i >= 0:
                print(f'您选择了{dic_refri[int(i)]}')
                break
            else:
                print('您输入有误，请输入i的范围在0-4之间')

        T = float(input('请输入温度T = ？ ℃\n'))
        p = float(input('请输入压力p = ？ MPa\n'))
        T=T+273.15 ; p=p*1000000

        if i<=3:
            ab_result = calc_ab(T,p,list_Tc[i],list_pc[i]*1000000,list_M[i],list_w[i])
            a = ab_result[0] ; b = ab_result[1] ; M = ab_result[2]
        else:
            ab_result = calc_abmix(T, p, list_Tc[0], list_pc[0] * 1000000, list_M[0], list_w[0],list_Tc[1], list_pc[1] * 1000000, list_M[1], list_w[1])
            a = ab_result[0] ; b = ab_result[1] ; M = ab_result[2]
        A = a * p / (R * R * T * T)
        B = b * p / (R * T)

        Z=1.1 ; Z1 = NewtonIteration(A,B,Z)
        Z=0.001 ; Z2 = NewtonIteration(A,B,Z)
        if abs(Z1-Z2)<1e-4:
            v = Z * R * T / p / M

            print('Z=%.6f  该温度压力下的工质的摩尔体积(m3/mol)为=%.6f' %(Z,v))
        else:
            vg = Z1 * R * T / p / M
            vl = Z2 * R * T / p / M
            print('Z1=%.6f  该温度压力下的气态工质的摩尔体积(m3/mol)为=%.6f' % (Z1,vg))
            print('Z2=%.6f  该温度压力下的液态工质的摩尔体积(m3/mol)为=%.6f' % (Z2,vl))
        break


    elif Option == 2:
        while True:
            i = int(input('请输入选择的制冷剂序号 0:R290,1:R600a,2:R1234yf,3:R1234ze,4:R290|R600a\n'))
            #i = 0
            if i <= 4 and i >= 0:
                print(f'您选择了{dic_refri[int(i)]}')
                break
            else:
                print('您输入有误，请输入i的范围在0-4之间')

        for T in np.arange(230,400,15):
            for p in np.arange(0.5,4,0.1):
                p = p * 1000000

                if i <= 3:
                    ab_result = calc_ab(T, p, list_Tc[i], list_pc[i] * 1000000, list_M[i], list_w[i])
                    a = ab_result[0] ; b = ab_result[1] ; M = ab_result[2]
                else:
                    ab_result = calc_abmix(T, p, list_Tc[0], list_pc[0] * 1000000, list_M[0], list_w[0], list_Tc[1],list_pc[1] * 1000000, list_M[1], list_w[1])
                    a = ab_result[0] ; b = ab_result[1] ; M = ab_result[2]
                A = a * p / (R * R * T * T)
                B = b * p / (R * T)
                Z=1.1 ; Z1 = NewtonIteration(A,B,Z)
                Z=0.001 ; Z2 = NewtonIteration(A,B,Z)
                v1 = Z1 * R * T / p / M
                v2 = Z2 * R * T / p / M
                x_plot[j].append(math.log2(v1))
                x2_plot[j].append(v2)
                y_plot[j].append(p/1000000)
                #print (f'T= %d K , p= %.4f MPa , v = %.6f m3/mol' %(T,p/1000000,v))

            plt.plot(x_plot[j], y_plot[j], linewidth=2, color='#007500', label='   ')
      #      plt.plot(x2_plot[j], y_plot[j], linewidth=2, color='#9F35FF', label='log1.5(x)')

            j += 1
        plt.show()
        break

    else:print('输入错误，请输入数字1或数字2')






