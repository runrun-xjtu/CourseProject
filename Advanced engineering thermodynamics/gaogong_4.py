import math

R=8.31451
"""用字典、列表记录各个制冷剂的热物性参数"""
dic_refri = {0:'R290',1:'R600a',2:'R1234yf',3:'R1234ze',4:'R290|R600a 50%|50%'}
list_M = [44.096e-3,58.122e-3,114.04e-3,114.04e-3]
list_Tc = [369.89,407.81,367.85,382.52]
list_pc = [4.2512,3.629,3.3822,3.6363]
list_w = [0.1521,0.184,0.276,0.313]
list_p0 = [0.47446,0.15696,0.31582,0.21648,0.32979]
list_c0 = [-95.80,-23.91,18.349,55.389]
list_c1 = [6.945,6.605,128.316,10.784]
list_c2 = [-3.597e-3,-3.176e-3,-33.354,99.25]
list_c3 = [7.29e-7,4.981e-7,2.086,-49.88]

def NewtonIteration(A,B,Z):
    """牛顿迭代法"""
    for i in range(0,100,1):
        Z0=Z
        f = Z * Z * Z - (1 - B) * Z * Z + (A - 3 * B * B - 2 * B) * Z - (A * B - B * B - B * B * B)
        df = 3 * Z * Z - 2 * (1 - B) * Z + (A - 3 * B * B - 2 * B)
        Z = Z - f / df
        if (abs(Z-Z0)<1e-5)and(f<1e-5): break
    return Z         #这次不返回i

def calc_ab(T,p,Tc,pc,M,w):
    """单质制冷剂的a、b计算"""
    k = 0.37464 + 1.54226 * w - 0.26992 * w * w
    Tr = T/Tc
    Alpha = pow((1+k*(1-pow(Tr,0.5))),2)
    a = 0.45727 * Alpha * R * R * Tc * Tc / pc
    b = 0.07780*R*Tc/pc
    return a,b,M

def calc_arsr(T,v,v1,a,b,beta):
    """ar sr计算"""
    ar = R*T*math.log((v-b)/v) - a/(2*1.4142*b)*math.log((v-0.414*b)/(v+2.414*b)) + R*T*math.log(v/v1)
    sr = -1*R*math.log((v-b)/v) + beta/(2*1.4142*b)*math.log((v-0.414*b)/(v+2.414*b)) - R*math.log(v/v1)
    return ar,sr

def calc_hr(T,p,Tc,pc,M,w):
    """hr计算"""
    ab_result = calc_ab(T,p,Tc,pc,M,w)
    a = ab_result[0] ; b = ab_result[1] ; M = ab_result[2]

    ab_result_calcbeta = calc_ab(T+0.1, p, Tc, pc, M, w) ; a_calcbeta = ab_result_calcbeta[0]
    beta = (a_calcbeta-a) / 0.1

    A = a * p / (R * R * T * T)
    B = b * p / (R * T)
    Z = 1.1 ; Z = NewtonIteration(A,B,Z)

    v = Z * R * T / p
    v1= R * T / p

    arsr_result = calc_arsr(T,v,v1,a,b,beta)
    ar = arsr_result[0]
    sr = arsr_result[1]
    hr = ar + T*sr + R*T*(1-Z)
    hr = hr / M /1000
    sr = sr / M /1000
    return hr,sr,M

def calc_abZmix(T,p,Tc1,pc1,M1,w1,Tc2,pc2,M2,w2):
    """计算混合物的a，b，A，B，Z"""
    ab_result1 = calc_ab(T,p,Tc1,pc1,M1,w1)
    a1 = ab_result1[0] ; b1 = ab_result1[1] ; M1 = ab_result1[2]
    ab_result2 = calc_ab(T,p,Tc2,pc2,M2,w2)
    a2 = ab_result2[0] ; b2 = ab_result2[1] ; M2 = ab_result2[2]

    k12 = 0.01
    x1 = 1 / (1 + M1/M2) ; x2 = 1 / (1 + M2/M1)
    a = x1 * x1 * a1 + x2 * x2 * a2 + 2 * x1 * x2 * (1 - k12) * math.sqrt(a1 * a2)
    b = x1 * b1 + x2 * b2
    M = x1 * M1 + x2 * M2
    A = a * p / (R * R * T * T)
    B = b * p / (R * T)
    Z = 1.1 ; Z = NewtonIteration(A,B,Z)
    return a,b,Z,M

def calc_hrmix(T,p,Tc1,pc1,M1,w1,Tc2,pc2,M2,w2):
    """计算混合物的beta，hr，sr"""
    ab_result = calc_abZmix(T,p,Tc1,pc1,M1,w1,Tc2,pc2,M2,w2)
    a = ab_result[0] ; b = ab_result[1] ; Z = ab_result[2] ; M = ab_result[3]
    ab_result_calcbeta = calc_abZmix(T+0.1, p, Tc2, pc2, M2, w2) ; a_calcbeta = ab_result_calcbeta[0]
    beta = (a_calcbeta-a) / 0.1

    v = Z * R * T / p
    v1= R * T / p
    arsr_result = calc_arsr(T,v,v1,a,b,beta)
    ar = arsr_result[0]
    sr = arsr_result[1]
    hr = ar + T*sr + R*T*(1-Z)
    hr = hr / M
    sr = sr / M
    return hr,sr,M

"""开始主程序"""
while True:
    i = int(input('请输入选择的制冷剂序号 0:R290,1:R600a,2:R1234yf,3:R1234ze,4:R290|R600a\n'))
    if i <= 4 and i >= 0:
        print(f'您选择了{dic_refri[int(i)]}')
        break
    else:
        print('您输入有误，请输入i的范围在0-4之间')

T = float(input('请输入温度T = ？ ℃\n'))
p = float(input('请输入压力p = ？ MPa\n'))
T=T+273.15 ; p=p*1000000 ; T0 = 273.15

if i<=3:
    hrsr_result = calc_hr(T,p,list_Tc[i], list_pc[i] * 1000000, list_M[i], list_w[i])
    hr = hrsr_result[0] ; sr = hrsr_result[1] ; M = hrsr_result[2]
    print(f'该温度及压力下的该工质的余焓为hr=%.6f kJ/kg' %hr)
    print(f'该温度及压力下的该工质的余熵为sr=%.6f kJ/kgK' %sr)
    hrsr_result0 = calc_hr(273.15,list_p0[i]*1000000,list_Tc[i], list_pc[i] * 1000000, list_M[i], list_w[i])
    hr0 = hrsr_result0[0] ; sr0 = hrsr_result0[1]
    print(f'0℃及对应压力下的该工质的余焓为hr=%.6f kJ/kg' %hr0)
    print(f'0℃及对应压力下的该工质的余熵为sr=%.6f kJ/kgK' %sr0)

    integral_h = list_c0[i]*T + 1/2*list_c1[i]*pow(T,2) + 1/3*list_c2[i]*pow(T,3) + 1/4*list_c3[i]*\
                 pow(T,4) - (list_c0[i]*T0+1/2*list_c1[i]*pow(T0,2)+1/3*list_c2[i]*pow(T0,3)+1/4*list_c3[i]*pow(T0,4))
    integral_s = list_c1[i]*T + 1/2*list_c2[i]*pow(T,2) + 1/3*list_c3[i]*pow(T,3) - list_c0[i]/pow(T,2) - \
                 (list_c1[i]*T0+1/2*list_c2[i]*pow(T0,2)+1/3*list_c3[i]*pow(T0,3)-list_c0[i]/pow(T0,2))
    h = 200*M + hr0*M + integral_h  - hr*M
    s = 1*M + sr0*M + integral_s  - R * math.log(p/list_p0[i]*1000000)  - sr*M
    h = h / M /1000 ; s = s / M /1000
else:
    x1 = 1 / (1 + list_M[0]/list_M[1]) ; x2 = 1 / (1 + list_M[1]/list_M[0])
    c0 = x1 * list_c0[0] + x2 * list_c0[1]
    c1 = x1 * list_c1[0] + x2 * list_c1[1]
    c2 = x1 * list_c2[0] + x2 * list_c2[1]
    c3 = x1 * list_c3[0] + x2 * list_c3[1]

    hrsr_result = calc_hrmix(T,p,list_Tc[0],list_pc[0]*1000000,list_M[0],list_w[0],list_Tc[1],list_pc[1]*1000000,list_M[1],list_w[1])
    hr = hrsr_result[0] ; sr = hrsr_result[1] ; M = hrsr_result[2]
    print(f'该温度及压力下的该工质的余焓为hr=%.6f kJ/kg' %hr)
    print(f'该温度及压力下的该工质的余熵为sr=%.6f kJ/kgK' %sr)
    hrsr_result0 = calc_hrmix(T0,list_p0[4]*1000000,list_Tc[0],list_pc[0]*1000000,list_M[0],list_w[0],list_Tc[1],list_pc[1]*1000000,list_M[1],list_w[1])
    hr0 = hrsr_result0[0] ; sr0 = hrsr_result0[1]
    print(f'0℃及对应压力下的该工质的余焓为hr=%.6f kJ/kg' %hr0)
    print(f'0℃及对应压力下的该工质的余熵为sr=%.6f kJ/kgK' %sr0)

    integral_h = c0*T + 1/2*c1*pow(T,2) + 1/3*c2*pow(T,3) + 1/4*c3*pow(T,4) - (c0*T0+1/2*c1*pow(T0,2)+1/3*c2*pow(T0,3)+1/4*c3*pow(T0,4))
    integral_s = c1*T + 1/2*c2*pow(T,2) + 1/3*c3*pow(T,3) - c0/pow(T,2) - (c1*T0+1/2*c2*pow(T0,2)+1/3*c3*pow(T0,3)-c0/pow(T0,2))
    h = 200*M + hr0*M + R*integral_h  - hr*M
    s = 1*M + sr0*M + R*integral_s  - R * math.log(p/(list_p0[i]*1000000))  - sr*M
    h = h/M/1000 ; s =s/M/1000

print(f'该温度及压力下工质的焓为h=%.6f kJ/kg' %h)
print(f'该温度及压力下工质的熵为s=%.6f kJ/kgK' %s)


























