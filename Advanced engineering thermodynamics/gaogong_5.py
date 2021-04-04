import math

R=8.31451
"""用字典、列表记录各个制冷剂的热物性参数"""
dic_refri = {0:'R290',1:'R600a',2:'R1234yf',3:'R1234ze',4:'R290|R600a 50%|50%'}
list_M = [44.096e-3,58.122e-3,114.04e-3,114.04e-3]
list_Tc = [369.89,407.81,367.85,382.52]
list_pc = [4.2512,3.629,3.3822,3.6363]
list_w = [0.1521,0.184,0.276,0.313]

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

def calc_mix(T,p,Tc1,pc1,M1,w1,Tc2,pc2,M2,w2):
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
    Z = 1.1 ; Z = NewtonIteration(A, B, Z)
    v = Z * R * T / p
    v1 = R * T / p
    mf1 = math.exp(b1/b*(Z-1) - math.log(Z-B) - A/(2*1.1412*B)*(2*(x1*a1+x2*(1-k12)*math.sqrt(a1*a2))/a - b1/b)*math.log((Z+2.414*B)/(Z-0.414*B)))
    mf2 = math.exp(b2/b*(Z-1) - math.log(Z-B) - A/(2*1.1412*B)*(2*(x2*a2+x1*(1-k12)*math.sqrt(a1*a2))/a - b2/b)*math.log((Z+2.414*B)/(Z-0.414*B)))
    lnfg = R*T*math.log((v-b)/v) + R*T*math.log(v/v1) - a/(2*1.4142*b)*math.log((v-0.414*b)/(v+2.414*b)) - (1-Z)
    return lnfg

T = float(input('请输入温度T = ？ ℃\n'))
p = 1.013e5 ; T = T+273.15
lnfg = calc_mix(T,p,list_Tc[0], list_pc[0] * 1000000, list_M[0], list_w[0],list_Tc[1], list_pc[1] * 1000000, list_M[1], list_w[1])
print(f'气相逸度系数的自然对数为： {lnfg},气相逸度系数为： {math.exp(lnfg)}')

