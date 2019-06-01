# -*- coding: utf-8 -*-
"""
@author: 1510099204
"""

import numpy as np
import matplotlib.pyplot as plt

#パラメータ設定
a=0.6
gamma = 5
c=0.5
alpha_deg = 5
beta_deg = 20
v0=1.0

alpha_rad = alpha_deg*np.pi/180
beta_rad = beta_deg*np.pi/180

#Kutta condition
gamma = 4*np.pi*v0*a*np.sin(alpha_rad+beta_rad)

#プリントで言う「Z面」の作成
n = 256
rrr = np.linspace(a, a+5, n)
theta = np.linspace(0, 2*np.pi, n)

RRR,THETA = np.meshgrid(rrr,theta)

X = RRR * np.cos(THETA)
Y = RRR * np.sin(THETA)

Z_c = a * np.exp(1j*(np.pi-beta_rad))+c
Z = X + 1j*Y +Z_c

#プリントで言う「z面」の作成
z = (Z-Z_c)*np.exp(1j*(-alpha_rad))

#プリントで言う「ζ面」の作成
zeta = Z + c**2/Z

#ジュコフスキー翼表面の情報
rrr_c = np.ones_like(theta)*a
theta_c = theta

x_c, y_c = rrr_c*np.cos(theta_c),rrr_c*np.sin(theta_c)
z_c = x_c + 1j*y_c + Z_c

zeta_c = z_c + c**2/z_c


#=========以下計算==========

#流線関数psiの計算
comp = v0*(z + a**2/z) + 1j*gamma*np.log(z)/(2*np.pi)
psi = comp.imag

#Cpの計算
cz2 = X + 1j*Y + Z_c
cz1 = (cz2-Z_c)*np.exp(1j*(-alpha_rad))

cf = np.exp(-1j*alpha_rad)*v0*((1-a**2/cz1**2) + 1j*gamma/(2*np.pi*cz1))/(1.0-c**2/cz2**2)
Cp = 1.0 - (cf.real**2 + cf.imag**2)/v0**2

#壁面圧力分布の計算
cz2_c = z_c
cz1_c = (z_c-Z_c)*np.exp(1j*(-alpha_rad))

cf = np.exp(-1j*alpha_rad)*v0*((1-a**2/cz1_c**2) + 1j*gamma/(2*np.pi*cz1_c))/(1.0-c**2/cz2_c**2)
Cp_c = 1.0 - (cf.real**2 + cf.imag**2)/v0**2
#Cp_c = -Cp_c

#圧力係数の積分

cxp = 0.
cyp = 0.
x_m = np.linspace(-1,1,n)
cm = 0
for i in range(0,n-1):
    dxw = zeta_c[i+1].real-zeta_c[i].real
    dyw = zeta_c[i+1].imag-zeta_c[i].imag
    dnx = dyw
    dny = -dxw
    cpm = (Cp_c[i+1]+Cp_c[i])/2.0
    
    cxp = cxp - cpm*dnx
    cyp = cyp - cpm*dny
    cm = cm - (x_m.reshape(n,1) - zeta_c[i].real)*cpm*dny - zeta_c[i].imag*cpm*dnx

cxp = cxp/(4*c)
cyp = cyp/(4*c)
cm = cm/((4*c)**2)
    
cdp = cxp*np.cos(alpha_rad)+cyp*np.sin(alpha_rad)
clp = cyp*np.cos(alpha_rad)-cxp*np.sin(alpha_rad)

x_cp = np.argmin(abs(cm.reshape(1,n)))/float(n)

print '揚力係数 C_L={0}'.format(clp)
print '抵抗係数 C_D={0}'.format(cdp)
print '風圧中心 {0}%'.format(x_cp/(4*c)*100)

print cm.reshape(1,n)
print x_m

#流線の描画
plt.figure(1)
plt.contour(zeta.real, zeta.imag, psi, levels=np.arange(-10, 10, 0.1))
plt.plot(zeta_c.real, zeta_c.imag, color="black", linewidth=1.0, linestyle="-")
plt.xlim([-3.,3.])
plt.ylim([-3.,3.])

#Cpの分布の描画
plt.figure(2)
plt.contour(zeta.real, zeta.imag, Cp, levels=np.arange(-100, 100, 0.2))
plt.plot(zeta_c.real, zeta_c.imag, color="black", linewidth=1.0, linestyle="-")
plt.xlim([-3.,3.])
plt.ylim([-3.,3.])


plt.xticks(())
plt.yticks(())


#壁面圧力分布の描画
plt.figure(3)
plt.plot(zeta_c.real, zeta_c.imag, color="black", linewidth=1.0, linestyle="-")
plt.plot(zeta_c.real, Cp_c, color="green", linewidth=1.0, linestyle="-")
plt.xlim([-1.5,1.5])
plt.ylim([-4.,4.])

#Cmの計算（迎角変化）
m=6
alpha = np.linspace(-5,5,m)
for k in range(0,m):
    alpha_deg = alpha[k]
    alpha_rad = alpha_deg*np.pi/180

    #Kutta condition
    gamma = 4*np.pi*v0*a*np.sin(alpha_rad+beta_rad)
    
    #プリントで言う「Z面」の作成
    n = 256
    rrr = np.linspace(a, a+5, n)
    theta = np.linspace(0, 2*np.pi, n)
    
    RRR,THETA = np.meshgrid(rrr,theta)
    
    X = RRR * np.cos(THETA)
    Y = RRR * np.sin(THETA)
    
    Z_c = a * np.exp(1j*(np.pi-beta_rad))+c
    Z = X + 1j*Y +Z_c
    
    #プリントで言う「z面」の作成
    z = (Z-Z_c)*np.exp(1j*(-alpha_rad))
    
    #プリントで言う「ζ面」の作成
    zeta = Z + c**2/Z
    
    #ジュコフスキー翼表面の情報
    rrr_c = np.ones_like(theta)*a
    theta_c = theta
    
    x_c, y_c = rrr_c*np.cos(theta_c),rrr_c*np.sin(theta_c)
    z_c = x_c + 1j*y_c + Z_c
    
    zeta_c = z_c + c**2/z_c
    
    #壁面圧力分布の計算
    cz2_c = z_c
    cz1_c = (z_c-Z_c)*np.exp(1j*(-alpha_rad))
    
    cf = np.exp(-1j*alpha_rad)*v0*((1-a**2/cz1_c**2) + 1j*gamma/(2*np.pi*cz1_c))/(1.0-c**2/cz2_c**2)
    Cp_c = 1.0 - (cf.real**2 + cf.imag**2)/v0**2
    #Cp_c = -Cp_c
    
    #圧力係数の積分
    cxp = 0.
    cyp = 0.
    x_m = np.linspace(-1,1,n)
    cm = 0
    for i in range(0,n-1):
        dxw = zeta_c[i+1].real-zeta_c[i].real
        dyw = zeta_c[i+1].imag-zeta_c[i].imag
        dnx = dyw
        dny = -dxw
        cpm = (Cp_c[i+1]+Cp_c[i])/2.0
    
        cxp = cxp - cpm*dnx
        cyp = cyp - cpm*dny
        cm = cm - (x_m.reshape(n,1) - zeta_c[i].real)*cpm*dny - zeta_c[i].imag*cpm*dnx
    
    cxp = cxp/(4*c)
    cyp = cyp/(4*c)
    cm = cm/((4*c)**2)
    
    #ピッチングモーメントの描画
    plt.figure(4)
    plt.plot(x_m, cm, color="black", linewidth=1.0)
    plt.xlim([-0.6,-0.3])
    plt.ylim([-0.7,-0.3])


plt.show()
