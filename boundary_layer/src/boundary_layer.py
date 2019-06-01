# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

rho0 = 1.225    # 空気密度
T0 = 293.0      # 絶対零度
visc = 1.458*10**(-6)*T0**1.5/(T0+110.4) # 粘度
width = 0.5     # 幅
v0 = 20.0       # 速度
re0 = rho0 * v0 * width/visc    # レイノルズ数
blt = width*5.3/np.sqrt(re0)

def init_print():
    print("Solve Boundary Layer Equation for Flat Plate")
    print("Flow Conditions")
    print("Velocity (m/s) =",v0)
    print("Visc (m^2/s)",visc)
    print("Reynolds Number =",re0)
    print("Length (m) =",width)
    print("Boundary Layer Thickness (mm) =",blt*1000.0)

nx = 60000
dx = width/float(nx)
ny = 300
ymax = blt*2.0
dy = ymax/float(ny)
y,x = np.meshgrid(np.arange(0,ymax,dy),np.arange(0,width,dx))
u,v = np.zeros([nx,ny]),np.zeros([nx,ny])

def get_u(u_last,v_last,DpDx):
    # 差分方程式
    u_now = np.zeros(ny)
    for i in range(1,ny-1):
        DuDy = (u_last[i+1]-u_last[i-1])/(2.0*dy)
        DuDyy = (u_last[i+1]-2.0*u_last[i]+u_last[i-1])/(dy**2)
        DuDx = (DuDyy*visc/rho0 - v_last[i]*DuDy - DpDx/rho0)/u_last[i]
        u_now[i] = u_last[i] + dx * DuDx
    # 境界条件
    u_now[0] = 0.0
    u_now[ny-1] = v0
    return u_now

def get_v(u_last,v_last,u_now):
    # 差分方程式
    v_now = np.zeros(ny)
    for i in range(1,ny-1):
        DuDx = (u_now[i-1] - u_last[i-1])/dx
        DuDxx = (u_now[i] - u_last[i])/dx
        v_now[i] = v_now[i-1] - (DuDx + DuDxx)*dy/2.0
    return v_now

def get_u_v(DpDx):
    for i in range(ny): #x=0での境界条件
        u[0][i] = v0
    #v[0][i] = 0.0
    for i in range(nx): #y=farでの境界条件
        u[i][ny-1] = v0
    for p in range(1,nx):       # xに関して回す
        u[p] = get_u(u[p-1],v[p-1],DpDx)     # i+1でのuを計算
        v[p][0] = 0.0
        v[p] = get_v(u[p-1],v[p-1],u[p])
    return np.array([u,v])

def get_fig(u_vector,n,graph_number):
    plt.figure(graph_number)
    y = np.arange(0,ymax,dy)
    plt.plot(u_vector + n*30,y*1000)
    return 0

def get_graphs(u,graph_number,kankaku):
    number = int(nx/kankaku)
    for i in range(number-1):
        get_fig(u[(i+1)*kankaku],i,graph_number)
    plt.xlabel("u[m/s]")
    plt.ylabel("y[mm]")
    return 0

def get_boundary_layer_thickness(u):
    delta = np.zeros(nx)
    for i in range(nx):
        j = 0
        while(u[i][j] < u[i][ny-1]*0.995):
            j += 1
        delta[i] = j*dy
    return delta


def get_displacement_thickness(u):
    delta = np.zeros(nx)
    for i in range(nx):
        for j in range(ny):
            delta[i] += (1 - u[i][j]/v0)*dy
    return delta

def get_momentum_thickness(u):
    delta = np.zeros(nx)
    for i in range(nx):
        for j in range(ny):
            delta[i] += u[i][j]/v0 *(1-(u[i][j]/v0)**2)*dy
    return delta

def get_energy_thickness(u):
    delta = np.zeros(nx)
    for i in range(nx):
        for j in range(ny):
            delta[i] += u[i][j]/v0 *(1-u[i][j]/v0)*dy
    return delta

def get_thicknesses_graph(u,i):
    delta_d = get_displacement_thickness(u)
    delta_m = get_momentum_thickness(u)
    delta_e = get_energy_thickness(u)
    delta_b = get_boundary_layer_thickness(u)
    xx = np.arange(0,width,dx)
    plt.figure(i)
    plt.plot(xx,delta_d*10**3)
    plt.plot(xx,delta_m*10**3)
    plt.plot(xx,delta_e*10**3)
    plt.plot(xx,delta_b*10**3)
    plt.legend(("displacement thickness","momentum thickness","energy thickness","boundary layer thickness"))
    plt.ylabel("thickness [mm]")
    plt.xlabel("x[m]")

def main():
    init_print()
    vel_DpDx_60 = get_u_v(40)
    
    vel_DpDx_30 = get_u_v(30)
    vel_DpDx_10 = get_u_v(50)
    #get_graphs(vel_DpDx_zero[0],1,5000)
    #plt.tick_params(labelbottom='off')
    #plt.title("Velocity distribution of boundary layer")
    # get_graphs(vel_DpDx_zero[0],2,13000)
    #get_graphs(vel_DpDx_pos[0],2,13000)
    #get_graphs(vel_DpDx_neg[0],2,13000)
    #plt.ylim(0,3)
    #plt.title("Velocity distribution of boundary layer with pressure gradient")
    get_thicknesses_graph(vel_DpDx_10[0],3)
    plt.title("dp/dx = 40")
    get_thicknesses_graph(vel_DpDx_30[0],4)
    plt.title("dp/dx = 30")
    get_thicknesses_graph(vel_DpDx_60[0],5)
    plt.title("dp/dx = 50")
    print("計算が終了しました。")
    plt.show()
    return 0

if __name__ == '__main__':
    main()
