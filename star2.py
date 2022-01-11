import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D


def f(xa,ya,za,xb,yb,zb):
    R=(xa-xb)**2+(ya-yb)**2+(za-zb)**2
    return (xb-xa)/(np.sqrt(R))**3     #力の定義

N=50     #星の数
E=0 
x_hist=[]
y_hist=[]
z_hist=[]  
R0_hist=[]  
v0_hist=[]    #初期値リスト

while E>=0:
    E=0
    x0_hist=[random.uniform(0,25) for i in range(N)]
    y0_hist=[random.uniform(0,25) for i in range(N)]
    z0_hist=[random.uniform(0,25) for i in range(N)]
    vx_hist=[random.uniform(-5,5) for i in range(N)]
    vy_hist=[random.uniform(-5,5) for i in range(N)]
    vz_hist=[random.uniform(-5,5) for i in range(N)]    #乱数で初期値を決める
    xg=sum(x0_hist)/N
    yg=sum(y0_hist)/N
    zg=sum(z0_hist)/N                                  #重心
    for j,k,l in zip(x0_hist,y0_hist,z0_hist):
        x=j-xg
        y=k-yg
        z=l-zg
        x_hist.append(x)
        y_hist.append(y)
        z_hist.append(z)                               #重心を原点と置いたときの座標
    for j,k,l,n,m,q in zip(x_hist,y_hist,z_hist,vx_hist,vy_hist,vz_hist):
         R00=np.sqrt(j**2+k**2+l**2)
         R0_hist.append(R00)
         v0=np.sqrt(n**2+m**2+q**2)
         v0_hist.append(v0)
         K=0
         P=0
         K += (1/2)*(np.sqrt(n**2+m**2+q**2))**2
         for o,p,r in zip(x_hist,y_hist,z_hist):
            RR=(j-o)**2+(k-p)**2+(l-r)**2
            if RR !=0:
               P += -1/(np.sqrt(RR))
    E = K+P  
    if E>0:
        x0_hist.clear()
        y0_hist.clear()
        z0_hist.clear()
        x_hist.clear()
        y_hist.clear()
        z_hist.clear()
        vx_hist.clear()
        vy_hist.clear()
        vz_hist.clear()                                #E＜０の確認

print(E)

dt=0.01  #時間間隔
T=150.0  #時間上限
t=0.0

x1_hist=[]
y1_hist=[]
z1_hist=[]
vx1_hist=[]
vy1_hist=[]
vz1_hist=[]  #t+dtにおけるリスト

v_hist=[]
R_hist=[]  
E_hist=[]  
Ma_hist=[]    #最終的な速度と（重心からの）距離のリスト

Fx_hist=[]
Fy_hist=[]
Fz_hist=[] 


#それぞれの星の座標と速度を初期値リスト（tリスト）からt+dtにおける各値を計算し、計算後、tリストにt+dtリストに代入し、t+dtリストを消去するループ

while t<T:                                                                      
    for i,j,k,l,n,m in zip(x_hist,y_hist,z_hist,vx_hist,vy_hist,vz_hist):
        Fx=0
        Fy=0
        Fz=0
        for o,p,q in zip(x_hist,y_hist,z_hist):
            R0=(i-o)**2+(j-p)**2+(k-q)**2
            if R0 !=0:
                Fx += f(i,j,k,o,p,q)
                Fy += f(j,i,k,p,o,q)
                Fz += f(k,i,j,q,o,p)                     #すべての星からの力の合計
        Fx_hist.append(Fx)
        Fy_hist.append(Fy)
        Fz_hist.append(Fz)
        x1=i+dt*l+(1/2)*Fx*(dt)**2
        y1=j+dt*n+(1/2)*Fy*(dt)**2
        z1=k+dt*m+(1/2)*Fz*(dt)**2                      #leap-frog座標
        x1_hist.append(x1)
        y1_hist.append(y1)
        z1_hist.append(z1)
    for i,j,k,l,n,m,r,s,u in zip(x1_hist,y1_hist,z1_hist,vx_hist,vy_hist,vz_hist,Fx_hist,Fy_hist,Fz_hist):
        Fx1=0
        Fy1=0
        Fz1=0
        for o,p,q in zip(x1_hist,y1_hist,z1_hist):
            R1=(i-o)**2+(j-p)**2+(k-q)**2
            if R1 !=0:
                Fx1 += f(i,j,k,o,p,q)
                Fy1 += f(j,i,k,p,o,q)
                Fz1 += f(k,i,j,q,o,p)                  #t+dtにおける力の合計
        vx1=l+(1/2)*dt*(r+Fx1)
        vy1=n+(1/2)*dt*(s+Fy1)
        vz1=m+(1/2)*dt*(u+Fz1)                         #leap-frog速度
        vx1_hist.append(vx1)
        vy1_hist.append(vy1) 
        vz1_hist.append(vz1)      
    for o in range(N):
        x_hist[o]=x1_hist[o]
        y_hist[o]=y1_hist[o]
        z_hist[o]=z1_hist[o]
        vx_hist[o]=vx1_hist[o]
        vy_hist[o]=vy1_hist[o]  
        vz_hist[o]=vz1_hist[o]                        #t+dtリストをtリストに代入        
    x1_hist.clear()
    y1_hist.clear()
    z1_hist.clear()
    vx1_hist.clear()
    vy1_hist.clear()  
    vz1_hist.clear()     
    Fx_hist.clear()
    Fy_hist.clear()
    Fz_hist.clear()                                  #t+dtリスト及び力リストを消去
    t=t+dt


for i,j,k,l,n,m in zip(vx_hist,vy_hist,vz_hist,x_hist,y_hist,z_hist):
    v=np.sqrt(i**2+j**2+k**2)
    RT=np.sqrt(l**2+n**2+m**2)
    E=1/2*v**2
    R_hist.append(RT)
    v_hist.append(v) 
    E_hist.append(E)                                 #最終的な（t=Tにおける)距離と速度

kT=sum(E_hist)/N

for i,j in zip(v_hist,E_hist):
    Ma=i**2*np.exp(-j/kT)
    Ma_hist.append(Ma)

plt.subplot(2,1,1)
plt.xlabel('v')
plt.ylabel('v**2exp(-E/kT)')
plt.scatter(v_hist,Ma_hist)


plt.subplot(2,1,2)
plt.hist(v_hist, bins=30)

plt.show()




