import matplotlib.pyplot as plt
import numpy as np
import random
from mpl_toolkits.mplot3d import Axes3D


fig=plt.figure()
ax=fig.add_subplot(111, projection='3d')

dt=0.01      
T=50     
t=0
N=10     
E=0

def f(xa,ya,za,xb,yb,zb):
    RR=(xa-xb)**2+(ya-yb)**2+(za-zb)**2
    return (xb-xa)/(np.sqrt(RR))**3       

while E>=0:
    x0=np.array([random.uniform(0,10) for i in range(N)])
    y0=np.array([random.uniform(0,10) for i in range(N)])
    z0=np.array([random.uniform(0,10) for i in range(N)])
    vx=np.array([random.uniform(-2,2) for i in range(N)])
    vy=np.array([random.uniform(-2,2) for i in range(N)])
    vz=np.array([random.uniform(-2,2) for i in range(N)])    
    v0=np.sqrt(vx**2+vy**2+vz**2)
    xg=x0.sum()/N
    yg=y0.sum()/N
    zg=z0.sum()/N
    x=x0-xg
    y=y0-yg
    z=z0-zg                                                 
    k=1/2*(vx**2+vy**2+vz**2)
    K=k.sum()
    u=0
    for i in range(N):
        xd=np.delete(x,i)
        yd=np.delete(y,i)
        zd=np.delete(z,i)
        for j in range(N-1):
            R0=(xd[j]-x[i])**2+(y[j]-y[i])**2+(z[j]-z[i])**2
            u+= -1/np.sqrt(R0)
    E=K+u
print(E)

xgu=x.reshape(N,1)
ygu=y.reshape(N,1)
zgu=z.reshape(N,1)   


while t<T:
    Fx=np.array([])
    Fy=np.array([])
    Fz=np.array([])
    for i in range(N):
        fx=0
        fy=0
        fz=0
        xd=np.delete(x,i)
        yd=np.delete(y,i)
        zd=np.delete(z,i)
        for j in range(N-1):
            fx += f(x[i],y[i],z[i],xd[j],yd[j],zd[j])
            fy += f(y[i],x[i],z[i],yd[j],xd[j],zd[j])
            fz += f(z[i],x[i],y[i],zd[j],xd[j],yd[j])
        Fx=np.append(Fx,[fx])
        Fy=np.append(Fy,[fy])
        Fz=np.append(Fz,[fz])
    x=x+dt*vx+(1/2)*Fx*dt**2
    y=y+dt*vy+(1/2)*Fy*dt**2
    z=z+dt*vz+(1/2)*Fz*dt**2
    xgu=np.append(xgu, x.reshape(N,1), axis=1)
    ygu=np.append(ygu, y.reshape(N,1), axis=1)
    zgu=np.append(zgu, z.reshape(N,1), axis=1)
    Fx1=np.array([])
    Fy1=np.array([])
    Fz1=np.array([])
    for i in range(N):
        fx1=0
        fy1=0
        fz1=0
        xd=np.delete(x,i)
        yd=np.delete(y,i)
        zd=np.delete(z,i)
        for j in range(N-1):
            fx1 += f(x[i],y[i],z[i],xd[j],yd[j],zd[j])
            fy1 += f(y[i],x[i],z[i],yd[j],xd[j],zd[j])
            fz1 += f(z[i],x[i],y[i],zd[j],xd[j],yd[j])
        Fx1=np.append(Fx1,[fx1])
        Fy1=np.append(Fy1,[fy1])
        Fz1=np.append(Fz1,[fz1])
    vx=vx+(1/2)*(Fx+Fx1)*dt
    vy=vy+(1/2)*(Fy+Fy1)*dt
    vz=vz+(1/2)*(Fz+Fz1)*dt
    t=t+dt

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")


for i in range(N):
    ax.plot(xgu[i],ygu[i],zgu[i])

plt.show()
