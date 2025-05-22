import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import math
ee=0.01
E=ee
gamma=4/3
n=3
rc=0.75/E
acs=1/(2*rc)
ac=math.sqrt(acs)
delr=0.1
i1=ac+(4*(ac**3)/(2*n+1))*(1+math.sqrt(n*(n-1.5)))*delr
i2=ac+(4*(ac**3)/(2*n+1))*(1-math.sqrt(n*(n-1.5)))*delr
i3=1.5*ac
i4=0.5*ac
r1=np.arange(rc+delr,1000,0.1)
r2=np.arange(rc-delr,1.0,-0.1)

def bondi(u,r):
    dudr=(1/(r**2)-(2/(r*n))*(E-(u**2)/2+1/r))/((E-(u**2)/2+1/r)/(u*n)-u)
    return dudr

u1=odeint(bondi,i1,r1)
u2=odeint(bondi,i2,r1)
u3=odeint(bondi,i3,r1)
u4=odeint(bondi,i4,r1)
u5=odeint(bondi,i1,r2)
u6=odeint(bondi,i2,r2)
u7=odeint(bondi,i3,r2)
u8=odeint(bondi,i4,r2)

m1=[(u/math.sqrt((E+1/r1[i]-(u**2)/2)*(1/n)))[0] for i,u in enumerate(u1)]
m2=[(u/math.sqrt((E+1/r1[i]-(u**2)/2)*(1/n)))[0] for i,u in enumerate(u2)]
m3=[(u/math.sqrt((E+1/r1[i]-(u**2)/2)*(1/n)))[0] for i,u in enumerate(u3)]
m4=[(u/math.sqrt((E+1/r1[i]-(u**2)/2)*(1/n)))[0] for i,u in enumerate(u4)]
m5=[(u/math.sqrt((E+1/r2[i]-(u**2)/2)*(1/n)))[0] for i,u in enumerate(u5)]
m6=[(u/math.sqrt((E+1/r2[i]-(u**2)/2)*(1/n)))[0] for i,u in enumerate(u6)]
m7=[(u/math.sqrt((E+1/r2[i]-(u**2)/2)*(1/n)))[0] for i,u in enumerate(u7)]
m8=[(u/math.sqrt((E+1/r2[i]-(u**2)/2)*(1/n)))[0] for i,u in enumerate(u8)]

plt.plot(r1, m1, 'm-', label="Sph Winds", linewidth=2.5)
plt.plot(r1, m2, 'r', label="Sph Accretion", linewidth=2.5)
plt.plot(r1, m3, 'b--')
plt.plot(r1, m4, 'b--')
plt.plot(r2, m5, 'r', linewidth=2.5)
plt.plot(r2, m6, 'm-', linewidth=2.5)
plt.plot(r2, m7, 'b--')
plt.plot(r2, m8, 'b--')
ax = plt.gca()
ax.set_xscale('log')
plt.xlabel('scaled r')
plt.ylabel('Mach number')
plt.title('E=%1.3f' % E)
plt.legend()
plt.text(rc+127, 0.93, 'Critical point')
plt.text(100, 2.60, 'Supersonic')
plt.text(45, -0.06, 'Subsonic')
plt.ylim(-0.31, 6.5)
plt.xlim(0, 1283)
plt.show()