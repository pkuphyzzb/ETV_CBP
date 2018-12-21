from numpy import cos,sin,pi,sqrt,arctan
import matplotlib.pyplot as pl 
import numpy as np
def K11m(e,omega):
	return 3/4*e**2+3/16*e**4+3/32*e**6+(e+0.5*e**3+0.25*e**5)*sin(omega)+(51/40*e**2+37/80*e**4+241/640*e**6)*cos(2*omega)-3/16*e**3*sin(3*omega)
	-(1/16*e**4-1/16*e**6)*cos(4*omega)-1/16*e**5*sin(5*omega)+3/64*e**6*cos(6*omega)


def K11p(e,omega):
	return 3/4*e**2+3/16*e**4+3/32*e**6-(e+0.5*e**3+0.25*e**5)*sin(omega)+(51/40*e**2+37/80*e**4+241/640*e**6)*cos(2*omega)+3/16*e**3*sin(3*omega)
	-(1/16*e**4-1/16*e**6)*cos(4*omega)+1/16*e**5*sin(5*omega)+3/64*e**6*cos(6*omega)

def K12m(e,omega):
	return -(e-1/2*e**3-1/4*e**5)*cos(omega)+(51/40*e**2+87/80*e**4+541/640*e**6)*sin(2*omega)-3/16*e**3*cos(3*omega)-(1/16*e**4+5/32*e**6)*sin(4*omega)
	+1/16*e**5*cos(5*omega)+3/64*e**6*sin(6*omega)


def K12p(e,omega):
	return (e-1/2*e**3-1/4*e**5)*cos(omega)+(51/40*e**2+87/80*e**4+541/640*e**6)*sin(2*omega)+3/16*e**3*cos(3*omega)-(1/16*e**4+5/32*e**6)*sin(4*omega)
	-1/16*e**5*cos(5*omega)+3/64*e**6*sin(6*omega)

def K1m(e,omega):
	return -e*sin(omega)+3/4*e**2*cos(2*omega)

def K1p(e,omega):
	return e*sin(omega)+3/4*e**2*cos(2*omega)

def f(e,omega):
	return 1+25/8*e**2+15/8*e**4+95/64*e**6

omega1grid=np.linspace(0,np.pi*2,13)
e1grid=np.linspace(0,0.95,20)

grid=[]
for w in omega1grid:
	line=[]
	for e in e1grid:
		line.append(abs(2*f(e,w)+3*K1m(e,w)+5*sin(2*w)*K12m(e,w))/abs(2*f(e,w)+3*K1p(e,w)+5*sin(2*w)*K12p(e,w)))
	grid.append(line)

grid=np.transpose(grid)
print(grid)
pl.matshow(np.log(grid),vmin=-3,vmax=3,cmap='jet',aspect=0.618)
ax = pl.gca();
ax.set_xticks(np.arange(0, 13, 3));
ax.set_yticks(np.arange(0, 20, 5));
ax.set_xticklabels(['0','$0.5\pi$','$\pi$','$1.5\pi$']);
ax.set_yticklabels(['0','0.25','0.5','0.75']);
pl.xlabel('$\omega_1$')
pl.ylabel('$e_1$')
#pl.title('The Simulated Amplitude Ratio of ETVs Induced by a Polar Planet on $e_1-\omega_1$ mesh grids')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
pl.colorbar()
pl.show()
