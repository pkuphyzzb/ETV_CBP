from numpy import cos,sin,pi,sqrt,arctan
import matplotlib.pyplot as pl 
import numpy as np

omega1grid=np.linspace(0,np.pi*2,13)
e1grid=np.linspace(0,0.1,21)

grid=[]
P1=27.79
P2=300
def Ratioc(e,w,P1,P2):
	A=sqrt(e**2+121/900*P1**2/P2**2+11/15*e*sin(w)*P1/P2)
	B=sqrt(e**2+121/900*P1**2/P2**2-11/15*e*sin(w)*P1/P2)
	return A/B
for w in omega1grid:
	line=[]
	for e in e1grid:
		line.append(Ratioc(e,w,P1,P2))
	grid.append(line)

grid=np.transpose(grid)
print(grid)
pl.matshow(np.log(grid),vmin=-3,vmax=3,cmap='jet',aspect=0.618)
ax = pl.gca();
ax.set_xticks(np.arange(0, 13, 3));
ax.set_yticks(np.arange(0, 20, 5));
ax.set_xticklabels(['0','$0.5\pi$','$\pi$','$1.5\pi$']);
ax.set_yticklabels(['0','0.025','0.05','0.075']);
pl.xlabel('$\omega_1$')
pl.ylabel('$e_1$')
#pl.title('The Simulated Amplitude Ratio of ETVs Induced by a Coplanar Planet on $e_1-\omega_1$ Mesh Grids At Small $e_1$ Region')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)
pl.colorbar()
pl.show()