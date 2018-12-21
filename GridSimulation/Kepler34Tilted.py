import rebound
from orbital.elements import KeplerianElements
import numpy as np 
import matplotlib.pyplot as pl 
from numpy import sin,cos,tan,arctan,arcsin,arccos
from scipy import optimize
from scipy.special import ellipkinc
from scipy.signal import periodogram
from scipy.signal import lombscargle
from scipy.stats import linregress
from mpl_toolkits.mplot3d import Axes3D

e1=0.52087

def GeneratePlanetaryParameters(i1,Omega1):
	omega2=2*np.random.rand()*np.pi
	Omega=np.pi/2+1.2468029277718442
	#i2=np.random.rayleigh(5)*np.pi/180
	i2=np.pi/2
	#i2=np.arccos(np.random.uniform(cos(5/180*np.pi),1))#*((np.random.rand()>0.5)*2-1)
	#i2=np.random.normal(90,5)*np.pi/180
	print(i2*180/np.pi)
	a=sin(i1+i2)*(cos(Omega)+sin(i1)**2*(1-cos(Omega)))+cos(i1+i2)*sin(i1)*cos(i1)*(1-cos(Omega))
	b=cos(i1)*sin(Omega)*sin(i1+i2)-sin(i1)*sin(Omega)*cos(i1+i2)
	c=cos(i1)*sin(i1)*(1-cos(Omega))*sin(i1+i2)+cos(i1+i2)*(cos(Omega)+cos(i1)*cos(i1)*(1-cos(Omega)))
	Truei=arccos(c)
	b/=sin(Truei)
	a/=sin(Truei)
	if a>0:
		TrueO=arctan(b/a)
	else:
		TrueO=arctan(b/a)+np.pi
	#print(TrueO)
	#print(Truei)
	print('mutualinc: '+str(i2/np.pi*180))
	print('mutualinc2: '+str(np.arccos(np.sin(i1)*np.sin(Truei)*np.cos(TrueO)+np.cos(i1)*np.cos(Truei))/np.pi*180))
	h=cos(i2)**2-0.5*e1**2*sin(i2)**2*(3-5*cos(2*Omega))
	print(h)
	print(e1**2)
	print(sin(i1)*sin(i2)*cos(TrueO-Omega1)+cos(i1)*cos(i2))
	return omega2,Truei,TrueO,h


#Kepler-34
sim=rebound.Simulation()
sim.G=4*np.pi**2
sim.add(m=1.0479)
#e=0.52087
#m=1/1047*0.220
sim.add(m=1.0208,P=27.7958103/365.2422,e=0.52087,omega=1.2468029277718442,M=300.1870/180*np.pi-1.2468029277718442,inc=89.8584/180*np.pi)
#sim.add(m=1/1047*0.220,a=1.0896,e=0.182,omega=0.12800602365751946,M=106.5/180*np.pi-0.12800602365751946,inc=0/180*np.pi,Omega=90/180*np.pi)
o2,i2,O2,h=GeneratePlanetaryParameters(89.8584/180*np.pi,0)
print(o2,i2,O2)
#0.00123006053303
sim.add(m=0.00123006053303,a=1.06484331396,e=0.014,omega=1,M=-1.85637434477,inc=i2,Omega=O2)
#sim.add(m=1/1047*0.220,P=288.822/365.2422,e=0.182,omega=0.13933934075776955,M=106.5/180*np.pi-0.13933934075776955,inc=90.355/180*np.pi,Omega=-1.74/180*np.pi)
R1=1.1618*0.004649130334214449
R2=1.0927*0.004649130334214449
R=R1+R2
#sim.add(m=0.00018854,a=1.0896,e=0.182,omega=2.9938,M=2.1110,inc=1.57408971576,Omega=0.56832494976)


#Kepler-35
'''
sim.add(m=0.8877)
sim.add(m=0.8094,a=0.17617,e=0.1421,omega=1.509931540233465,M=89.1784/180*np.pi-1.509931540233465,inc=90.4238/180*np.pi)
sim.add(m=1/1047*0.127,a=0.60347,e=0.042,omega=1.1186424645091211,M=136.4/180*np.pi-1.1186424645091211,inc=90.76/180*np.pi,Omega=-1.24/180*np.pi)
R1=1.0284*0.004649130334214449
R2=0.7861*0.004649130334214449
R=R1+R2
'''

#Kepler-16
'''
sim.add(m=0.6897)
sim.add(m=0.20255,a=0.22431,e=0.15944,omega=4.5983150407372015,M=92.3520/180*np.pi+4.5983150407372015,inc=90.3401/180*np.pi)
sim.add(m=1/1047*0.333,a=0.7048,e=0.000685,omega=-0.72077013022502856,M=136.4/180*np.pi+0.72077013022502856,inc=90.0322/180*np.pi,Omega=0.003/180*np.pi)
R1=0.6489*0.004649130334214449
R2=0.22623*0.004649130334214449
R=R1+R2
'''


M01=1.0479+1.0208
beta=1.0479*1.0208/(1.0479+1.0208)
a2=1.0896
a1=0.22882
k2=5*e1**2/(1-e1**2)*(1-h)/(h+4*e1**2)
print(k2)
P=8/3/np.pi*M01/beta*(a2/a1)**3.5*(1-0.182**2)**2/((1-e1**2)*(h+4*e1**2))**0.5*26.78/365.2422
print(P)



sim.move_to_com()

#fout=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/ReboundCodes/orbit.txt','w')


time=np.linspace(0,3.98,4000)
Period=(27.79578193-5/26/1440)/365.242199174

dyAB=0
olddyAB=0
dyAb=0
olddyAb=0
dyBb=0
olddyBb=0

ttvAb=[]
ttvBb=[]

ttvAB=[]
ttvBA=[]
transitABlist=[]
transitBAlist=[]
transitAblist=[]
transitBblist=[]


def DyAb(t):
	sim.integrate(t)
	p=sim.particles
	r=np.array([p[0].x-p[2].x,p[0].y-p[2].y])
	v=np.array([p[0].vx-p[2].vx,p[0].vy-p[2].vy])
	return np.dot(r,v)
	return p[0].x-p[2].x

def DyBb(t):
	sim.integrate(t)
	p=sim.particles
	r=np.array([p[1].x-p[2].x,p[1].y-p[2].y])
	v=np.array([p[1].vx-p[2].vx,p[1].vy-p[2].vy])
	return np.dot(r,v)
	return p[1].x-p[2].x

def DyAB(time):
	sim.integrate(time)
	p0=sim.particles
	r0=np.array([p0[0].x-p0[1].x,p0[0].y-p0[1].y])
	v0=np.array([p0[0].vx-p0[1].vx,p0[0].vy-p0[1].vy])
	return np.dot(r0,v0)
	#return p[0].x-p[1].x

#pl.figure()



oldt=time[0]
orbit1=[]
prbit2=[]
tim=[]
a1=[]
e1=[]
inc1=[]
Omega1=[]
omega1=[]
f1=[]
a2=[]
e2=[]
inc2=[]
Omega2=[]
omega2=[]
f2=[]
priduration=[]
secduration=[]
for i in range(len(time)):
	t=time[i]
	sim.integrate(t)
	p=sim.particles
	parameters=sim.calculate_orbits()
	
	#orbit1=sim.calculate_orbits(heliocentric=True)[1]
	if i%10==0:
		print(t)
		#print(parameters[0])
		a1.append(parameters[0].a)
		e1.append(parameters[0].e)
		inc1.append(parameters[0].inc*180/np.pi)
		Omega1.append(parameters[0].Omega*180/np.pi)
		omega1.append(parameters[0].omega*180/np.pi)
		f1.append(parameters[0].f*180/np.pi)
		a2.append(parameters[1].a)
		e2.append(parameters[1].e)
		inc2.append(parameters[1].inc*180/np.pi)
		Omega2.append(parameters[1].Omega*180/np.pi)
		omega2.append(parameters[1].omega*180/np.pi)
		f2.append(parameters[1].f*180/np.pi)
		#orbit2.append(arameters[1])
		tim.append(t)
		pl.scatter(-p[0].x,-p[0].z,s=1,c='k')
		pl.scatter(-p[1].x,-p[1].z,s=1,c='k')
		pl.scatter(-p[2].x,-p[2].z,s=1,c='r')
		position=np.array([p[2].x,p[2].y,p[2].z])
		momentum=np.array([p[2].vx,p[2].vy,p[2].vz])
		angularmomentum=np.cross(position,momentum)#/np.sqrt(np.dot(np.cross(position,momentum),np.cross(position,momentum)))
		position2=np.array([p[2].x,p[2].y,p[2].z])
		momentum2=np.array([p[2].vx,p[2].vy,p[2].vz])
		#fout.write(str(angularmomentum[0])+' '+str(angularmomentum[1])+' '+str(angularmomentum[2])+'\n')
		#fout.write(str(p[0].x)+' '+str(p[0].y)+' '+str(p[0].z)+' '+str(p[1].x)+' '+str(p[1].y)+' '+str(p[1].z)+' '+str(p[2].x)+' '+str(p[2].y)+' '+str(p[2].z)+'\n')
	olddyAB=DyAB(oldt)
	olddyAb=DyAb(oldt)
	olddyBb=DyBb(oldt)
	dyAB=DyAB(t)#p[0].x-p[1].x
	dyAb=DyAb(t)#p[0].x-p[2].x
	dyBb=DyBb(t)#p[1].x-p[2].x	
	if dyAb>0 and olddyAb<0 and p[2].z>0 and abs(p[2].y)<0.1 :
		#t0=oldt-(t-oldt)/(olddyAb-dyAb)*olddyAb
		transitAb=optimize.bisect(DyAb,oldt,t)
		sim.integrate(transitAb)
		p=sim.particles
		#if np.sqrt((p[0].x-p[2].x)*(p[0].x-p[2].x)+(p[0].y-p[2].y)*(p[0].y-p[2].y))<1/200:
		#pl.scatter(-p[0].x,-p[0].z,s=50,c='g')
		#pl.scatter(-p[2].x,-p[2].z,s=50,c='r')
		transitAblist.append(transitAb-30.8/365.242199174)
	if dyBb>0 and olddyBb<0 and p[2].z>0 and abs(p[2].y)<0.1:
		#t0=oldt-(t-oldt)/(olddyBb-dyBb)*olddyBb
		transitBb=optimize.bisect(DyBb,oldt,t)
		sim.integrate(transitBb)
		p=sim.particles
		#if np.sqrt((p[1].x-p[2].x)*(p[1].x-p[2].x)+(p[1].y-p[2].y)*(p[1].y-p[2].y))<1/200:
		#pl.scatter(-p[1].x,-p[1].z,s=50,c='g')
		#pl.scatter(-p[2].x,-p[2].z,s=50,c='r')
		transitBblist.append(transitBb-30.8/365.242199174)
	if dyAB>0 and olddyAB<0 and p[0].z<0 and np.abs(p[0].x)<0.05:
		#t0=oldt-(t-oldt)/(olddyAB-dyAB)*olddyAB
		try:
			transitAB=optimize.bisect(DyAB,oldt,t)
		except:
			print('oops!')
			pass
		else:
			transitABlist.append(transitAB)
			sim.integrate(transitAB)
			chord=np.sqrt((p[0].x-p[1].x)**2+(p[0].y-p[1].y)**2)
			duration=2*np.sqrt(R**2-chord**2)/np.abs(np.sqrt((p[0].vx-p[1].vx)**2+(p[0].vy-p[1].vy)**2))*365.2422*24
			priduration.append(duration)
		#sim.status()
			#p=sim.particles
		#print(p[0].x-p[1].x)
		#print(p[0].x)
		#print(p[1].x)
			#if 1.4<transitAB<1.5:
			#pl.scatter(-p[0].x,-p[0].y,s=50,c='r')
			#pl.scatter(-p[1].x,-p[1].y,s=50,c='m')
			#pl.scatter(-p[2].x,-p[2].y,s=50,c='b')
			#print(transitAB*365.242199174-31)
		#transitAblist.append(p[0].y-p[2].y)
		# B occults A
	if dyAB>0 and olddyAB<0 and p[0].z>0 and np.abs(p[0].x)<0.1:
		#t0=oldt-(t-oldt)/(olddyAB-dyAB)*olddyAB
		try:
			transitBA=optimize.bisect(DyAB,oldt,t)
		except:
			print('oops!')
			pass
		else:
			transitBAlist.append(transitBA)
			sim.integrate(transitBA)
			chord=np.sqrt((p[0].x-p[1].x)**2+(p[0].y-p[1].y)**2)
			duration=2*np.sqrt(R**2-chord**2)/np.abs(np.sqrt((p[0].vx-p[1].vx)**2+(p[0].vy-p[1].vy)**2))*365.2422*24
			secduration.append(duration)
		#print(transitBA)
		#sim.status()
			#p=sim.particles
		#pl.scatter(p[0].x,p[0].z,s=50,c='b')
		#pl.scatter(p[1].x,p[1].z,s=50,c='k')
		# A occults B
	oldt=t
pl.xlim([-1.5,1.5])
pl.ylim([-1.5,1.5])
pl.axis('equal')
pl.savefig('orbit')
pl.show()	

print(priduration)
print(secduration)
pl.subplot(2,1,1)
pl.xlim([0,4])
pl.title('Eclipsing Duration Variations for the Best Fitting Polar Model of Kepler-34 System (Simulation)')
pl.scatter(transitABlist,np.array(priduration)*3600,c='k',s=5)
#pl.ylim([5.55625*3600,5.5566*3600])
pl.ylabel('Primary Duration [s]')
pl.grid('on')
pl.subplot(2,1,2)
#pl.ylim([16.352*3600,16.3535*3600])
pl.grid('on')
pl.scatter(transitBAlist,np.array(secduration)*3600,c='k',s=5)
pl.ylabel('Secondary Duartion [s]')
pl.xlabel('time [y]')
pl.xlim([0,4])
pl.show()

print(str(omega1[-1]-omega1[0]))

#print(transitABlist)
pl.figure()
pl.subplot(6,2,1)
ax=pl.subplot(6,2,1)
ax.set_yticks(np.linspace(0.228819,0.228823,3))
ax.set_yticklabels(['0.228819','0.228821','0.228823'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.title('Orbital Parameter Variations in the Polar Model')
pl.plot(tim,a1)
pl.ylabel('a1 [AU]')
pl.grid('on')

pl.subplot(6,2,2)
ax=pl.subplot(6,2,2)
pl.grid('on')
ax.set_yticks(np.linspace(1.052,1.068,5))
ax.set_yticklabels(['1.052','1.056','1.060','1.064','1.068'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.plot(tim,a2)
pl.ylabel('a2[AU]')
pl.subplot(6,2,3)
pl.plot(tim,e1)
ax=pl.subplot(6,2,3)
pl.grid('on')
ax.set_yticks(np.linspace(0.52085,0.52088,4))
ax.set_yticklabels(['0.52085','0.52086','0.52087','0.52088'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.ylabel('e1')
pl.subplot(6,2,4)
ax=pl.subplot(6,2,4)
pl.grid('on')
ax.set_yticks(np.linspace(0,0.03,4))
ax.set_yticklabels(['0','0.01','0.02','0.02'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.plot(tim,e2)
pl.ylabel('e2')
pl.subplot(6,2,5)
pl.plot(tim,inc1)
pl.ylabel('inc1 [deg]')
ax=pl.subplot(6,2,5)
pl.grid('on')
ax.set_yticks(np.linspace(89.857,89.8585,2))
ax.set_yticklabels(['89.857','89.8575','89.8585'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.subplot(6,2,6)
pl.plot(tim,inc2)
pl.ylabel('inc2 [deg]')
ax=pl.subplot(6,2,6)
pl.grid('on')
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.subplot(6,2,7)
ax=pl.subplot(6,2,7)
pl.grid('on')
ax.set_yticks(np.linspace(0,0.001,3))
ax.set_yticklabels(['0','0.0005','0.001'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.plot(tim,Omega1)
pl.ylabel('Omega1 [deg]')
pl.subplot(6,2,8)
pl.plot(tim,Omega2)
pl.ylabel('Omega2 [deg]')
ax=pl.subplot(6,2,8)
pl.grid('on')
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.subplot(6,2,9)
pl.plot(tim,omega1)
ax=pl.subplot(6,2,9)
pl.grid('on')
ax.set_yticks(np.linspace(71.2,71.4,3))
ax.set_yticklabels(['71.2','71.3','71.4'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.ylabel('omega1 [deg]')
pl.subplot(6,2,10)
pl.plot(tim,omega2)
ax=pl.subplot(6,2,10)
pl.grid('on')
ax.set_yticks(np.linspace(-200,200,5))
ax.set_yticklabels(['-200','-100','0','100','200'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels([' ',' ',' ',' ',' '])
pl.ylabel('omega2 [deg]')
pl.subplot(6,2,11)
pl.plot(tim,f1)
ax=pl.subplot(6,2,11)
pl.grid('on')
ax.set_yticks(np.linspace(-300,100,5))
ax.set_yticklabels(['-300','-200','-100','0','100'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels(['0 ','1 ','2 ','3 ','4 '])
pl.ylabel('f1 [deg]')
pl.xlabel('Time [y]')
pl.subplot(6,2,12)
pl.plot(tim,f2)
pl.ylabel('f2 [deg]')
ax=pl.subplot(6,2,12)
pl.grid('on')
ax.set_yticks(np.linspace(-300,200,6))
ax.set_yticklabels(['-300','-200','-100','0','100','200'])
ax.set_xticks(np.linspace(0,4,5))
ax.set_xticklabels(['0 ','1 ','2 ','3 ','4 '])
pl.xlabel('Time [y]')
#pl.savefig('orbitalelements.pdf')
pl.show()

pl.figure()
for t in np.linspace(0,2,10000):
	sim.integrate(t)
	p=sim.particles
	if transitABlist[0]<t<transitABlist[1]:
		pl.scatter((t-transitABlist[0])*365.242199174/27.79578193,p[1].vz*4.753314680261785,c='r')
		pl.scatter((t-transitABlist[0])*365.242199174/27.79578193,p[0].vz*4.753314680261785,c='b')
pl.savefig('RVdiagram.pdf')
pl.show()

#Period=27.79534096605/365.242199174
Period=27.795392620850002/365.242199174
print(transitABlist)

ttvAB=np.array([(transitABlist[i]-transitABlist[0]-i*Period)*24*60*365.242199174 for i in range(len(transitABlist))])
ttvBA=np.array([(transitBAlist[i]-transitBAlist[0]-i*Period)*24*60*365.242199174 for i in range(len(transitBAlist))])

transitABlist=np.array([transitABlist[i]*365.242199174 for i in range(len(transitABlist))])
transitBAlist=np.array([transitBAlist[i]*365.242199174 for i in range(len(transitBAlist))])
transitAblist=np.array([transitAblist[i]*365.242199174 for i in range(len(transitAblist))])
transitBblist=np.array([transitBblist[i]*365.242199174 for i in range(len(transitBblist))])


truettvAB=ttvAB-linregress(transitABlist,ttvAB)[0]*transitABlist-linregress(transitABlist,ttvAB)[1]
truettvBA=ttvBA-linregress(transitBAlist,ttvBA)[0]*transitBAlist-linregress(transitBAlist,ttvBA)[1]
#print(linregress(transitBAlist,ttvBA)[0])
#truettvAB=ttvAB-27.79534096605*transitABlist-linregress(transitABlist,ttvAB)[1]
#truettvBA=ttvBA-27.79534096605*transitBAlist-linregress(transitBAlist,ttvBA)[1]

pl.figure()
#pl.plot(transitABlist-30.8,ttvAB-np.min(ttvAB),c='b')
#pl.scatter(transitABlist-30.8,ttvAB-np.min(ttvAB),c='b')
pl.plot(transitABlist,truettvAB,c='b')
pl.scatter(transitABlist,truettvAB,c='b',s=1)
#pl.xlabel('time [d]')
#pl.ylabel('TTV [min]')
#pl.title('Transit Timing Variations of the Secondary Transits in Kepler-34 System')
#pl.savefig('SecondaryTTV.pdf')
#pl.show()

#pl.figure()
#pl.plot(transitBAlist-30.8,ttvBA-np.max(ttvBA),c='r')
#pl.scatter(transitBAlist-30.8,ttvBA-np.max(ttvBA),c='r')
pl.plot(transitBAlist,truettvBA,c='r')
pl.scatter(transitBAlist,truettvBA,c='r',s=1)
pl.xlabel('time [d]')
pl.ylabel('TTV [min]')
pl.title('Transit Timing Variations of the Primary Transits in Kepler-34 System')
pl.savefig('TTV.pdf')
pl.show()

#pl.plot(truettvAB)
#pl.show()
pl.figure()
pl.subplot(2,1,1)
f,Pxx_den=periodogram(truettvAB,1/27.79578193)
print(f)
pl.loglog(1/np.array(f),Pxx_den)
pl.grid('on')
pl.subplot(2,1,2)
f,Pxx_den=periodogram(truettvBA,1/27.79578193)
pl.loglog(1/np.array(f),Pxx_den)
pl.grid('on')
pl.show()

fout2=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/ReboundCodes/TTVs.txt','w')
#fout2=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/ReboundCodes/TTVs.txt','w')
#fout2.write(transitABlist)
print(transitABlist+2454969.2-2454900)
#print(truettvAB)
print(transitBAlist+2454969.2-2454900)
#print(truettvBA)
print(transitAblist)
print(transitBblist)
for i in range(len(transitABlist)):
	fout2.write(str(transitABlist[i]+2454969.2-2454900)+' '+str(truettvAB[i])+' '+str(transitBAlist[i]+2454969.2-2454900)+' '+str(truettvBA[i])+'\n')




