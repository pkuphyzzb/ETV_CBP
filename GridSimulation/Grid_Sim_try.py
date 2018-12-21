import rebound
import numpy as np 
import matplotlib.pyplot as pl 
from numpy import sin,cos,tan,arctan,arcsin,arccos,pi,sqrt
from scipy import optimize
from scipy.stats import linregress
from scipy import interpolate

def GeneratePlanetaryParameters(i1,Omega1,Omega):
	omega2=2*np.random.rand()*np.pi
	Omega=np.pi/2+Omega
	i2=np.pi/2
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
	return omega2,Truei,TrueO

def DyAB(time):
	sim.integrate(time)
	p0=sim.particles
	r0=np.array([p0[0].x-p0[1].x,p0[0].y-p0[1].y])
	v0=np.array([p0[0].vx-p0[1].vx,p0[0].vy-p0[1].vy])
	return np.dot(r0,v0)

fout1=open('Griddat_pol.dat','w')
fout2=open('Griddat_cop.dat','w')
fout1.write('M2\te1\tP2\tomega1\ts1p\ts2p\ts3p\tphp1\tphp2\tphp3\ts1s\ts2s\ts3s\tphp1\tphp2\tphp3\n')#+str(M2)+' e1='+str(e1)+' P2='+str(P2)+' omega1='+str(omega1))
fout2.write('M2\te1\tP2\tomega1\ts1p\ts2p\ts3p\tphp1\tphp2\tphp3\ts1s\ts2s\ts3s\tphp1\tphp2\tphp3\n')
fmass_radius=open('/Users/apple/Desktop/UnderGraduateResearch/LAMOST-Kepler/codes/Simulation/mass_logg.txt')
mass=[]
radius=[]
for line in fmass_radius.readlines():
	string=line.split(' ')
	m=float(string[0])
	logg=float(string[4])
	G=6.67e-8
	r=np.sqrt(G*m*1.989e33/np.power(10,logg))/6.95e10
	mass.append(m)
	radius.append(r)

def MassToRadius(Mass):
	f=interpolate.interp1d(mass,radius,'cubic')
	if Mass<mass[0] or Mass>np.max(mass):
		return -1
	return f(Mass)

def harmonic(var,T0,f,dT0,P,A1,B1,A2,B2,A3,B3):
	t,n=var
	return T0+P*n+dT0+A1*cos(2*pi*f*t)+B1*sin(2*pi*f*t)+A2*cos(4*pi*f*t)+B2*sin(4*pi*f*t)+A3*cos(6*pi*f*t)+B3*sin(6*pi*f*t)#+A4*cos(8*pi*f*t)+B4*sin(8*pi*f*t)

P1=27.79
M1=1
P2grid=np.linspace(100,2000,96)
M2grid=np.linspace(0.1,1,19)
omega1grid=np.linspace(0,np.pi*2,13)
e1grid=np.linspace(0,0.5,11)


for P2 in P2grid:
	for M2 in M2grid:
		for omega1 in omega1grid:
			for e1 in [0,0.5]:#e1grid:
				P2=300
				e1=0
				omega1=np.pi/2
				M2=0.5
				sim=rebound.Simulation()
				sim.G=4*np.pi**2
				sim.add(m=M1)
				sim.add(m=M2,P=P1/365.2422,e=e1,omega=omega1,M=np.random.rand()*np.pi*2,inc=np.pi/2)
				o2,i2,O2=GeneratePlanetaryParameters(np.pi/2,0,omega1)
				sim.add(m=1/1047*0.220,P=P2/365.2422,e=0,omega=np.random.rand()*np.pi*2,M=np.random.rand()*np.pi*2,inc=i2,Omega=O2)
				R1=0.004649130334214449
				R2=MassToRadius(M2)*0.004649130334214449
				R=R1+R2
				sim.move_to_com()
				parameters=sim.calculate_orbits()
				a1=parameters[0].a
				a2=parameters[1].a
				mu=M2/(M1+M2)
				criticala=(1.6+4.12*mu+5.1*e1-4.27*mu*e1-2.22*e1*e1-5.09*mu*mu+4.61*e1**2*mu**2)*a1
				if criticala>a2:
					print("unstable orbits! "+str(a2)+" lower stable limit: "+str(criticala))
					fout1.write(str(M2)+'\t'+str(e1)+'\t'+str(P2)+'\t'+str(omega1)\
					+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)\
					+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\n')
					fout2.write(str(M2)+'\t'+str(e1)+'\t'+str(P2)+'\t'+str(omega1)\
					+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)\
					+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\t'+str(0)+'\n')
					continue
				Period=P1/365.242199174
				time=np.linspace(0,P2/365.2422*5,4000)
				dyAB=0
				olddyAB=0
				ttvAB=[]
				ttvBA=[]
				transitABlist=[]
				transitBAlist=[]
				oldt=time[0]
				for i in range(len(time)):
					t=time[i]
					sim.integrate(t)
					p=sim.particles
					olddyAB=DyAB(oldt)
					dyAB=DyAB(t)
					if dyAB>0 and olddyAB<0 and p[0].z<0 and np.abs(p[0].x)<0.05:
						transitAB=optimize.bisect(DyAB,oldt,t)
						transitABlist.append(transitAB)
						sim.integrate(transitAB)
					if dyAB>0 and olddyAB<0 and p[0].z>0 and np.abs(p[0].x)<0.1:
						transitBA=optimize.bisect(DyAB,oldt,t)
						transitBAlist.append(transitBA)
						sim.integrate(transitBA)
					oldt=t

				transitABlist=np.array(transitABlist)
				transitBAlist=np.array(transitBAlist)
				
				timing=transitABlist*365.2422
				index=range(len(transitABlist))
				try:
					r=optimize.curve_fit(harmonic,[timing,index],timing,[0,1/P2,0,P1,0,0,0,0,0,0],method='dogbox')
				except:
					print("Oops, cannot fit parameters from standard ETV model")
					#continue
					strength1p=0
					strength2p=0
					strength3p=0
					phase1p=0
					phase2p=0
					phase3p=0
				else:
					#for i in range(len(r[0])):
					#	print (str(r[0][i])+'+-'+str(np.sqrt(r[1][i][i])))
					T0,f,dT0,P,A1,B1,A2,B2,A3,B3=r[0]
					t=np.linspace(0,2,20000)
					model=A1*cos(2*pi*t)+B1*sin(2*pi*t)+A2*cos(4*pi*t)+B2*sin(4*pi*t)+A3*cos(6*pi*t)+B3*sin(6*pi*t)
					#pl.scatter((timing-np.floor(timing*0.5*f)/0.5/f)*f,np.array(timing)-np.array(index)*P-T0-dT0,color='g')
					#pl.ylim([-0.0002,0.0002])
					#pl.plot(t,model)
					#pl.grid('on')
					#pl.show()
					strength1p=sqrt(A1**2+B1**2)
					strength2p=sqrt(A2**2+B2**2)
					strength3p=sqrt(A3**2+B3**2)
					phase1p=arctan(A1/B1)+pi*np.sign(B1)
					phase2p=arctan(A2/B2)+pi*np.sign(B2)
					phase3p=arctan(A3/B3)+pi*np.sign(B3)

				timing=transitBAlist*365.2422
				index=range(len(transitBAlist))
				try:
					r=optimize.curve_fit(harmonic,[timing,index],timing,[0,1/P2,0,P1,0,0,0,0,0,0],method='dogbox')
				except:
					print("Oops, cannot fit parameters from standard ETV model")
					#continue
					strength1s=0
					strength2s=0
					strength3s=0
					phase1s=0
					phase2s=0
					phase3s=0
				else:
					#for i in range(len(r[0])):
					#	print (str(r[0][i])+'+-'+str(np.sqrt(r[1][i][i])))
					T0,f,dT0,P,A1,B1,A2,B2,A3,B3=r[0]
					strength1s=sqrt(A1**2+B1**2)
					strength2s=sqrt(A2**2+B2**2)
					strength3s=sqrt(A3**2+B3**2)
					phase1s=arctan(A1/B1)+pi*np.sign(B1)
					phase2s=arctan(A2/B2)+pi*np.sign(B2)
					phase3s=arctan(A3/B3)+pi*np.sign(B3)
					t=np.linspace(0,2,20000)
					model=A1*cos(2*pi*t)+B1*sin(2*pi*t)+A2*cos(4*pi*t)+B2*sin(4*pi*t)+A3*cos(6*pi*t)+B3*sin(6*pi*t)
					pl.scatter((timing-np.floor(timing*0.5*f)/0.5/f)*f,np.array(timing)-np.array(index)*P-T0-dT0,color='g')
					pl.ylim([-0.0002,0.0002])
					pl.plot(t,model)
					pl.grid('on')
					pl.show()
				
				print('M2='+str(M2)+' e1='+str(e1)+' P2='+str(P2)+' omega1='+str(omega1))
				print('Polar Orbit')
				print(strength1p,strength2p,strength3p,phase1p,phase2p,phase3p)
				print(strength1s,strength2s,strength3s,phase1s,phase2s,phase3s)
				print('\n')
				fout1.write(str(M2)+'\t'+str(e1)+'\t'+str(P2)+'\t'+str(omega1)\
					+'\t'+str(strength1p)+'\t'+str(strength1p)+'\t'+str(strength2p)+'\t'+str(strength3p)+'\t'+str(phase1p)+'\t'+str(phase2p)+'\t'+str(phase3p)\
					+'\t'+str(strength1s)+'\t'+str(strength1s)+'\t'+str(strength2s)+'\t'+str(strength3s)+'\t'+str(phase1s)+'\t'+str(phase2s)+'\t'+str(phase3s)+'\n')


				sim=rebound.Simulation()
				sim.G=4*np.pi**2
				sim.add(m=M1)
				sim.add(m=M2,P=P1/365.2422,e=e1,omega=omega1,M=np.random.rand()*np.pi*2,inc=np.pi/2)
				o2,i2,O2=GeneratePlanetaryParameters(np.pi/2,0,omega1)
				sim.add(m=1/1047*0.220,P=P2/365.2422,e=0,omega=0,M=np.random.rand()*np.pi*2,inc=np.pi/2)
				R1=0.004649130334214449
				R2=MassToRadius(M2)*0.004649130334214449
				R=R1+R2
				sim.move_to_com()
				Period=P1/365.242199174
				time=np.linspace(0,P2/365.2422*10,4000)
				dyAB=0
				olddyAB=0
				ttvAB=[]
				ttvBA=[]
				transitABlist=[]
				transitBAlist=[]
				oldt=time[0]
				for i in range(len(time)):
					t=time[i]
					sim.integrate(t)
					p=sim.particles
					olddyAB=DyAB(oldt)
					dyAB=DyAB(t)
					#if i%10==0:
					#	print(t)
					#	pl.scatter(-p[0].x,-p[0].z,s=1,c='k')
					#	pl.scatter(-p[1].x,-p[1].z,s=1,c='k')
					#	pl.scatter(-p[2].x,-p[2].z,s=1,c='r')
					if dyAB>0 and olddyAB<0 and p[0].z<0 and np.abs(p[0].x)<0.05:
						transitAB=optimize.bisect(DyAB,oldt,t)
						transitABlist.append(transitAB)
						sim.integrate(transitAB)
					if dyAB>0 and olddyAB<0 and p[0].z>0 and np.abs(p[0].x)<0.1:
						transitBA=optimize.bisect(DyAB,oldt,t)
						transitBAlist.append(transitBA)
						sim.integrate(transitBA)
					oldt=t
				#pl.xlim([-1.5,1.5])
				#pl.ylim([-1.5,1.5])
				#pl.axis('equal')
				#pl.show()

				transitABlist=np.array(transitABlist)
				transitBAlist=np.array(transitBAlist)
				
				timing=transitABlist*365.2422
				index=range(len(transitABlist))
				try:
					r=optimize.curve_fit(harmonic,[timing,index],timing,[0,1/P2,0,P1,0,0,0,0,0,0],method='dogbox')
				except:
					print("Oops, cannot fit parameters from standard ETV model")
					#continue
					strength1p=0
					strength2p=0
					strength3p=0
					phase1p=0
					phase2p=0
					phase3p=0
				else:
					#for i in range(len(r[0])):
					#	print (str(r[0][i])+'+-'+str(np.sqrt(r[1][i][i])))
					T0,f,dT0,P,A1,B1,A2,B2,A3,B3=r[0]
					t=np.linspace(0,2,20000)
					model=A1*cos(2*pi*t)+B1*sin(2*pi*t)+A2*cos(4*pi*t)+B2*sin(4*pi*t)+A3*cos(6*pi*t)+B3*sin(6*pi*t)
					pl.scatter((timing-np.floor(timing*0.5*f)/0.5/f)*f,np.array(timing)-np.array(index)*P-T0-dT0,color='g')
					#pl.ylim([-0.0002,0.0002])
					pl.plot(t,model)
					pl.grid('on')
					pl.show()
					strength1p=sqrt(A1**2+B1**2)
					strength2p=sqrt(A2**2+B2**2)
					strength3p=sqrt(A3**2+B3**2)
					phase1p=arctan(A1/B1)+pi*np.sign(B1)
					phase2p=arctan(A2/B2)+pi*np.sign(B2)
					phase3p=arctan(A3/B3)+pi*np.sign(B3)

				timing=transitBAlist*365.2422
				index=range(len(transitBAlist))
				try:
					r=optimize.curve_fit(harmonic,[timing,index],timing,[0,1/P2,0,P1,0,0,0,0,0,0],method='dogbox')
				except:
					print("Oops, cannot fit parameters from standard ETV model")
					#continue
					strength1s=0
					strength2s=0
					strength3s=0
					phase1s=0
					phase2s=0
					phase3s=0
				else:
					#for i in range(len(r[0])):
					#	print (str(r[0][i])+'+-'+str(np.sqrt(r[1][i][i])))
					T0,f,dT0,P,A1,B1,A2,B2,A3,B3=r[0]
					t=np.linspace(0,2,20000)
					model=A1*cos(2*pi*t)+B1*sin(2*pi*t)+A2*cos(4*pi*t)+B2*sin(4*pi*t)+A3*cos(6*pi*t)+B3*sin(6*pi*t)
					pl.scatter((timing-np.floor(timing*0.5*f)/0.5/f)*f,np.array(timing)-np.array(index)*P-T0-dT0,color='g')
					pl.ylim([-0.0002,0.0002])
					pl.plot(t,model)
					pl.grid('on')
					pl.show()
					strength1s=sqrt(A1**2+B1**2)
					strength2s=sqrt(A2**2+B2**2)
					strength3s=sqrt(A3**2+B3**2)
					phase1s=arctan(A1/B1)+pi*np.sign(B1)
					phase2s=arctan(A2/B2)+pi*np.sign(B2)
					phase3s=arctan(A3/B3)+pi*np.sign(B3)
				print('Coplanar Orbit')
				print(strength1p,strength2p,strength3p,phase1p,phase2p,phase3p)
				print(strength1s,strength2s,strength3s,phase1s,phase2s,phase3s)
				print('\n')
				fout2.write(str(M2)+'\t'+str(e1)+'\t'+str(P2)+'\t'+str(omega1)\
					+'\t'+str(strength1p)+'\t'+str(strength1p)+'\t'+str(strength2p)+'\t'+str(strength3p)+'\t'+str(phase1p)+'\t'+str(phase2p)+'\t'+str(phase3p)\
					+'\t'+str(strength1s)+'\t'+str(strength1s)+'\t'+str(strength2s)+'\t'+str(strength3s)+'\t'+str(phase1s)+'\t'+str(phase2s)+'\t'+str(phase3s)+'\n')