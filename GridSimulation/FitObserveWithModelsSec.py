import rebound
from orbital.elements import KeplerianElements
import numpy as np 
import matplotlib.pyplot as pl 
from scipy import optimize
from scipy.stats import linregress
from mpl_toolkits.mplot3d import Axes3D
from numpy import cos,sin,pi

#fread=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/FourierAnalysis/koi34.tt.dan.pri.trans')
'''
fread=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/ReboundCodes/TTVs.txt')
index=[]
timing=[]
uncertainties=[]
variables=[]
i=0
for line in fread.readlines():
	#line=fread.readlines()
	string=line.split()
	#print(string)
	#index.append(float(string[0]))
	index.append(i+1)
	#timing.append(float(string[1]))
	#timing.append(float(string[2]))
	
	uncertainties.append(0.000001)
	#uncertainties.append(float(string[2]))
	i+=1
'''
fread=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/FourierAnalysis/koi34.tt.dan.sec.trans')
#fread=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/FourierAnalysis/NoisedTTV.txt')
index=[]
timing=[]
uncertainties=[]
variables=[]
i=0
for line in fread.readlines():
	#line=fread.readlines()
	string=line.split()
	#print(string)
	index.append(float(string[0]))
	#index.append(i+1)
	timing.append(float(string[1]))
	#timing.append(float(string[0]))
	variables.append([float(string[1]),float(string[0])])
	#uncertainties.append(0.000001)
	uncertainties.append(float(string[2]))
	#i+=1
timing=[69.17983426,
124.7706616,
152.5661017,
180.3614993,
208.1568031,
235.9520152,
263.7474246,
291.5430842,
319.3383165,
347.1335515,
374.9289228,
402.724358,
430.5197998,
458.3152062,
486.110533,
513.9057617,
541.701077,
569.4967448,
597.2920715,
680.6780433,
708.4734872,
764.0642472,
791.8594983,
819.6547589,
875.2458261,
930.8363068,
958.6317315,
986.4271798,
1014.222603,
1042.017965,
1069.813239,
1125.403969,
1153.199578,
1180.994746,
1208.790017,
1236.58543,
1264.38088,
1292.176309,
1319.971682,
1375.562216,
1403.357607,
1431.153272,
1486.74372]
def harmonic(var,T0,f,dT0,P,A1,B1,A2,B2,A3,B3):
	#print(var)
	t,n=var
	#print(T0+P*n+dT0+A1*cos(2*pi*f*t)+B1*sin(2*pi*f*t)+A2*cos(4*pi*f*t)+B2*sin(4*pi*f*t)+A3*cos(6*pi*f*t)+B3*sin(6*pi*f*t))
	return T0+P*n+dT0+A1*cos(2*pi*f*t)+B1*sin(2*pi*f*t)+A2*cos(4*pi*f*t)+B2*sin(4*pi*f*t)+A3*cos(6*pi*f*t)+B3*sin(6*pi*f*t)#+A4*cos(8*pi*f*t)+B4*sin(8*pi*f*t)

r=optimize.curve_fit(harmonic,[timing,index],timing,[0,0.00348,0,30,0,0,0,0,0,0])#,uncertainties)
#print(r[0])
for i in range(len(r[0])):
	print (str(r[0][i])+'+-'+str(np.sqrt(r[1][i][i])))
timing=np.array([timing])

T0,f,dT0,P,A1,B1,A2,B2,A3,B3=r[0]
print(T0+dT0)
print(f)
print('\n')
t=np.linspace(0,2,20000)

P0=np.copy(P)

model=A1*cos(2*pi*t)+B1*sin(2*pi*t)+A2*cos(4*pi*t)+B2*sin(4*pi*t)+A3*cos(6*pi*t)+B3*sin(6*pi*t)#+A4*cos(8*pi*t)+B4*sin(8*pi*t)


pl.scatter((timing-np.floor(timing*0.5*f)/0.5/f)*f,np.array(timing)-np.array(index)*P-T0-dT0,color='g')
pl.ylim([-0.0002,0.0002])
pl.plot(t,model)
pl.grid('on')
pl.show()
#print(timing)
#timing=np.delete(timing,[21,31,33,48])
#print(timing)
#index=np.delete(index,[21,31,33,48])
timing0=(timing-np.floor(timing*0.5*f)/0.5/f)*f
timing1=np.array(timing)-np.array(index)*P-T0-dT0

for t in timing[0]:
	print(t)

fread=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/FourierAnalysis/koi34.tt.dan.sec.trans')
#fread=open('/Users/apple/Desktop/SeniorThesis/Fabrycky/FourierAnalysis/NoisedTTV.txt')
index=[]
timing=[]
uncertainties=[]
variables=[]
i=0
for line in fread.readlines():
	#line=fread.readlines()
	string=line.split()
	#print(string)
	index.append(float(string[0]))
	#index.append(i+1)
	timing.append(float(string[1]))
	#timing.append(float(string[0]))
	variables.append([float(string[1]),float(string[0])])
	#uncertainties.append(0.000001)
	uncertainties.append(float(string[0]))
	#i+=1



variables=np.reshape(variables,[len(index),2])


def harmonic2(var,T0,f,dT0,P,A1,B1,A2,B2,A3,B3):
	#print(var)
	t,n=var
	#print(T0+P*n+dT0+A1*cos(2*pi*f*t)+B1*sin(2*pi*f*t)+A2*cos(4*pi*f*t)+B2*sin(4*pi*f*t)+A3*cos(6*pi*f*t)+B3*sin(6*pi*f*t))
	return T0+P*n+dT0+A1*cos(2*pi*f*t)+B1*sin(2*pi*f*t)+A2*cos(4*pi*f*t)+B2*sin(4*pi*f*t)+A3*cos(6*pi*f*t)+B3*sin(6*pi*f*t)#+A4*cos(8*pi*f*t)+B4*sin(8*pi*f*t)

r=optimize.curve_fit(harmonic2,[timing,index],timing,[0,0.00348,0,30,0,0,0,0,0,0])#,uncertainties)
#print(r[0])
for i in range(len(r[0])):
	print (str(r[0][i])+'+-'+str(np.sqrt(r[1][i][i])))
timing=np.array([timing])



T0,f0,dT0,P=r[0][0:4]
print(T0+dT0)
print(1/f)

ttv=np.array(timing)-np.array(index)*P-T0-dT0
ttv=ttv[0]
#print(ttv)
timing=timing[0]

t=np.linspace(0,2,20000)

def Model(t):
	return (A1*cos(2*pi*(t))+B1*sin(2*pi*(t))+A2*cos(4*pi*(t))+B2*sin(4*pi*(t))+A3*cos(6*pi*(t))+B3*sin(6*pi*(t)))

model=Model(t)
residuals=(ttv-Model((timing-np.floor(timing*0.5*f)/0.5/f)*f))
print(np.mean((residuals/uncertainties)**2))
print(np.mean((timing1-ttv)**2/np.array(uncertainties)**2))
print('\n')
pl.errorbar((timing-np.floor(timing*0.5*f)/0.5/f)*f,ttv,yerr=uncertainties,fmt='.g',ms=10)
pl.scatter(timing0,timing1,color='r',s=10)
#pl.ylim([-0.0002,0.0002])
pl.plot(t,model)
pl.grid('on')
pl.show()


def harmonic3(var,t0,C):
	#print(var)
	t,n=var
	#print(T0+P*n+dT0+A1*cos(2*pi*f*t)+B1*sin(2*pi*f*t)+A2*cos(4*pi*f*t)+B2*sin(4*pi*f*t)+A3*cos(6*pi*f*t)+B3*sin(6*pi*f*t))
	t-=t0
	result=A1*cos(2*pi*f*t)+B1*sin(2*pi*f*t)+A2*cos(4*pi*f*t)+B2*sin(4*pi*f*t)+A3*cos(6*pi*f*t)+B3*sin(6*pi*f*t)
	#print(T0+P*n+dT0+A1*cos(2*pi*f*t)+B1*sin(2*pi*f*t)+A2*cos(4*pi*f*t)+B2*sin(4*pi*f*t)+A3*cos(6*pi*f*t)+B3*sin(6*pi*f*t))
	return result*C

r=optimize.curve_fit(harmonic3,[timing,index],ttv,[0,1])#,uncertainties)
#print(r[0])
for i in range(len(r[0])):
	print (str(r[0][i])+'+-'+str(np.sqrt(r[1][i][i])))
timing=np.array([timing])
t0,C=r[0]
t=np.linspace(0,2,20000)
timing=timing[0]

def Model(t):
	return C*(A1*cos(2*pi*(t-f*t0))+B1*sin(2*pi*(t-f*t0))+A2*cos(4*pi*(t-f*t0))+B2*sin(4*pi*(t-f*t0))+A3*cos(6*pi*(t-f*t0))+B3*sin(6*pi*(t-f*t0)))

model=Model(t)
residuals=(ttv-Model((timing-np.floor(timing*0.5*f)/0.5/f)*f))
print(np.mean((residuals/uncertainties)**2))
pl.errorbar((timing-np.floor(timing*0.5*f)/0.5/f)*f,ttv,yerr=uncertainties,ms=10,fmt='.g')
pl.ylim([-0.0002,0.0002])
pl.plot(t,model)
pl.grid('on')
pl.show()

'''
pl.subplot(2,1,1)
pl.scatter((timing-np.floor(timing*0.5*f)/0.5/f)*f,Model((timing-np.floor(timing*0.5*f)/0.5/f)*f))
pl.errorbar((timing-np.floor(timing*0.5*f)/0.5/f)*f,ttv,yerr=uncertainties,ms=10,fmt='.g')
pl.ylim([-0.0002,0.0002])
pl.subplot(2,1,2)
pl.errorbar((timing-np.floor(timing*0.5*f)/0.5/f)*f,(ttv-Model((timing-np.floor(timing*0.5*f)/0.5/f)*f))/uncertainties,yerr=1,fmt='.g')
pl.show()
'''



