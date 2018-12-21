load('/Users/apple/Desktop/SeniorThesis/GridSimulation/OmegapiP300.mat')
%e1grid=%reshape(e1,11,13);
%omega1grid=reshape(omega1,11,13);
P2grid=linspace(100,2000,96);
M2grid=linspace(0.1,1,19);
omega1grid=linspace(0,pi*2,13);
e1grid=linspace(0,0.5,11);
%%draw coplanar ratio 
s2pc=reshape(s2pc,11,19);
s2sc=reshape(s2sc,11,19);
ratioc=s2pc./s2sc;
figure();
imagesc(M2grid,e1grid,ratioc);
colorbar;
%%draw polar ratio 
s2pp=reshape(s2pp,11,19);
s2sp=reshape(s2sp,11,19);
ratiop=s2pp./s2sp;
figure();
imagesc(M2grid,e1grid,ratiop);
colorbar;
%%draw coplanar shift
ph2pc2=reshape(ph2pc,11,19);
ph2sc2=reshape(ph2sc,11,19);
shiftc=mod(ph2pc2-ph2sc2,2*pi);
shiftc=shiftc/2/pi;
figure();
imagesc(M2grid,e1grid,shiftc);
colorbar;
%%draw polar shift
ph2pp2=reshape(ph2pp,11,19);
ph2sp2=reshape(ph2sp,11,19);
shiftp=mod(ph2pp2-ph2sp2,2*pi);
shiftp=shiftp/2/pi;
figure();
imagesc(M2grid,e1grid,shiftp);
colorbar;