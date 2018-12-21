%load('/Users/apple/Desktop/SeniorThesis/GridSimulation/Omegapim0.5.mat')
%e1grid=%reshape(e1,11,13);
%omega1grid=reshape(omega1,11,13);
P2grid=linspace(100,2000,96);
M2grid=linspace(0.1,1,19);
omega1grid=linspace(0,pi*2,13);
e1grid=linspace(0,0.95,20);
%%draw coplanar ratio 
s2pc=reshape(s2pc,20,96);
s2sc=reshape(s2sc,20,96);
ratioc=s2pc./s2sc;
figure();
ratioc(find((ratioc>100)+(ratioc<0.1)+(ratioc==NaN)))=NaN;
b=imagesc(P2grid,e1grid,log(ratioc));
caxis([-3 3]);
set(b,'AlphaData',~isnan(ratioc))
colorbar;
%%draw polar ratio 
s2pp=reshape(s2pp,20,96);
s2sp=reshape(s2sp,20,96);
ratiop=s2pp./s2sp;
figure();
ratiop(find((ratiop>100)+(ratiop<0.1)+(ratioc==NaN)))=NaN;
b=imagesc(P2grid,e1grid,log(ratiop));
caxis([-3 3]);
set(b,'AlphaData',~isnan(ratioc))
colorbar;
%%draw coplanar shift
ph2pc2=reshape(ph2pc,20,96);
ph2sc2=reshape(ph2sc,20,96);
shiftc=mod(ph2pc2-ph2sc2,2*pi);
shiftc(find(shiftc==0))=NaN;
shiftc=0.5-abs(shiftc/2/pi-0.5);
figure();

b=imagesc(P2grid,e1grid,shiftc);
caxis([0 0.5])
set(b,'AlphaData',~isnan(ratioc))

colorbar;
%%draw polar shift
ph2pp2=reshape(ph2pp,20,96);
ph2sp2=reshape(ph2sp,20,96);
shiftp=mod(ph2pp2-ph2sp2,2*pi);
shiftp(find(shiftp==0))=NaN;
shiftp=0.5-abs(shiftp/2/pi-0.5);
figure();
b=imagesc(P2grid,e1grid,shiftp);
caxis([0 0.5])
set(b,'AlphaData',~isnan(ratioc))
colorbar;