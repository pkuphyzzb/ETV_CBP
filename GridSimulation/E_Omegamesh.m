%load('/Users/apple/Desktop/SeniorThesis/GridSimulation/M0.5P300finer.mat')
%%e1grid=%reshape(e1,11,13);
%omega1grid=reshape(omega1,11,13);
P2grid=linspace(100,2000,96);
M2grid=linspace(0.1,1,19);
omega1grid=linspace(0,pi*2,16);
e1grid=linspace(0,0.95,21);
%%draw coplanar ratio 
s2pc=reshape(s2pc,20,16);
s2sc=reshape(s2sc,20,16);
ratioc=s2pc./s2sc;
ratioc(find((ratioc==0)+(ratioc==inf)+(ratioc==NaN)+(ratioc>100)))=NaN;
figure();
%subplot(2,2,1)
colormap('jet'); 

b=imagesc(omega1grid,e1grid,log(ratioc));
caxis([-3 3]);
set(b,'AlphaData',~isnan(ratioc))

colorbar;

%%draw polar ratio 
s2pp=reshape(s2pp,20,16);
s2sp=reshape(s2sp,20,16);
ratiop=s2pp./s2sp;
ratiop(find((ratiop==0)+(ratiop==inf)+(ratiop==NaN)+(ratiop>100)))=NaN;
figure();
%subplot(2,2,2)
colormap('jet'); 
b=imagesc(omega1grid,e1grid,log(ratiop));
caxis([-3 3]);
set(b,'AlphaData',~isnan(ratiop))
colorbar;
%%draw coplanar shift
ph2pc2=reshape(ph2pc,20,16);
ph2sc2=reshape(ph2sc,20,16);
shiftc=mod(ph2pc2-ph2sc2,2*pi);
shiftc(find(shiftc==0))=NaN;
shiftc=0.5-abs(shiftc/2/pi-0.5);
figure();
%subplot(2,2,3)
colormap('cool')
b=imagesc(omega1grid,e1grid,shiftc);
caxis([0 .5])
set(b,'AlphaData',~isnan(shiftc))

colorbar;
%%draw polar shift
ph2pp2=reshape(ph2pp,20,16);
ph2sp2=reshape(ph2sp,20,16);
shiftp=mod(ph2pp2-ph2sp2,2*pi);
shiftp(find(shiftp==0))=NaN;
shiftp=0.5-abs(shiftp/2/pi-0.5);
figure();
%subplot(2,2,4)
colormap('cool');
b=imagesc(omega1grid,e1grid,shiftp);
caxis([0 .5])
set(b,'AlphaData',~isnan(shiftp))

colorbar;