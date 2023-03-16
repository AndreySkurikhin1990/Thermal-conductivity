function rkv = rasKeramoVerm()
a=OprParam1mm(); 
tsr=a(1); koeftep=a(2); Tna=a(3); Tko=a(4); y0=a(5); q=a(6);
kopove=knusreddvvsrVer(tsr,1);
np=nsreddvVer(tsr);
R1=Refl_Shamot(Tna);
R2=Refl_Shamot(Tko);
tau0=kopove*y0;
RasTePo=rasTemPolvKeramover(tau0,koeftep,R1,R2,np,kopove,Tna,Tko,y0,q)
rkv=postgraf(0:tau0/(length(RasTePo)-1):tau0,RasTePo);
end

function [ r ] = oprKTPVer()
format longg; te0=273.15; %teh=[445 626 720]+te0; tec=[65 89 80]+te0; to=[50 60 70]*1e-3; tsp=(teh+tec)/2;
temvh=arrTemHigh207()+te0;
temvc=arrTemCold207()+te0; 
F=13.85; F=F*((1e-2)^2);
y0=3e1*1e-3;
qv=arrTepPot207()/F; 
ts=(temvc+temvh)/2;
p=length(ts);
for k=1:p
gT=abs(temvh(k)-temvc(k))/y0;
tepv(k)=qv(k)/gT;
end
koefktp=polyfit(ts,tepv',2);
tsr=mean(ts);
tgs=mean(temvh);
ths=mean(temvc);
ktp=polyval(koefktp,tsr);
qs=ktp*abs(tgs-ths)/y0;
ktpg=polyval(koefktp,tgs);
r=[tsr,ktp,tgs,ths,y0,qs,ktpg];
end

function [ r ] = OprParam1mm()
a=oprKTPVer();
koeftep=a(7);
Tna=a(3); 
Tko=a(4); 
y0=a(5); 
q=a(6);
ksi=1e0*1e-3;
Tk=abs(Tna-Tko)*ksi/y0;
Tk=Tko-Tk;
tsr=(Tna+Tk)/2;
r=[tsr,koeftep,Tna,Tk,ksi,q];
end

function p = postgraf(tau,tem)
pl=plot(tau,tem(1,:),'-m',tau,tem(2,:),'-r',tau,tem(3,:),'-b',tau,tem(4,:),'-k');
set(pl,'LineWidth',3); 
hold on; grid on; 
xlabel({'Координата (безразмерная)'}); 
ylabel({'Температура, К'}); 
title({'График зависимости температуры от безразмерной координаты'});
legend('Начальное приближение','Первая итерация','Вторая итерация', 'Третья итерация');
p=0;
end