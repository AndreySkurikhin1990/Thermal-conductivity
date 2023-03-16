format long g; teta=0; te0=273.15; T0=1000+te0; tol=3e4; y0=1e3; l=2e2;
npp=Kramers_n(); dl=dlvoVer53101; p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k); 
end;
temvc=[109 235 109 238]; temvc=temvc+te0;
temvh=[585 1000 603 1000]; temvh=temvh+te0;
tepo=[3.9596 8.6377 3.5614925 7.123]; qv=(1e4/13.85)*tepo;
for k=1:4
tepv(k)=qv(k)*(tol*1e-6)/(temvh(k)-temvc(k));
end
temvs=(1/2)*(temvc+temvh);
tepv1=(tepv(1)+tepv(3))/2;
tepv2=(tepv(2)+tepv(4))/2;
tem1=(temvs(1)+temvs(3))/2;
tem2=(temvs(2)+temvs(4))/2;
k1=(tepv2-tepv1)/(tem2-tem1);
k2=tepv2-k1*tem2; 
Tz= 1e2+1:5e0:1e3-1;
Tz=Tz+te0;
p=length(Tz); koeftep=0;
for k=1:p
    koeftep(k) = abs(k1*Tz(k)+k2); 
end
sigma=5.668e-8; lamt=0; lamizl=0;
for k=1:p
    alf=1/sprko(Tz(k));
    np=real(nsreddv(dl,npp,Tz(k)));
    t=abs(16*(np^2)*sigma*(Tz(k)^3)/(3*alf*1e6));
    lamizl(k)=t;
    lamt(k)=koeftep(k)-lamizl(k);
    disp(k);
end
pl=plot(Tz-te0,lamizl,'-b',Tz-te0,lamt,'-k',Tz-te0,koeftep,'-c');
set(pl,'LineWidth',3); hold on; grid on; 
xlabel({'Температура, °С'}); 
ylabel({'Составляющие коэффициента теплопроводности, Вт/(м*К)'}); 
title({'График зависимости составляющих коэффициента теплопроводности от температуры'});
legend('Излучательная компонента', 'Кондуктивная часть КТП', 'Общий КТП');