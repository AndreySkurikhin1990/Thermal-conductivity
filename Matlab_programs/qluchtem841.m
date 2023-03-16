format long g; teta=0; te0=273.15; y0=3e4*1e-6; y0l=5e3*1e-6;%temvh=[558 950 560 950 545 927 530 925 600 1000 587 1000 590 1000 540 920 540 920 580 1000 580 1000 560 900]; temvh=temvh+te0;%temvc=[149 350 149 357 168 369 168 369 107 242 111 240 109 233 166 375 171 384 123 294 127 291 184 396]; temvc=temvc+te0;%tepo=[4.7766 11.5016 3.55732 11.4997 3.9805 9.5532 3.6447 9.4779 3.5144 8.6593 3.352 9.218 3.0313 7.7946 3.4023 10.2068 3.92812 11.17333 2.977 8.0448 3.0671 6.1342 4.4624 11.6021]; 
alfs=0; alfs=1e6*SredGraf(); npp=0; npp=Kramers_n(); te0=273.15; ep=1e-20; 
temvh=arrTemHigh(); temvh=temvh+te0; temvc=arrTemCold(); temvc=temvc+te0; 
tepo=arrTepPot84(); qv=(1e4/13.85)*tepo; pt=length(temvc); ts=sredtemp(temvh,temvc,pt); dl=dlivoln(); tepv=koeftepv(qv,y0,temvh,temvc,pt);
temvs=(1/2)*(temvc+temvh);
%tepv1=0; tepv2=0; tem1=0; tem2=0; for k=1:p if (rem(k,2)==1) tem1=tem1+temvs(k); tepv1=tepv1+tepv(k); else tem2=tem2+temvs(k); tepv2=tepv2+tepv(k); end; end; tem1=tem1/(p/2);tem2=tem2/(p/2); tepv1=tepv1/(p/2); tepv2=tepv2/(p/2); 
km=danPoTemTepl6(temvs,tepv);
k1=km(1);
k2=km(2);
%tem1=temvc(1); tem2=temvh(2); for k=1:p     if (tem2<temvh(k)) tem2=temvh(k);     end;    if (tem1>temvc(k))        tem1=temvc(k);    end; end;
Tz=2e2+1:2e0:1e3-1; Tz=Tz+te0;
p=length(Tz); koeftep=0;
for k=1:p
    koeftep(k) = abs(k1*Tz(k)+k2); 
end
sigma=5.67e-8; lamizl=0; lamte=0;
for k=1:p
    alf=1e6/sredRosSieg(Tz(k));
    np=nsreddv(dl,npp,Tz(k));
    t=abs(16*(np^2)*sigma*(Tz(k)^3)/(3*alf));
    lamizl(k)=t;
    lamte(k)=koeftep(k)-t;
    disp(k);
end
p=plot(Tz,lamizl,'-b',Tz,lamte,'-k',Tz,koeftep,'-c');
set(p,'LineWidth',3); hold on; grid on; 
xlabel({'Температура, К'}); 
ylabel({'Составляющие коэффициента теплопроводности, Вт/(м*К)'}); 
title({'График зависимости составляющих коэффициента теплопроводности от температуры'});
legend('Излучательная компонента', 'Кондуктивная часть КТП', 'Общий КТП');