format long g;
te0=273.15; T0=585+te0; Tk=109+te0; Tz=T0; 
temva = [25 300 600]; temva = temv+te0; tepva = [0.039 0.127 0.147]; kotepva = koef(temva,tepva,length(tepva))'; koeftepa = polyval(kotepva,Tz)
a=114e-3; F=13.85e-4; d=(4*F/pi)^0.5; h=30e-3;
%temv = [20 50 100 250 500]; temv = temv+te0; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))'; koeftep = polyval(kotepv,Tz)
%до обжига
te1=(583+111)/2+(580+127)/2+(571.5+124)/2+(603+101)/2
te2=(1000+199)/2+(1000+291)/2+(1000+262)/2+(1000.8+266)/2
la1=(0.106+0.103+0.1+0.147)/4 
la2=(0.145+0.187+0.203+0.227)/4
lam=0; tn=1e2; tk=1e3; Nt=5e3; h=(tk-tn)/Nt; te=tn; Ntr=1e4;
ko=(la2-la1)/(te2-te1); koo=0:1:3e4; koo=koo*1e-6;
for j=1:Nt
    lam(j)=la1+(te-te1)*ko;
    t(j)=te;
    te=te+h;
end
p=plot(t,lam,'-b');
set(p,'LineWidth',3); 
hold on; grid on; 
xlabel({'Температура, °С'}); 
ylabel({'Коэффициент теплопроводности, Вт/(м2*К)'}); 
title({'График зависимости lambda (T)'});