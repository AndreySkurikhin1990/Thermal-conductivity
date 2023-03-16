l=5e-2;
lv=1e-3;
t1=65+273.15;
t0=445+273.15;
epsil=(0.6e-8)*4184/3600; 
temv = [20 50 100 250 500];
temv = temv+273.15;
tepv = [0.045 0.047 0.051 0.065 0.096];
tepv=tepv*4184/3600;
kotepv = koef(temv,tepv,length(tepv))';
tv=350:50:550;
tv=tv+273.15;
tepvo=[0.0491 0.0521 0.0548 0.0574 0.0598];
kotepvo=koef(tv,tepvo,length(tepvo))';
for k=1:5 
    polj(k)=kotepv(k)/(6-k);
end
polj(6)=0;
jt1=polyval(polj,t1)/l;
jv=polyval(polj,t0)/l-jt1;
lamt=[0.1383 0.1426 0.1631 0.1925];
%75:25:375 0.0495 0.0532 0.0569 0.0608 0.0647 0.0689 0.0735 0.0785 0.0841 0.0907 0.0986 0.1084 0.1211 
h=5;
n=length(lamt);
tv=0;
t0=t0-273.15;
tv=t0-rem(t0,h)-h*(n-2):h:t0-rem(t0,h);
tv(n)=t0;
t0=t0+273.15;
tv=tv+273.15;
ksi=(epsil/jv)^0.25;
te0=0;
te0=401:1:445;
te0=te0+273.15;
lam=polyval(kotepv,te0);
n=length(lam);
x0=0;
for k=1:n
x0(k)=-(lam(k)/(4*jv*ksi))*log(abs((1-ksi*te0(k))/(1+ksi*te0(k))))+0.5*lam(k)*atan(ksi*te0(k))/(jv*ksi);
end
dt=0.1;
te0=te0-dt;
lam=polyval(kotepv,te0);
x1=0;
for k=1:n
x1(k)=-(lam(k)/(4*jv*ksi))*log(abs((1-ksi*te0(k))/(1+ksi*te0(k))))+0.5*lam(k)*atan(ksi*te0(k))/(jv*ksi);
end
dx=x1-x0;
la=-jv*dx/dt;
tv=tv+dt-273.15;
te0=te0-273.15;
plot(te0,la,'-m',te0,lam,'-r');
grid on;
xlabel({'Температура, °С'});
ylabel({'Коэффициент теплопроводности, Вт/(м*К)'});
title({'Зависимость коэффициента теплопроводности от температуры'});