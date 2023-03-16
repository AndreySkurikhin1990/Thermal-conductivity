npp=(1.538+1.542)/2; y0=1000; l=1;
Nto=round(y0/l); Nit=100;
te0=273.15; T0=510+te0; Tk=490+te0; Tz=(T0+Tk)/2; alfs=IntegrDlVol; teta=0; tem=0; koo=0;
for k=1:Nto+1     teta(k)=((T0^4+(Tk^4-T0^4)*(k-1)/Nto)^0.25)/Tz; tem(k)=teta(k); koo(k)=(k-1)*l; end; tau=alfs*koo;
int4tau(1)=1/3; int4tau0m(Nto+1)=1/3;
    for k=2:Nto+1  tm=integroexpon(4,tau(k)); int4tau(k)=tm; int4tau0m(Nto+2-k)=tm; end;
int3taus(1)=1/2; int3tau0Mins(Nto+1)=1/2;
   for k=2:Nto+1  tm=integroexpon(3,tau(k)); int3taus(k)=tm; int3tau0Mins(Nto+2-k)=tm; end;
%disp(tau); disp(int3taus);
for i=1:Nit    q=1; m=0; r=0; y=0; disp(i);
for k=1:Nto+1
m = Gfunc(y, Tz, T0/Tz, Tk/Tz, y0, T0, Tk, npp, l,alfs,int4tau(k),int4tau(Nto+1),int4tau0m(k));
r = integrTau(y, npp, teta, y0, Tz, Nto,alfs,int3taus,int3tau0Mins, k);
m = m + r ;
if (isnan(m))     y=y+l; k=k+1;    continue; elseif (isinf(m))     y=y+l;  k=k+1;   continue; 
elseif (abs(m)>2)     y=y+l;   k=k+1;  continue;
else    teta(q)=real(m);     koo(q)=y;     q=q+1;     m=0;     r=0;
for g=1:Nto         teta(g)=(teta(g)+teta(g+1))/2;     end; end; y=y+l;     
end; end;
%q=length(koo);
for k=1:Nto+1  if (teta(k)<0) || (abs(teta(k))>1.4) teta(k)=tem(k); end;     teta(k)=teta(k)*Tz; end;
m=real(RaschTepPot(y0,l,teta,T0,Tk,Tz));
%koo=sort(koo); teta=-sort(-teta);
teta=teta-te0;
%for k=1:q disp(koo(k)); end;for k=1:q disp(teta(k)); end;
k=Nto+1;
disp(real(m*((koo(k-Nto/10)-koo(Nto/10))*(1e-6))/(-teta(k-Nto/10)+teta(Nto/10))))
p=plot(koo(Nto/10:k-Nto/10),teta(Nto/10:k-Nto/10),'-b');set(p,'LineWidth',3); hold on; grid on; xlabel({'Координата, мкм'}); ylabel({'Температура'}); title({'График зависимости температуры от координаты'});