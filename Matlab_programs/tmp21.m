function tm = tmp21()
te0=273.15;
no=15;
y0=30e-3;
temvh=arrTemHigh(); temvh=temvh+te0; 
temvc=arrTemCold(); temvc=temvc+te0; 
ts=(temvc+temvh)/2;
tepo=arrTepPot84(); 
qv=(1/13.85e-4)*tepo; 
pt=length(temvc); 
tepv=koeftepv(qv,y0,temvh,temvc,pt);
T0=temvh(no);
Tk=temvc(no);
tol=1e-3
k1=abs(T0-Tk)/y0; 
ktp1=tepv(no);
qv=ktp1*k1;
ktp1=tepv(15);
ktp2=tepv(16);
ts1=ts(15);
ts2=ts(16);
kktp1=abs(ktp2-ktp1)/abs(ts2-ts1); 
kktp2=ktp1-kktp1*ts1;
koeftep3=kktp1*T0+kktp2
k1=qv/koeftep3;
Tk=T0-abs(k1*tol)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
npp=Kramers_n();
dl=RasshDiapDlinVoln();
p=length(dl);
Nt=1e1
alfs=RasMasKoAbs();
hksi=tol/Nt;
ksi=0:hksi:tol;
C1=2*pi*PP*(c0^2)
C2=PP*c0/PB
sig=2*C1*(pi^5)/(15*(C2^4))/(2*pi)
Nt=length(ksi);
k1=(Tk-T0)/tol;
k2=T0;
Nt=length(ksi);
for k=1:p
    dv(k)=dl(k)/npp(k);
end
for k=1:Nt
    temp(k)=k1*ksi(k)+k2;
end
temp=1*temp'
%1-----------
I1q=0;
for k=1:Nt
eb6=0;
eb6m=0;
tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
eb6m(m)=c1t/(lam^5)/me;
eb6(m)=alfs(m)*eb6m(m);
    end
    %eb6m=1*eb6m
    I1q(k)=4*pi*integ2ma(dv,eb6m);
end
I1q=1*I1q'
%3-----------
nu=0;
I1qq=0;
for j=1:Nt
    Ib=0;
    Ibm=0;
    tem=temp(j);
for k=1:p
    vl=c0/npp(k);
    nu(k)=vl/dv(k);
me=exp(PP*nu(k)/PB/tem)-1;
Ibm(k)=2*PP*(nu(k)^3)/(vl^2)/me;
Ib(k)=alfs(k)*Ibm(k);
end
for k=1:floor(p/2)
    me=nu(k);
    nu(k)=nu(p-k+1);
    nu(p-k+1)=me;
end
I1qq(j)=4*pi*trapz(nu,Ibm);
end
I1qq=2*I1qq'
%3-----------
inu=0; 
zsb=0;
n=0;
for j=1:Nt
    tem=temp(j);
    Ib=0;
for k=1:p
    vl=c0/npp(k);
me=exp((PP*nu(k))/(PB*tem))-1;
Ib(k)=2*PP*(nu(k)^3)/(vl^2)/me;
end
inu(j)=trapz(nu,Ib);
n(j)=nsreddvVer(tem);
    zsb(j)=inu(j)/sig/(tem^4)/n(j)^2;
    zsb(j)=zsb(j)*pi;
end
%4------------
I3g=0;
zsbl=0;
for k=1:Nt
    eb3=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
eb3(m)=c1t/(lam^5)/me;
    end
    I3g(k)=integ2ma(dv,eb3);
    zsbl(k)=I3g(k)/(sig*(tem^4)*(n(k)^2));
end
I3g=1*I3g'
%------------
inu=1*inu'
zsb=1*zsb'
pl=plot(temp-te0,I1q,'-b',temp-te0,I1qq,'-k');
set(pl,'LineWidth',2); 
hold on; grid on; 
xlabel({'Температура, К'});
ylabel({'Коэффициент СБ'}); 
title({'График зависимости КСБ от температуры'});
legend('I1(lambda)', 'I1(nu)');
tm=0;
end

function inte = integ2ma(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end