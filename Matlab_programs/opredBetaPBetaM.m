function [ bplbmi ] = opredBetaPBetaM(T0,Tk,Tz,d,koeftep1,koeftep2,ts1,ts2)
Tna=T0
Tko=Tk
delT=-1e-1;
Te=T0:delT:Tk;
wmg=217e-3; wal=132e-3; wsi=368e-3;
koal=rasKoefAlpha(te,wmg,wsi,wal);
%kktp1=(koeftep2-koeftep1)/(ts2-ts1); 
%kktp2=koeftep1-kktp1*ts1; 
%k1=abs(T0-Tk)/d; 
%qv=koeftep1*k1;
%tol=1e-3
tol=30e-3
%Nte=1e2;
%te=raTePo(qv,T0,Nte,d,kktp1,kktp2)
%te=raTePoLin(T0,Tk,Nte,tol)
%Ntol=ceil(d/tol);
%koeftep3=kktp1*T0+kktp2;
%koeftep3=kktp1*(Tna+Tko)/2+kktp2;
%koeftep3=koeftep1;
%k1=qv/koeftep3;
%Tk=T0-abs(k1*tol)
%Tk=ts1
 s=0; n=length(Te);
for k=1:n
    s=s+koal(k);
end
s=s/n;
dko=0.444905; dko2=65e-2; dko3=1.3; s=s*dko*dko2*dko3;
sig=5.67e-8;
alf=s*sredkopoPatch(tol,Tna,Tko)
%reiz=resuIzlu(tol,T0,Tk,wmg,wsi,wal);
alf=knusreddvvsrVer((Tk+T0)/2)
crp=koeftep1*alf/4/sig/Tz^3
tau0=alf*tol
Nt=5e0*ceil(tau0)
tau00=tau0/Nt
Ntau=1e2*ceil(tau0)
htau=tau0/Ntau
taus=0:htau:tau0;
tetatau0=Tk/Tz;
teta0=T0/Tz;
tetataus=0;
E04m=0; E2=0; E3m=0; E2m=0; E01m=0;
npr=nsreddvVer((Tk+T0)/2);
for k=1:length(taus)
E2(k)=ProvAd(integroexpon(2,taus(k)));
tetataus(k)=teta0+taus(k)*(tetatau0-teta0)/tau0;
%np(k)=nsreddv(dv,npp,tetataus(k)*Tz);%np(k)=nsreddvvsr(dv,npp,tetataus(k)*Tz);
np(k)=npr;
E01m(k)=(np(k)^2)*(tetataus(k)^4);
E04m(k)=E01m(k)*ProvAd(integroexpon(1,abs(tau00-taus(k))));
E2a(k)=E2(k)*E01m(k);
E2m(k)=(ProvAd(integroexpon(2,tau0-taus(k)))+E2(k))*E01m(k);
E2(k)=(2-E2(k))*E01m(k);
E2t(k)=ProvAd(integroexpon(2,tau0-taus(k)))*E01m(k);
E3m(k)=ProvAd(integroexpon(3,tau0-taus(k)))*E01m(k);
end
tetatau00=teta0+tau00*(tetatau0-teta0)/tau0
npr=nsreddvVer(tetatau00*Tz);
E01=(npr^2)*(tetatau00^4)
E02=(-1/2)*ProvAd(integroexpon(2,tau00))
E03=(-1/2)*ProvAd(integroexpon(2,tau0-tau00))
E04=trapz(taus,E04m);
E01=-E01-E04
%--------
E20=trapz(taus,E2m)
E320=1/2-ProvAd(integroexpon(3,tau0))
E20=E20/E320;
%1--------
bpbm=reshenbpbm(1,1,E02,E03,E20,E01)
bp=bpbm(1)
bm=bpbm(2)
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau0,E2a,taus,E2t))
%--------
E300=integn2teta4doteta(np,tetataus,tau00,taus)
E301=(ProvAd(integroexpon(3,tau00))-(1/2))/2
E302=(ProvAd(integroexpon(3,tau0-tau00))-ProvAd(integroexpon(3,tau0)))/(-2)
E303=(-1/2)*trapz(taus,E2);
E304=(1/2)*integE2tmtsdot(np,tetataus,tau00,taus);
E300=-(E300+E303+E304)
%2-------
bpbm=reshenbpbm(1,1,E301,E302,E20,E300)
bp=bpbm(1)
bm=bpbm(2)
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau0,E2a,taus,E2t))
%3-------
bpbm=reshenbpbm(E02,E03,E301,E302,E20,E300)
bp=bpbm(1)
bm=bpbm(2)
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau0,E2a,taus,E2t))
%-------
E331=trapz(taus,E3m)/2/Nt;
E332=(-1/2)*intexpt3(tau00,np,tetataus,taus);
PO=Pervoobrn2teta4(np,tetataus,taus);
VPO=pervVtorRaz(taus,PO,tau00)
E331=-(E331+E332+VPO)
C11=((1/3)-ProvAd(integroexpon(4,tau0/Nt)))*(1-1/Nt)/2;
C12=(ProvAd(integroexpon(4,tau0))*(1-1/Nt)-ProvAd(integroexpon(4,tau0*(1-1/Nt)))+((1/3)/Nt))/2;
%4-------
bpbm=reshenbpbm(1,1,C11,C12,E20,E331)
bp=bpbm(1)
bm=bpbm(2)
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau0,E2a,taus,E2t))
%5-------
bpbm=reshenbpbm(E301,E302,C11,C12,E300,E331)
bp=bpbm(1)
bm=bpbm(2)
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau0,E2a,taus,E2t))
%6-------
bpbm=reshenbpbm(E02,E03,C11,C12,E01,E331)
bp=bpbm(1)
bm=bpbm(2)
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau0,E2a,taus,E2t))
%-------
izT0=epssr(T0)
izT0=izT0*sig*(Tz^4)
bplbmi=qr;
end

function mk = ProvAd(m)
ep=1e-40;
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
if (abs(m)<ep)  
    m=0;  
end; 
mk=real(abs(m));
end

function [ perobr ] = Pervoobrn2teta4(np,teta,tau)
Fp=0;
fu=0;
n=length(np);
taus=0;
for k=2:n
    for j=1:k
        fu(j)=(np(j)^2)*(teta(j)^4);
        taus(j)=tau(j);
    end
    Fp(k)=trapz(taus,fu);
    taus=0;
    fu=0;
end
Fp(1)=0;
perobr=Fp;
end

function [ vtoper ] = pervVtorRaz(tau,perPervoob,taukon)
VtPer=0;
perPer=0;
n=length(tau);
taum=0;
q=1;
for k=1:n
    if (tau(k)<taukon)
    taum(q)=tau(k);
    perPer(q)=perPervoob(k);
    q=q+1;
    end
end
VtPer=trapz(taum,perPer);
vtoper=VtPer;
end

function E32 = intexpt3(taukon,npp,tetam,taum)
E3a=0;
n=length(npp);
q=1;
taumm=0;
for k=1:n
    if (taum(k)<taukon)
    E3a(q)=ProvAd(integroexpon(3,taukon-taum(k)))*(npp(k)^2)*(tetam(k)^4);
    taumm(q)=taum(k);
    q=q+1;
    end
end
E3=trapz(taumm,E3a);
E32=E3;
end

function inte = integE2tmtsdot(np,teta,taukon,tau)
n=length(tau);
taus=0;
q=1;
integ=0;
for k=1:n
    if (tau(k)<taukon)
        integ(q)=ProvAd(integroexpon(2,taukon-tau(k)))*(np(k)^2)*(teta(k)^4);
        taus(q)=tau(k);
        q=q+1;
    end
end
inte=(1/2)*trapz(taus,integ);
end

function inte = integn2teta4doteta(np,teta,taukon,tau)
n=length(tau);
taus=0;
q=1;
integ=0;
for k=1:n
    if (tau(k)<taukon)
        integ(q)=(np(k)^2)*(teta(k)^4);
        taus(q)=tau(k);
        q=q+1;
    end
end
inte=trapz(taus,integ);
end

function [ luchtepl ] = rasqrlamtep(bp,bm,Tz,sig,tau0,E2,tau,E2t)
qr=0;
E3=ProvAd(integroexpon(3,tau0));
E2i=trapz(tau,E2);
E2ti=trapz(tau,E2t);
qr0=2*sig*(Tz^4)*(bp/2-bm*E3-E2i);
qrtau0=2*sig*(Tz^4)*(bp*E3-bm/2+E2ti);
qr(1)=qr0;
qr(2)=qrtau0;
luchtepl=qr;
end

function epsi = epssr(tem)
dv=RasshDiapDlinVoln();
npp=Kramers_n();
alf=RasMasKoAbs();
n=length(npp);
np=0;
hi=0;
ep=0;
for k=1:n
hi=dv(k)*alf(k)/4/pi;
np=(npp(k)^2+hi^2)^0.5;
ep(k)=emdiel(np);
%ep(k)=emmet(npp(k),hi);
end
epsi=epssred(dv,ep,tem,npp);
end

function sredpatch = sredkopoPatch(d,T0,Tk)
Ns=2e1;
h=d/Ns;
s=0:h:d;
alf=0;
Tn=T0;
k1=abs(T0-Tk)/d;
for k=1:length(s)
    te=Tn-h*k1;
alf(k)=ProvAd(koefpoglPatch(s(k),Tn));
Tn=te;
disp(k);
end
sredpatch=trapz(s,alf)/d
end

function [ ratepol ] = raTePo(qv,T0,Nt,tol,k1,k2)
ht=tol/Nt;
xl=0:ht:tol;
n=length(xl);
tem=T0; te=0; delt=0;
for k=1:n
    te(k)=tem;
    koeftep=k1*tem+k2;
    delt=abs(qv*ht/koeftep);
    tem=tem-delt;
end
ratepol=te;
end

function [ ratepol ] = raTePoLin(T0,Tk,Nt,tol)
ht=tol/Nt;
xl=0:ht:tol;
n=length(xl);
kn=abs(T0-Tk)/tol;
tem=T0; te=0; delt=0;
for k=1:n
    te(k)=tem;
    delt=ht*kn;
    tem=tem-delt;
end
ratepol=te;
end

function [ res ] = reshenbpbm(a11,a12,a21,a22,b1,b2)
matkobpbm(1,1)=a11; 
matkobpbm(1,2)=a12; 
matkobpbm(2,1)=a21; 
matkobpbm(2,2)=a22; 
matsvchbpbm(1)=b1; 
matsvchbpbm(2)=b2;
bpbm=(inv(matkobpbm)*matsvchbpbm')'
bpbm(1)=(b1*a22-b2*a12)/(a11*a22-a21*a12);
bpbm(2)=(b2*a11-b1*a21)/(a11*a22-a21*a12);
res=bpbm;
end

function [ rezizl ] = resuIzlu(tol,T0,Tk,wmg,wsi,wal)
npp=Kramers_n();
dl=RasshDiapDlinVoln();
p=length(dl);
Nt=2e1
alfs=RasMasKoAbs();
tau=0; 
hksi=tol/Nt; 
ksi=0:hksi:tol; 
k1=(Tk-T0)/tol;
k2=T0;
Nt=length(ksi);
for k=1:p
    dv(k)=dl(k)/npp(k);
end
for k=1:Nt
    temp(k)=k1*ksi(k)+k2;
end
koal=rasKoefAlpha(te,wmg,wsi,wal);
temp=1*temp'
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=2*pi*PP*(c0^2);
C2=PP*c0/PB;
%1----------
taut=0; I1g=0; I2g=0;
for j=1:Nt
    ksi1=0; 
    I11g=0;
for k=1:j
    eb1=0; it1=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
altek=alfs(m)*koal(k);
altej=alfs(m)*koal(j);
eb1(m)=altek*c1t/(lam^5)/me;
it1(m)=eb1(m)*ProvAd(integroexpon(2,altej*ksi(j)-altek*ksi(k)));
    end
    ksi1(k)=ksi(k);
    I11g(k)=integ2ma(dv,it1);
end
I1g(j)=2*pi*integ2ma(ksi1,I11g);
%2----------
    ksi2=0; I22g=0; q=1;
for k=j:Nt
    eb2=0; 
    it2=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
altek=alfs(m)*koal(k);
altej=alfs(m)*koal(j);
eb2(m)=altek*c1t/(lam^5)/me;
taut=altek*ksi(k)-altej*ksi(j);
it2(m)=eb2(m)*ProvAd(integroexpon(2,taut));
    end
    ksi2(q)=ksi(k);
    I22g(q)=integ2ma(dv,it2);
    q=q+1;
end
I2g(j)=-2*pi*integ2ma(ksi2,I22g);
end
%3----------
I3g=0;
for k=1:Nt
    eb3=0; 
    it3=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
alte=alfs(m)*koal(k);
eb3(m)=alte*c1t/(lam^5)/me;
taut=alte*ksi(k);
it3(m)=eb3(m)*ProvAd(integroexpon(2,taut));
    end
    I3g(k)=2*pi*integ2ma(dv,it3);
end
I3gg=integ2ma(ksi,I3g);
%4-----------
I41gm=0;
I42gm=0;
I4gm=0; 
for k=1:Nt
    it4=0;
    tem=temp(k);
    for m=1:p
        it4(m)=ProvAd(integroexpon(3,alte*ksi(k)));
    end
    I41gm(k)=2*pi*integ2ma(dv,it4);
    I42gm(k)=-pi;
    I4gm(k)=I41gm(k)+I42gm(k);
end
%5-----------
I51gp=0;
I52gp=0;
I5gp=0;
for k=1:Nt
    it51=0;
    it52=0; 
    tau0=0; 
    taut=0;
    tem=temp(k);
    for m=1:p
        alte=alfs(m)*koal(Nt);
        tau0=alte*ksi(Nt);
        it52(m)=ProvAd(integroexpon(3,tau0));
        alte=alfs(m)*koal(k);
        taut=alte*ksi(k);
        it51(m)=ProvAd(integroexpon(3,tau0-taut));
    end
    I51gp(k)=-2*pi*integ2ma(dv,it51);
    I52gp(k)=2*pi*integ2ma(dv,it52);
    I5gp(k)=I51gp(k)+I52gp(k);
end
%6-----------
I1q=0; 
for k=1:Nt
    eb6=0; 
    it6=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
eb6(m)=c1t/(lam^5)/me;
alte=alfs(m)*koal(k);
it6(m)=eb6(m)*alte;
    end
    I1q(k)=4*pi*integ2ma(dv,it6);
end
%7---dq/dt-----
I1dq=0;
for k=1:Nt
    eb7=0; 
    it7=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
eb7(m)=c1t/(lam^5)/me;
alte=alfs(m)*koal(k);
it7(m)=eb7(m)*alte;
    end
    I1dq(k)=4*pi*integ2ma(dv,it7);
end
%8---dq/dt-----
I41gmdq=0;
for k=1:Nt
    it8=0;
    for m=1:p
alte=alfs(m)*koal(k);
taut=alte*ksi(k);
it8(m)=ProvAd(integroexpon(2,taut));
    end
    I41gmdq(k)=-2*pi*integ2ma(dv,it8);
end
%9---dq/dt-----
I51gpdq=0;
for k=1:Nt
    it9=0;
    for m=1:p
        alte=alfs(m)*koal(Nt);
tau0=alte*ksi(Nt);
alte=alfs(m)*koal(k);
taut=alte*ksi(k);
it9(m)=ProvAd(integroexpon(2,tau0-taut));
    end
    I51gpdq(k)=-2*pi*integ2ma(dv,it8);
end
%10---dq/dt-----
ep=1e-5;
I12gdq=0;
for j=1:Nt
    I12ggdq=0;
    for k=1:Nt
    eb10=0; 
    it10=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
altej=alfs(m)*koal(j);
altek=alfs(m)*koal(k);
eb10(m)=altek*c1t/(lam^5)/me;
tauj=altej*ksi(j);
tauk=altek*ksi(k);
tauj=abs(tauj-tauk);
it10(m)=eb10(m)*ProvAd(integroexpon(1,tauj));
    end
    I12ggdq(k)=integ2ma(dv,it10);
end
I12gdq(j)=-2*pi*integ2ma(ksi,I12ggdq);
end
I12gdq=1*I12gdq'
%Total-----------
Ipl=0; Imi=0; Ite=0; q=1;
for k=2:Nt-1
    Ite=oprIpIm(I41gm, I1g, I2g, I51gp, k, Nt);
    Ipl(q)=Ite(1);
    Imi(q)=Ite(2);
    q=q+1;
end
qr=0; qr0=0; ksi4=0;
for k=1:length(Imi)
    qr0(k)=-(I3gg+I52gp(k)*Imi(k)+I42gm(k)*Ipl(k));
    qr(k)=I1q(k)-(I1g(k)+I2g(k)+I41gm(k)*Ipl(k)+I51gp(k)*Imi(k))-qr0(k);
    ksi4(k)=ksi(k);
end
qr=1*qr'
qr0=1*qr0'
Ipl=1*Ipl'
Imi=1*Imi'
I4gm=1*I4gm'
I5gp=1*I5gp'
I41gm=1*I41gm'
I52gp=1*I52gp'
I1q=1*I1q'
I1g=1*I1g'
I2g=1*I2g'
I3gg=1*I3gg
I1q=1*I1q
n=nsreddvVer((T0+Tk)/2)
Tz=T0;
disp((5.67e-8)*(n^2)*(Tz^4));
pl=plot(1e3*ksi4,qr,'-b');
set(pl,'LineWidth',2); 
hold on; grid on; 
xlabel({'Координата x, мм'}); 
ylabel({'Плотность теплового потока, Вт/м2'}); 
title({'График зависимости ПТП от x'});
rezizl=I3g';
end

function ie3m = integE3min(tau0,tau)
i1=-ProvAd(integroexpon(3,tau0-tau));
i2=ProvAd(integroexpon(3,tau0));
%i2=0;
ie3m=2*pi*(i1+i2);
end

function ie3p = integE3plu(tau)
i1=2*pi*ProvAd(integroexpon(3,tau));
i2=-pi;
%i2=0;
ie3p=(i1+i2);
end

function inte = integ2ma(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ mas ] = oprIpIm(I41gm, I1g, I2g, I51gp, k, Nt)
a11=I41gm(1)-I41gm(k);
a12=I51gp(1)-I51gp(k);
b1=-(I1g(1)-I1g(k)+I2g(1)-I2g(k));
a21=I41gm(1)-I41gm(Nt);
a22=I51gp(1)-I51gp(Nt);
b2=(I1g(Nt)-I1g(1)+I2g(Nt)-I2g(1));
de=a11*a22-a21*a12;
de1=b1*a22-b2*a12;
de2=a11*b2-a21*b1;
mas=[de1/de de2/de];
end

function [ ko ] = rasKoefAlpha(te,wmg,wsi,wal)
tem=3e2:1e2:16e2;
mgo = [1 0.94 0.87 0.801 0.736 0.676 0.635 0.59 ...
    0.567 0.543 0.53 0.525 0.515 0.507];
al2o3 = [1 0.98 0.946 0.898 0.854 0.797 0.753 ...
    0.709 0.676 0.642 0.618 0.598 0.582 0.57];
sio2 = [1 0.984 0.953 0.917 0.854 0.808 0.756 ...
    0.711 0.578 0.523 0.495 0.468 0.448 0.429];
wo=wmg+wsi+wal;
kumm=vychKoe(tem,te,mgo);
kuma=vychKoe(tem,te,al2o3);
kums=vychKoe(tem,te,sio2);
n=length(te);
for k=1:n
    koo(k)=kumm(k)*wmg/wo+kuma(k)*wal/wo+kums(k)*wsi/wo;
end
ko=koo;
end

function epsi = epsilam(n)
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(abs(n))/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log(abs(n-1)/abs(n+1))*((n^2-1)^2)/((n^2+1)^3);
epsi=eps;
end

function [ koe ] = vychKoe(tem,te,epra)
for q=1:length(te)
    f=1;
for k=2:length(tem)
    if (te(q)<tem(k))
        if (f>0)
        j=k;
        f=0;
        end
    end
end
ep(q)=epra(j-1)+(epra(j)-epra(j-1))*(te(q)-tem(j-1))/(tem(j)-tem(j-1));
end
koe=ep;
end