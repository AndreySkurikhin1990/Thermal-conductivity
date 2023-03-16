%Расчет по методу Виканта-Гроша (1962)
function znfu = RaschTetaViskantaGroshN()
alfs=1e6*SredGraf(); 
npp=Kramers_n_uk();
dl=dlinyvoln(npp);
y0=3e1*1e-3;
l=1e0*1e-3;
koo=0:l:y0;
Nt=length(koo);
Nit=5;
te0=273.15; 
vmi=1; %данные 2020 года
vsv=0; %выбор состояния - исходный
ktpv=koefkoeficteplopver(vmi,te0,vsv,0)
Tgv=koefkoeficteplopver(vmi,te0,vsv,1);
Thv=koefkoeficteplopver(vmi,te0,vsv,2);
Tz=Tgv;
npr=nsreddv(dl,npp,Tz);
Tsr=(Tgv+Thv)/2;
alfs=alfs*OslabKP(Tsr);
alsr=knusreddv(alfs,npp,dl,Tz);
tau=alsr*koo;
n=length(tau);
tau0=tau(n);
for k=1:Nt
    a=tau(k);
    tm=integroexpon(3,a);
    int3taus(k)=tm;
    a=tau0-a; 
    tm=integroexpon(3,a); 
    int3tau0mins(k)=tm; 
    teta(k)=oprTemTeta(Tgv,Thv,Tz,k,Nt);
end;
ra=1e2;
tetam=teta;
for j=1:Nit    
for k=1:Nt
m=Gfunc(Tz,Tgv/Tz,Thv/Tz,alsr,ktpv,tau0,tau(k));
r=integrTau(npr,teta,Tz,Nt,alsr,int3taus,int3tau0mins,ktpv,tau0,npp,alfs,y0,Tgv,Thv,dl,tau(k),tau);
m=PrAd(m+r);
tetam(k)=real(m);
end; 
for k=1:Nt-1
teta(k)=(tetam(k)+tetam(k+1))/2; 
end; 
ra=abs(trapz(tetam)-trapz(teta));
if (ra<1e-6)
    break;
end
tetam=teta'
jv=RaschTepPot(npr,tau,y0,teta*Tz,Tgv,Thv,ktpv,tau0,alfs,npp,dl);
end;
jv
teta
znfu=0;
end

function [ dlv ] = dlinyvoln(npp)
format long g;
dl=dlvoVer53101(); 
dv=dl;
p=length(dl); 
fid = fopen('Chastota_izlucheniya_ver.txt','w');
vl0=299792458;
vl=vl0;
for k=1:p  
    dl(k)=1e-2/dl(k);
    dl(k)=dl(k)/npp(k);
    vl=vl/npp(k); 
    fr(k)=vl/dl(k);
    vl=vl0;
    fprintf(fid,'%0.20f\n',fr(k));
end;%knuSha=preobMas(); ronu=0; nsh=0; nsha=0; nsh=Kramers_n_Sha(knuSha,dl); nv=0; hi=0; hiSha=0; np=0; nsham=0;%for k=1:p  hiSha(k)=dl(k)*knuSha(k)/(4*pi); hi(k)=dl(k)*alfs(k)/(4*pi); nsha(k)=((nsh(k)^2+hiSha(k)^2)/(npp(k)^2+hi(k)^2))^0.5;ronu(k)=0.5+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log((nsha(k)-1)/(nsha(k)+1))/((nsha(k)^2+1)^3);end;%p=length(dl); su=0; for j=1:p-1    su=su+(dl(j+1)-dl(j))*(nsha(j)+nsha(j+1))/2; end; %su=0; for j=1:p-1 su=su+(dl(j+1)-dl(j))*(ronu(j)+ronu(j+1))/2; end; Refl=su/(dl(p)-dl(1)); Refl=abs(real(Refl)); %su=0; for j=1:p-1   su=su+(dl(j+1)-dl(j))*(npp(j)+npp(j+1))/2; end; %su=0; for j=1:p-1 su=su+(dl(j+1)-dl(j))*(alfs(j)+alfs(j+1))/2; end; alf=su/(dl(p)-dl(1)); alf=real(alf);
fclose(fid);
dlv = dv;
end

function w = koefkoeficteplopver(vmi,te0,vsv,vvz)
F=13.85*1e-4;
htp=30*1e-3;
if (vmi==0)
ktpver=arrKTP_VVF2();
Ts=arrTem_VVF2()+te0;
kktpv=polyfit(Ts,ktpver,2);
elseif (vmi==1)
    ktpv=arrKTP_VVF1();
    Tg=arrTem1VVF1()+te0;
    Th=arrTem2VVF1()+te0;
    Tc=arrTem3VVF1()+te0;
    Ts=(Tg+Th)/2;
    kktpv=polyfit(Ts,ktpv,2);
elseif (vmi==2)
Th=arrTemCold207()+te0;
Tg=arrTemHigh207()+te0;
Ts=(Th+Tg)/2;
tepo=arrTepPot207();
for k=1:length(tepo)
    qv=tepo(k)/F;
    gT=abs(Tg(k)-Th(k))/htp;
    ktpv(k)=qv/gT;
end
if (vsv==0)
kktp1=danPoTemTepl2071(Ts,ktpv); %исходный
kktpv=kktp1;
elseif (vsv==1)
kktp2=danPoTemTepl2072(Ts,ktpv); %повторы
kktpv=kktp2;
elseif (vsv==2)
kktp3=danPoTemTepl2073(Ts,ktpv); %после обжига
kktpv=kktp3;
end
end
if (vvz==0)
tsr=mean(Ts);
w=polyval(kktpv,tsr);
elseif (vvz==1)
w=mean(Tg);
elseif (vvz==2)
w=mean(Th);
else
w=Ts;
end
end

function [ w ] = RasMatrTeta(Ntt,Nto)
for k=1:Ntt 
    koo(k,1)=ko; ko=ko+y0l; 
end;
for k=1:Ntt-1 
    koo(k,Nto+1)=koo(k+1,1); 
end; 
ko=0;
for k=1:Ntt 
    for j=2:Nto 
        ko=ko+l; 
        koo(k,j)=ko; 
    end; ko=koo(k,1); 
end; 
koo(Ntt,Nto+1)=y0; 
Tna=T0; Tko=Tk; Tsr=Tz; teta(1,1)=T0/Tz; teta(Ntt,Nto)=Tk/Tz;
for k=2:Ntt 
    teta(k,1)=((T0^4+(Tk^4-T0^4)*(k-1)/Ntt)^0.25)/Tz; 
end;
for k=1:Ntt-1 
    teta(k,Nto+1)=teta(k+1,1); 
end;
for j=1:Ntt 
    Tna=teta(j,1)*Tz; 
    Tko=teta(j,Nto+1)*Tz; 
    Tsr=Tz; 
    for k=2:Nto
    teta(j,k)=((Tna^4+(Tko^4-Tna^4)*(k-1)/Nto)^0.25)/Tsr; 
    end; 
end;
format longg; 
for k=1:Nto+1 
tetam(k)=teta(1,k); 
koom(k)=koo(1,k); 
end; %nsham=real(trapz(dl,nsha)/(dl(p)-dl(1)));
w=teta;
end

function Gzn = Gfunc(Tz,teta0,tetatau0,alsred,koeftep,tau0,tau)
sigma=5.67e-8; %temv = [20 50 100 250 500]; temv = temv+273.15; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))';koeftep = polyval(kotepv,Tz);
N=(koeftep*alsred)/(4*sigma*(Tz^3));
beta0=1/pi;
betatau0=1/pi;
a=tau; 
tm=integroexpon(4,a); 
a=tau0; 
tr=integroexpon(4,a);
a=tau0-tau; 
tu=integroexpon(4,a); 
G1 = beta0*(-tm+(tau/tau0)*tr+(1-tau/tau0)/3);
G2 = betatau0*((1-tau/tau0)*tr-tu+tau/3/tau0);
G3 = 2*N*(teta0+(tau/tau0)*(tetatau0-teta0));
Gzn = (G1+G2+G3)/(2*N);
end

function ns = nsreddv(dv,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=PP*c0^2;
c2=PP*c0/PB;
c01=c1;
c02=c2;
for k=1:length(npp)
    lambda=dv(k);
    lambda=dv(k)/npp(k);
    c1=c1/(npp(k)^2);
    c2=c2/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
    %Ib(k)=Ib(k)/(npp(k)^2);    
Ibn(k)=npp(k)*Ib(k);
c1=c01;
c2=c02;
end
nc=trapz(dv,Ibn);
ns=abs(real(nc/trapz(dv,Ib)));
end

function knus = knusreddv(knu,npp,dv,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c10=PP*c0^2;
c20=PP*c0/PB;
c1=c10;
c2=c20;
for k=1:length(npp)
    c1=c1/(npp(k)^2);
    c2=c2/npp(k);
    lambda=dv(k)/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
Ibc(k)=knu(k)*(npp(k)^2)*Ib(k);
Ibz(k)=(npp(k)^2)*Ib(k);
c1=c10;
c2=c20;
end
nc=trapz(dv,Ibc);
nz=trapz(dv,Ibz);
knus=nc/nz;
end

function w = postgraf(tetam)
p=plot(1e-3*koom(n0+1:q-n0),tetam(n0+1:q-n0),'-b');
set(p,'LineWidth',2); 
hold on; grid on; 
xlabel({'Координата, мм'}); 
ylabel({'Температура, °С'}); 
title({'График зависимости температуры от координаты'});
w=0;
end

function intgr = integrTau(npr,teta,Tz,Nt,alsr,int3taus,int3tau0Mins,koeftep,tau0,npp,knuVer,y0,T1,T2,dv,tauk,tau)
sigma=5.67e-8;
qo=koeftep*abs(T1-T2)/y0;
eps=epssr([T1,T2],dv,npp,knuVer);
koeftep=PoiskChisKTP(qo,y0,T1,T2,sigma,alsr,eps);
htaush=tau0/Nt;
for k=1:Nt
taush=tau(k);
E3=integroexpon(3,abs(tauk-taush)); 
tm=int3taus(k); 
E3=tm-E3; 
tr=tm; 
tm=int3tau0Mins(k);
E3=E3+(tauk/tau0)*(tm-tr); 
E3=(npr^2)*E3*(teta(k)^4); 
inte(k)=real(E3);
bk(k)=taush;
end;
su=trapz(bk,inte);
N=(koeftep*alsr)/(4*sigma*(Tz^4));
intgr=su/(2*N);
end

function mk = PrAd(m)
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
mk=abs(real(m));
end

function t = PoiskChisKTP(qo,L,T1,T2,sig,betav,eps)
a=0; b=1e9; ra=abs(b-a); ep=1e-4; q=0; ni=1e4;
while ((ra>ep) && (q<ni))
    c=(a+b)/2e0;
    fa=a*abs(T1-T2)/L+(sig/(3*betav*L/4+1/eps(1)+1/eps(2)-1))*abs(T1^4-T2^4)-qo;
    fb=b*abs(T1-T2)/L+(sig/(3*betav*L/4+1/eps(1)+1/eps(2)-1))*abs(T1^4-T2^4)-qo;
    fc=c*abs(T1-T2)/L+(sig/(3*betav*L/4+1/eps(1)+1/eps(2)-1))*abs(T1^4-T2^4)-qo;
if ((fc*fb>0) && (fa*fc<0)) 
    b=c; 
end
    if ((fc*fa>0) && (fb*fc<0)) 
        a=c; 
    end
    q=q+1; ra=abs(b-a);
end
t=c;
end

function [ epsi ] = epssr(tem,dv,npp,alf)
n=length(npp); p=length(tem);
np=0; hi=0; ep=0;
for j=1:p
    eps(j)=0;
for k=1:n
hi=dv(k)*alf(k)/4/pi/npp(k);
np=(npp(k)^2+hi^2)^0.5;
ep(k)=emdiel(np);
%ep(k)=emmet(npp(k),hi);
end
eps(j)=epssred(dv,ep,tem(j),npp);
end
epsi=eps;
end

function koefte = RaschTepPot(nv,tau,y0,tem,T0,Tk,koeftep,tau0,alf,npp,dv)
cons=abs(Tk-T0)/y0; jv=koeftep*cons;
sigma=5.67e-8;
stchm=epssr([T0,Tk],dv,npp,alf);
stch=stchm(1);
R0=stch*sigma*(T0^4);
stch=stchm(2);
Rtau0=stch*sigma*(Tk^4);
stch=mean(stchm);
ppi=denTherFlux(nv,tem,Tk,tau0,stch,sigma,tau,R0,Rtau0);
jv=jv+ppi;
koefte=jv;
end

function qtp = denTherFlux(nv,tem,Tk,tau0,stch,sigma,tau,R0,Rtau0)
E3tau0=integroexpon(3,tau0);
E4tau0=integroexpon(4,tau0);
qss=(1-stch)*E3tau0+E4tau0/tau0-1/3/tau0;
qss=R0*qss;
qssa=(-E4tau0)/tau0-1/2+1/3/tau0;
qssa=Rtau0*qssa;
qss=qss+qssa+(stch/2)*sigma*(Tk^4);
for k=1:length(tau)
taus=tau(k);
a=tau0-taus;
tm=integroexpon(2,a);
tr=integroexpon(3,a);
tu=integroexpon(3,taus);
inte=(nv^2)*((1-stch)*tm+(tr-tu)/tau0)*sigma*(tem(k)^4);
integ(k)=PrAd(inte);
end;
integr=2*(trapz(tau,integ)+qss);
integr=abs(real(integr));
qtp=integr;
end

function w = OslabKP(Tsr)
dko2m=[4.68,4.69,5.65,13.17,20.2,27.81]/1e2; 
dko2t=6e2:2e2:16e2; 
for k=1:length(dko2m)
    dko2m(k)=1-dko2m(k);
end
dko2=opredKTPTKTochVer(dko2m, dko2t, Tsr);
wmg=217*1e-3; wal=132*1e-3; wsi=368*1e-3;
koal=rasKoefAlpha([Tsr],wmg,wsi,wal); s=0; n=length(koal);
for k=1:n
    s=s+koal(k);
end
s=s/n;
dko=0.444905; dko3=1.3; s=s*dko*dko2*dko3;
w=s;
end

function [ ko ] = rasKoefAlpha(te,wmg,wsi,wal)
Tna=3e2; Tko=16e2; delT=1e2; tem=Tna:delT:Tko;
mgo = [1 0.94 0.87 0.801 0.736 0.676 0.635 0.59 ...
    0.567 0.543 0.53 0.525 0.515 0.507];
al2o3 = [1 0.98 0.946 0.898 0.854 0.797 0.753 ...
    0.709 0.676 0.642 0.618 0.598 0.582 0.57];
sio2 = [1 0.984 0.953 0.917 0.854 0.808 0.756 ...
    0.711 0.578 0.523 0.495 0.468 0.448 0.429];
wo=(wmg+wsi+wal);
kumm=vychKoe(tem,te,mgo);
kuma=vychKoe(tem,te,al2o3);
kums=vychKoe(tem,te,sio2);
n=length(te);
for k=1:n
    koo(k)=kumm(k)*wmg/wo+kuma(k)*wal/wo+kums(k)*wsi/wo;
end
ko=koo';
end

function [ koe ] = vychKoe(tem,te,epra)
n=length(tem);
m=length(te);
for q=1:m
    f=1; j=1;
for k=1:n
    if (te(q)<tem(k))
        if (f>0)
        j=k;
        f=0;
        end
    end
end
if ((j==1) && (f==0)) 
    j=2;
end
if ((j==1) && (f==1))
    j=n;
end
ep(q)=epra(j-1)+(epra(j)-epra(j-1))*(te(q)-tem(j-1))/(tem(j)-tem(j-1));
end
koe=ep;
end

function t = opredKTPTKTochVer(ktptks, te, temp)
n=length(te); f=1; p=0;
for k=1:n
if ((te(k)>=temp) && (f>0)) 
        p=k; f=0;
end
end
    if ((p==1) && (f==0)) 
        p=2; 
    end
    if ((f==1) && (p==1)) 
        p=n; f=0;
    end
    if (f==0)
	ko=(ktptks(p)-ktptks(p-1))/(te(p)-te(p-1)); 
    ktp=ktptks(p-1)+ko*(temp-te(p-1));
    else ktp=0; ko=0; 
    end
    if ((ktp<0) || (ktp>1))
        ktp=1;
    end
t=ktp;
end

function w = oprTemTeta(Tgv,Thv,Tz,k,Nt)
T=Tgv-abs(Tgv-Thv)*(k/Nt);
w=T/Tz;
end