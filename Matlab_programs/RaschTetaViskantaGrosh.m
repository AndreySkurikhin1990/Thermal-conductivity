function g = RaschTetaViskantaGrosh()
format long g;
te0=273.15;
alfs=1e6*SredGraf(); 
knuSha=preobMasSha_n(); 
npp=Kramers_n();
dl=dlVo();
nsh=Kramers_n_Sha_loc(knuSha,dl);
ronu=RaschRonu(npp,nsh,dl,knuSha,alfs,0);
nsha=RaschRonu(npp,nsh,dl,knuSha,alfs,1);
nsham=chislInteg(dl,nsha)
Refl=chislInteg(dl,ronu)
np=chislInteg(dl,npp)
alf=chislInteg(dl,alfs)
%-------------
a=OprParam1mm(); koeftep=a(2); T0=a(3); Tk=a(4); y0=a(5);
Nto=1e2; l=y0/Nto; 
taum=0:l:y0; koom=taum; taum=alf*taum; 
Nto=length(taum); Tz=T0;
for k=1:Nto
    tetam(k)=((T0^4+(Tk^4-T0^4)*(k-1)/Nto)^(1/4))/Tz; 
end;
%-------------
int4tau=int34tau(0,4,taum);
int4tau0m=int34tau(1,4,taum);
int3taus=int34tau(0,3,taum); 
int3tau0mins=int34tau(1,3,taum); 
%-------------
n0=1; nite=1;
for j=1:nite
for k=1:Nto
    Tna=tetam(n0)*Tz; Tko=tetam(Nto-n0+1)*Tz;
    m=Gfunc_loc(koom(k),Tz,Tna/Tz,Tko/Tz,y0,alf,Refl, koeftep);
    r=integrTau_loc(koom(k),np,tetam,y0,Tz,Nto,alf,int3taus,int3tau0mins,k,koeftep);
    m=m+r;
    m=ProvAd(m);
    tetam(k)=real(m);
end; 
for g=1:Nto
tetam(g)=(tetam(g)+tetam(g+1))/2; 
end;
Tna=tetam(n0)*Tz; Tko=tetam(Nto-n0+1)*Tz; 
jvk=RaschTepPot(np,taum,y0,l,tetam,Tna,Tko,Tz,alf, int3taus,int4tau(Nto),int3tau0mins,Nto,koeftep,Refl);
ktk(j)=jvk*(koom(Nto-n0)-koom(n0)); 
ktk(j)=ktk(j)/abs(Tko-Tna); 
end;
Tna=tetam(n0)*Tz; 
Tko=tetam(Nto-n0+1)*Tz; 
%q=length(ktk); kt=trapz(ktk)/(q-1)
jv=RaschTepPot(np,taum,y0,l,tetam,Tna,Tko,Tz,alf,int3taus,int4tau(Nto+1),int3tau0mins,Nto,koeftep,Refl)
q=length(tetam); 
for k=1:q 
    tetam(k)=tetam(k)*Tz; 
    tetam(k)=tetam(k)-te0; 
end; 
kt=jv*abs(koom(q-n0)-koom(n0)); 
kt=kt/abs(tetam(q-n0+1)-tetam(n0))
p=plot(koom,tetam,'-b');
set(p,'LineWidth',2); hold on; grid on; 
xlabel({'Координата, мкм'}); 
ylabel({'Температура, °С'}); 
title({'График зависимости температуры от координаты'});
g=0;
end

function m = chislInteg(dl,nsha)
p=length(dl); su=0; 
for j=1:p-1 
    su=su+(dl(j+1)-dl(j))*(nsha(j)+nsha(j+1))/2; 
end
m=real(su)/(dl(p)-dl(1));
end

function [ dv ] = dlVo()
dl=dlvoVer53101; 
p=length(dl); 
for k=1:p 
    dl(k)=1e-2/dl(k); 
end;
dv=dl;
end

function [ rn ] = RaschRonu(npp,nsh,dl,knuSha,alfs,iden)
p=length(dl);
for k=1:p
    hiSha(k)=dl(k)*knuSha(k)/4/pi/nsh(k); 
    hi(k)=dl(k)*alfs(k)/4/pi/npp(k);
    nsha(k)=(nsh(k)^2+hiSha(k)^2)/(npp(k)^2+hi(k)^2);
    nsha(k)=nsha(k)^(1/2);
    ronu(k)=(1/2)+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2);
    ronu(k)=ronu(k)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); 
    ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));
    ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log((nsha(k)-1)/(nsha(k)+1))/((nsha(k)^2+1)^3);
end
if (iden>0)
    rn=nsha;
else
rn=ronu;
end
end

function p = ProvAd(m)
if (isnan(m))
    m=0;
elseif (isinf(m))     
    m=0;
elseif (abs(m)>1e0)   
    m=0;
elseif (m<0)   
    m=0;
end
p=m;
end

function [ nsh ] = Kramers_n_Sha_loc(alsr,Sp)
c0=299792458;
leSp=length(Sp);
for k=1:leSp
    nu(k)=2*pi*c0/Sp(k);
    hi(k)=Sp(k)*alsr(k)/4/pi;
end
enu=1e-4;
ma=izmMasChast(nu,enu); %определение массива nus на enu от шага nu
na=PokazPrelomAl(nu,ma,alsr,c0); %главная часть - определение ПП по КК
nnu=PokazPrelomPribl(nu,ma,na); %определение второго массива ПП
nsh=nnu-1+1.56;
end

function [ nus ] = izmMasChast(nu,enu)
lenu=length(nu);
for k=1:lenu-1
    nust(k)=enu*(nu(k+1)-nu(k))+nu(k);
end
nust(lenu)=2*nust(lenu-1)-nust(lenu-2);
nus=nust;
end

function [ pop ] = PokazPrelomAl(om,mao,ka,vl)
p=length(om);
for k=1:p
    q=1;
    for j=1:p-1
        dkpodo=(ka(j+1)-ka(j))/(om(j+1)-om(j));
        podln=(om(j)+mao(k))/(om(j)-mao(k));
        podln=abs(podln);
        fn(q)=dkpodo*log(podln);
        q=q+1;
    end
    fn(p)=fn(p-1);
    np(k)=1+(vl/pi)*integ2ma(om,fn)/2/mao(k);
end
pop=np;
end

function inte = integ2ma(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ ns ]  = PokazPrelomPribl(nu,nus,na)
nnu=0;
p=length(nu);
for k=1:p-1
    de=(nus(k)-nu(k))*(na(k+1)-na(k))/(nu(k+1)-nu(k));
    nnu(k)=de+na(k);
end
nnu(p)=(nus(p)-nu(p))*(na(p)-na(p-1))/(nu(p)-nu(p-1))+na(p);
ns = nnu;
end

function [ int4 ] = int34tau(mi,n,taum)
ep=1e-10; 
Nto=length(taum);
tau0=taum(Nto);
for k=1:Nto
    a=taum(k);
    y=0;
    if (mi>0)
        a=tau0-a; 
    end
    if (a>ep)
    y=abs(exp(-a)/a); 
    else
        tm=1/(n-1);
    end
    if (y<ep) 
        tm=0; 
    else
        tm=integroexpon(n,a); 
    end
    int4tau(k)=tm; 
end
int4=int4tau;
end

function [ r ] = oprKTPVer()
format longg; te0=273.15; %teh=[445 626 720]+te0; tec=[65 89 80]+te0; to=[50 60 70]*1e-3; tsp=(teh+tec)/2;
temvh=arrTemHigh207()+te0;
temvc=arrTemCold207()+te0; 
F=13.85; F=F*((1e-2)^2);
y0=3e1*1e-3;
qv=arrTepPot207()/F; 
ts=(temvc+temvh)/2;
for k=1:length(ts)
gT=abs(temvh(k)-temvc(k))/y0;
tepv(k)=qv(k)/gT;
end
tepv=tepv';
koefktp=polyfit(ts,tepv,2);
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

function alf = alfsBBIn(Tz,te0,Nto,dl)
no=BBIn(Tz-te0); j=1; 
tet=tempvtst(T0-te0,Tk-te0,y0,Nto); 
for k=no(1):no(2) 
    al(j)=alfs(k); 
    dlv(j)=dl(k); 
    j=j+1; 
end;
alf=chislInteg(dl,al);
end

function Gzn = Gfunc_loc(y,Tz,teta0,tetatau0,y0,alsred,Ref,koeftep)
sigma=5.668e-8; %te0=273.15; temv = [20 50 100 250 500]; temv = temv+273.15; tepv = [0.045 0.047 0.051 0.065 0.096]; tepv = tepv*4184/3600; kotepv = koef(temv,tepv,length(tepv))';koeftep = polyval(kotepv,Tz);
tau0=alsred*y0;
tau=alsred*y;
N=(koeftep*alsred)/(4*sigma*(Tz^3)); %Tr=1;for k=1:Ns Tr=Tr*exp(-alsred*tol); Tr=Tr*(1-Rfl); end;Tr=Tr*exp(-tau0);Ab=1-Tr;Ab=0;lam0=2898/(T0+te0);pi=3.1415926535897932;kap=(lam0*alsred)/(4*pi);Rfl=((nv-1)^2+kap^2)/((nv+1)^2+kap^2);eps=(0.6e-8)*4184/3600; stch=eps/sigma;Ns=y0/tol;
stch=1-Ref;
R0=stch*sigma*(teta0^4)+Ref*sigma*(tetatau0^4);
R0=R0*(Tz^4);
Rtau0=stch*sigma*(tetatau0^4)+Ref*sigma*(teta0^4);
Rtau0=Rtau0*(Tz^4);
tm=sigma*(Tz^4);
beta0=R0/tm;
betatau0=Rtau0/tm; %taus=alsred*y;
a=tau; 
tm=integroexpon(4,a); 
if (a==0)
    tm=1/3; 
end;
a=tau0; 
tr=integroexpon(4,a);
a=tau0-tau;
tu=integroexpon(4,a); 
G1 = beta0*(-tm+(tau/tau0)*tr+(1-tau/tau0)/3);
G2 = betatau0*((1-tau/tau0)*tr-tu+tau/(3*tau0));
G3 = 2*N*(teta0+(tau/tau0)*(tetatau0-teta0));
Gzn = (G1+G2+G3)/2/N;
end

function intgr = integrTau_loc(y,npp,teta,y0,Tz,Nt,alsr,int3taus,int3tau0Mins,nom,koeftep)
sigma=5.668e-8; %temv = [20 50 100 250 500]; %te0=273.15; %temv = temv+te0; %tepv = [0.045 0.047 0.051 0.065 0.096];%tepv = tepv*4184/3600; %kotepv = koef(temv,tepv,length(tepv))'; %koeftep = polyval(kotepv,Tz); 
N=(koeftep*alsr)/(4*sigma*(Tz^4));
tau0=alsr*y0; 
htaush=tau0/Nt; 
tau=alsr*y; 
for k=1:Nt+1
taush=(k-1)*htaush;
if (nom == k) 
    E3=1/2; 
else
    E3=-(integroexpon(3,abs(tau-taush))); 
end;
tm=int3taus(k); 
E3=E3+tm; 
tr=tm; 
tm=int3tau0Mins(k);
E3=E3+(tau/tau0)*(tm-tr); 
E3=(npp^2)*E3*(teta(q)^4); 
inte(q)=real(E3);
ko(q)=taush; 
q=q+1;
E3=0;
end;
integri=0; 
su=0; 
q=length(inte);
for k=1:q-1
su=su+(ko(k+1)-ko(k))*(inte(k)+inte(k+1))/2;
end
integri=su;
integri=integri/(2*N);
%disp(integri);
intgr=integri;
end