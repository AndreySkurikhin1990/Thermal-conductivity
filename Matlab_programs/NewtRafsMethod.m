function nf = NewtRafsMethod()
por = 0.6635; %2-0.7 mm
v=1; d=15e-6; N=20; te0=273.15; %zpo=0.17;  por=por+zpo;
%zpo=0.17; m22=(0.3+0.45)/2;
dan=vyborTempTepl(v,d,N);
la_voz=opredTeploprovVozd(T);
la_e=dan(3);
ts=dan(2)-te0;
%eps=epsisred(ts+te0)
eps=0.899851588244671;
%q_r=dan(1);
%tns=tn;
%la_co=(DulnevKoefTep(tns)+DulnevKoefTep207(ts, dan(4), dan(5), dan(6), lae1))/2;
%q_p=dan(1);
%d=30e-6;
%PoiKorn(epsSr(tn,d), opredTeploprovVozd(tn), dan(2), la_e, q_p);
%nf=PloPotIzl(1);
%nf = opredKoeTepSrInSp(la_e,la_voz,por);
nf=opredDulnLam1(por,ts+te0,eps,la_voz,la_e);
end

function [ tem ] = RasTempZer(Temp,d,N)
c=Temp(1); h=30e-3; M=0;
M(1,1)=h^2; M(1,2)=h; M(2,1)=(h/2)^2; M(2,2)=h/2; b=0;
b(1)=Temp(3)-c; b(2)=Temp(2)-c;
y=0; y=inv(M)*b'; a=y(1); b=y(2);
xn=h/2-N*d; xk=h/2+N*d;
Tn=a*(xn^2)+b*xn+c; Tk=a*(xk^2)+b*xk+c;
tem = [Tn Tk a b c];
end

function [ la ] = vyborTempTepl(vyb,d,N)
te0=273.15;
T1 = [(585*2+6e2)/3   (377+(383*2+384*2)/4+396)/3 ((120*3+119)/4+129+138)/3];
T1=T1+te0;
T2 = [1e3   ((697*3+698)/4+703)/2 ((261*3+262)/4+273)/2];
T2=T2+te0;
T3 = [8e2   548 2e2];
T3=T3+te0;
ts=[(T1(1)+T1(3))/2 (T3(1)+T3(3))/2 (T2(1)+T2(3))/2];
t_in_cal_600 = [(28.7*2+28.75*2)/4  (28.2*3+28.25)/4    (28.25*3+28.3)/4];
t_in_cal_1000 = [29.7   29.7];
t_in_cal_800 = 29.2;
t_out_cal_600 = [(29.15*2+29.2*2)/4  (28.65*3+28.7)/4    (28.7*3+28.75)/4];
t_out_cal_1000 = [30.65 30.65];
t_out_cal_800 = 29.9;
delta_t_600 = (t_out_cal_600 - t_in_cal_600);
delta_t_1000 = (t_out_cal_1000 - t_in_cal_1000);
delta_t_800 = (t_out_cal_800 - t_in_cal_800);
m_tau_600 = 1e-3 * [120 120 120] / 60;
m_tau_1000 = 1e-3 * [122 120] / 60;
m_tau_800 = 1e-3 * 120 / 60;
h=30e-3;
s_600=0;
t_600=0;
for k = 1:3
Q_600(k) = 4.19e3 * m_tau_600(k) * delta_t_600(k);
lam_600(k)=Q_600(k)*h/13.85e-4/(-T1(3)+T1(1));
s_600 = s_600+lam_600(k)*T1(k);
t_600=t_600+T1(k);
end
s_600=s_600/t_600;
s_1000=0;
t_1000=0;
for k = 1:2
Q_1000(k) = 4.19e3 * m_tau_1000(k) * delta_t_1000(k);
lam_1000(k)=Q_1000(k)*h/13.85e-4/(-T2(3)+T2(1));
s_1000=s_1000+lam_1000(k)*T2(k);
t_1000=t_1000+T2(k);
end
s_1000=s_1000/t_1000;
Q_800 = 4.19e3 * m_tau_800 * delta_t_800;
s_800=Q_800*h/13.85e-4/(-T3(3)+T3(1));
q=0; tsr=0; da=0;
switch (vyb) 
    case 1
    q = s_600*(T1(1)-T1(3))/h;
    tsr=(T1(1)+T1(3))/2;
    da=RasTempZer(T1,d,N);
    lam=s_600;
    case 2
    q = s_1000*(T2(1)-T2(3))/h;
    tsr=(T2(1)+T2(3))/2;
    da=RasTempZer(T2,d,N);
    lam=s_1000;
    case 3
    q = s_800*(T3(1)-T3(3))/h;
    tsr=(T3(1)+T3(3))/2;
    da=RasTempZer(T3,d,N);
    lam=s_800;
end
format long g;
ko=koefOprTemp(ts(1),ts(2),ts(3),s_600,s_800,s_1000);
la = [q,tsr,lam,ko(1),ko(2),ko(3),da(1),da(2),da(3),da(4),da(5)];
end

function [ tem ] = RasTemVnut()
tem=0;
end

function kor = PoiKorn(Aef, lamvo, tn, dlv, jv)
sig=5.668e-8;
ta=tn-50;
tb=tn+50;
eps=1;
h=1;
while (eps>1e-7)
tc=(ta+tb)/2;    
fc=sig*Aef*(tn^4-tc^4)+lamvo*(tn-tc)/dlv-jv;
fa=sig*Aef*(tn^4-ta^4)+lamvo*(tn-ta)/dlv-jv;
fb=sig*Aef*(tn^4-tb^4)+lamvo*(tn-tb)/dlv-jv;
if ((fa*fc)<0)
elseif ((fb*fc)>0)
    tb=tc; 
end
if ((fa*fc)>0) 
elseif ((fb*fc)<0) 
    ta=tc; 
end
eps=abs(fa-fb);
h=h+1;
if (h>1e2) 
    break; end
end
kor=tc;
end

function es = epsSr(T,d)
npp=Kramers_n();
dv=RasshDiapDlinVoln();
alsr=skleiv2Mas(soglMas(1e6*SredGraf(),1e6*SredGraSt()),1e6*SredGraf());
for k=1:length(npp)
    ronu(k)=dependRefln(npp(k));
    eps(k)=1-exp(-d*alsr(k))-ronu(k);
end
es = epssred(dv,eps,T,npp);
end

function [ obma ] = skleiv2Mas(ar1,ar2)
q=1; obm=0;
for k=1:(length(ar1)+length(ar2))
    if (k<(length(ar1)+1))
        obm(k)=ar1(k);
    else
        obm(k)=ar2(q);
        q=q+1;
    end 
end
obma=obm;
end

function [ dlv ] = soglMas(ar1,ar2)
dl=0; dl=dlvoi(); 
dlvv=0; dlvv=dlvoVer53101();
kon=dlvv(1);
ar3=0; p=length(dl); 
q=1; arr=0;
for k=1:p  
    if (dl(k)>kon)
        q=q+1;
    end
end
t2=ar2(q);
t1=ar1(1);
t1=(t1-t2);
j=1;
for k=1:q+p
    if (k<(q+1))
    arr(k)=ar2(k)+t1;
    end
end
dlv = arr;
end

function qr = PoiskHo(vy,d,N)
kpp=KoefPoglPla();
t=0; da=0;
t=1e2:1e-1:11e2;
da=vyborTempTepl(vy,d);
m=length(da);
a=da(m-2);
b=da(m-1);
c=da(m);
h=30e-3;
xn=h/2-N*d;
xk=h/2+N*d;
xi=xn+d/2:d:xk-d/2;
p=length(xi);
te=0;
for k=1:p
    te(k)=a*(xi^2)+b*xi+c;
end
pkp=kopo(kpp,t,te);
qr=0;
end

function [ koe ] = koefOprTemp(x1,x2,x3,y1,y2,y3)
A=0; b=0; A(1,1)=x1^2; A(1,2)=x1; A(1,3)=1;
A(2,1)=x2^2; A(2,2)=x2; A(2,3)=1;
A(3,1)=x3^2; A(3,2)=x3; A(3,3)=1;
b(1)=y1; b(2)=y2; b(3)=y3;
koe = inv(A)*b';
end

function [ al ] = kopo(kpp,t,te)
n=length(t);
m=length(te);
alp=0;
for q=1:m
    f=1;
for k=1:n
    if (t(k)>te(q))
        if (f>0)
            a=(kpp(k)-kpp(k-1))/(t(k)-t(k-1));
            b=kpp(k-1)-a*t(k-1);
            alp(q)=a*te(q)+b;
            f=0;
            break;
        end
    end
end
end
al=alp;
end

function [ qrr ] = PloPotIzl(vt)
if (vt==6e2)
    qr = [];
end
if (vt==8e2)
    qr = [];
end
if (vt==1e3)
    qr = [];
end
qrr=1*qr';
end

function ak = opredKoefAkkomodLandau(T)
PP=6.6260755e-34/2/pi;
kB=1.380658e-23;
NA=6.0221409e23;
R=kB*NA;
mu=29e-3;
a=3e-9;
gamv=7/5;
m=mu/NA;
ro=opredPlotnVozd(T-273.15);
c=sqrt(gamv*R*T/mu);
alpha=(a*c)^2;
alpha=(kB*T)/alpha;
alpha=alpha^(3/2);
alpha=alpha/6/ro;
alpha=alpha/sqrt(2*pi*m);
%alpha=(kB*T/PP/c)^3;
%alpha=1.7*m*alpha/ro;
ak=alpha;
end

function lam4 = DulnevSloFor(m2,T,eps,lavo,lam1)
r=(2+0.7)/2/2; d=2*r;
Nk=sqrt(m2^2-10*m2+9);
Nk=(m2+3+Nk)/2/m2;
%y2=3.3e-3*(1-m2)^(-2/9);
%pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2);
%y2=y2*(pud+9.8*ro1*(1-m2)*hsl)^(1/3);
%y1=1e-20; y2=y1;
eta=1e-4; 
y1=10e-4; y2=y1/sqrt(eta); %y1=(10+50)*1e-4/2;
y3=2*sqrt(Nk-1)/Nk;
y4=y3/((1-m2)^(1/3));
%hsh=2e-3/2;
hsh=1e-20;
%hshsr=hsh/r;
fb=opredFiBol(y3/y2,y2);
A=y2^2-y1^2;
F=sqrt(1-y2^2);
D=sqrt(1-y3^2);
E=y4^2-y3^2;
delsrsz=d*(hsh+1/Nk);
lamszm=MoleSostTeplVozd(T,delsrsz,lavo);
epspr=eps/(2-eps);
sig=5.668e-8;
lamszl=4*epspr*sig*delsrsz*(T^3);
lamsz=lamszm+lamszl;
nusz=lamsz/lam1;
w=(lavo/lamszm-nusz*D)/(1-nusz);
y12=(y1^2)/(hsh/2+(1-hsh/2)*fb);
DF=abs(w-D)/abs(w-F);
nug=nusz;
numz=MoleSostTeplVozd(T,hsh*r,lavo)/lam1;
nu2sp=nug;
DF=(D-F+w*log(DF))*2*nug/(1-nug);
AF=A/(1-hsh/2-F+hsh/2/numz);
ADF=1/(AF+DF);
ADF=1/(D/(y3^2)+ADF);
ADF=ADF+nu2sp*E+y12;
lam4=ADF*lam1/(y4^2);
end

function la = MoleSostTeplVozd(T,de,lamg)
gam=7/5;
H=101325;
H0=1e5;
kB=1.380658e-23;
n=H/kB/T;
d=(0.6*0.21+0.65*0.79)*1e-10;
sig=pi*(d^2)/4;
dli=1/sqrt(2)/n/sig;
Pr=opredPrVozd(T);
%a=opredKoefAkkomodLandau(T);
a=0.9;
cz=8.42e-3;
Ty=113;
la0=cz/H0/(1+Ty/T);
Kn=la0*H0/H/de;
B=4*gam/(gam+1)*(2-a)/a*Kn/Pr;
la=lamg/(1+B);
end

function  fb = opredFiBol(bb,rr)
nb=2e1;
sn=0;
for k=0:nb
    n=2*k+1;
    I1b=besseli(1,n*pi*bb);
    I1r=besseli(1,n*pi*rr);
    K1b=besselk(1,n*pi*bb);
    K1r=besselk(1,n*pi*rr);
   sn=sn+I1b/n^2/I1r*(I1r*K1b-K1r*I1b); 
end
fb=1-16*sn/pi^2;
end

function [ vy ] = VybCAB(fa,fb,fc,xa,xb,xc)
if ((fa*fc)<0)
if ((fb*fc)>0)
    xb=xc; 
end
end
if ((fa*fc)>0) 
if ((fb*fc)<0) 
    xa=xc; 
end
end
vy = [xa xb xc];
end

function ol1 = opredDulnLam1(po,T,eps,lavo,lae)
lada=1e5;
ladb=1e-5;
ra=1e2;
ep=1e-7;
h=0; la=0;
while (ra>ep)
ladc=(lada+ladb)/2;    
fa=DulnevSloFor(po,T,eps,lavo,lada);
fa=fa-lae;
fb=DulnevSloFor(po,T,eps,lavo,ladb);
fb=fb-lae;
fc=DulnevSloFor(po,T,eps,lavo,ladc)
fc=fc-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
if (h>1e2) 
    break; end
end
ol1=ladc;
end

function es = epsisred(T)
npp=Kramers_n();
dv=RasshDiapDlinVoln();
lenp=length(npp);
for k=1:lenp
    eps(k)=epsilnu(npp(k));
end
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*c0^2;
c2=PP*c0/PB;
iz0=0; iz1=0; dl=0;
c1t=c1; c2t=c2;
for j=1:lenp
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
iz0(j)=c1t/((la^5)*(exp(c2t/(la*T))-1));
iz1(j)=eps(j)*iz0(j);
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
es=chi/zna;
end