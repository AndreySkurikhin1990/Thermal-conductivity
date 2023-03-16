%Расчет экономического эффекта
function t = tmp68()
format longg;
sigma=5.67e-8;
te0=273.15;
T=1e3+te0;
tkomn=23+te0;
skna=5e2;
vrna=(T-tkomn)/skna;
nsic=5;
koe=1e-3;
r=4*koe;
h=65*koe;
h1=25*koe;
a=114*koe; 
b=a; l=a; 
lt=150*koe;
dT=1e2; 
hv=34*koe;
rmt=30*koe;
Qv=RasKonvVnesh(h,h1,a,hv,lt)
pout=oprConvSiC(r, dT, a, T, Qv, nsic, b, h, h1, l, 0); %Расчет теплового потока от источника
k=1; Qc=pout(k); 
k=length(pout); 
Tpo=pout(k); %Fvp=pout(k-1); %площадь внутренней поверхности стенки, кроме образца
Tsvo=pout(k-2); %температура воздуха
ptpv=pout(5); %ptpn=pout(4);
Qob=pout(k-4);
kpd=l/lt
tsena=5.92;
kvch=1e3*36e2;
vr=8*36e2; kona=vrna/vr; koi=1e0-kona;
y=365-53*2-9-1-1-1-1-1-1-1-1
chid=1e0;
Qnagr=RasQNagr(ptpv,Tpo,Tsvo,a,h,h1,hv,rmt)
stelen=(Qnagr*chid+koi*Qob*vr/kpd+kona*Qob*vr/kpd/2)*tsena*y/kvch
t=Qob;
end

function p = RasF21(xn,xk,h,arg,l,r,d)
x=xn;
dA1=h*l;
A2=2*pi*r*l;
A1=l*d;
F12=0;
while (x<=xk)
arg(4)=x;
Fd12=viewfactor28(arg);
F12=F12+Fd12*dA1;
x=x+h;
end
F12=F12/A1;
p=F12*A1/A2;
end

function v = viewfactor28(ARG)
H=ARG(1);
R=ARG(2);
S=ARG(3);
X=ARG(4);
HH=H/R;
SS=S/R;
XX=X/R;
C=SS^2+XX^2;
CC=sqrt(C);
A=HH^2+C-1;
B=HH^2-C+1;
v=(SS/C)*(1-(acos(B/A)-sqrt((A^2)+4*(HH^2))*acos(B/A/CC)/2/HH-B*asin(1/CC)/2/HH)/pi-A/HH/4);
end

function v = viewfactor38(ARG)
A=ARG(1);
B=ARG(2);
C=ARG(3);
if (abs(C)>0)
X=A/C;
Y=B/C;
else 
    X=0;
    Y=0;
end
eps=1e-15;
if ((X<eps) || (Y<eps))
v=0;
else        
RTX=sqrt(1+X^2);
RTY=sqrt(1+Y^2);
RT=sqrt(1+X^2+Y^2);
v=(log(RTX*RTY/RT)+X*RTY*atan(X/RTY)+Y*RTX*atan(Y/RTX)-X*atan(X)-Y*atan(Y))*2/(pi*X*Y);
end
end

function e = epsver(tem)
dv=RasshDiapDlinVoln();
npp=Kramers_n();
eps=epsilonnu(npp);
e=epssredVer(dv,npp,tem,eps);
end

function ns = epssredVer(dv,npp,tem,eps)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458; vl=c0; dl=0;
for k=1:length(npp)
    vl=vl/npp(k);
    c1=2*pi*PP*(vl^2);
    c2=PP*vl/PB;
    lambda=dv(k)/npp(k);
    Ib(k)=c1/(lambda^5)/(exp(c2/lambda/tem)-1);
    Ibn(k)=eps(k)*Ib(k);
    vl=c0;
    dl(k)=lambda;
end
nc=trapz(dl,Ibn);
nz=trapz(dl,Ib);
ns=abs(real(nc/nz));
end

function [ epsnu ] = epsilonnu(npp)
ep=0;
for k=1:length(npp)
n=npp(k);
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(n)/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log((n-1)/(n+1))*((n^2-1)^2)/((n^2+1)^3);
ep(k)=eps;
end
epsnu=ep;
end

function [ vykh ] = oprConvSiC(r, dT0, a, T, Qv, nsic, b, h, h1, l, vyb)
n=nsic+2;
d=2*r;
F=oprviewf(nsic, d, n, a, r, b, h);
A=oprPlo(a, b, h, h1, r, nsic, l);
eps=opreps(T, n, n-1);
dl=a; %dl=150*1e-3;
Fp=A(n);
qs=Qv/Fp;
ep=1e-3;
k=0;
kit=1e2;
ra=1e2;
Ta=T-dT0+ep;
Tb=T+dT0; %температура SiC
T0=T;
sigma=5.67e-8;
while ((ra>ep) && (k<kit))
Tc=(Ta+Tb)/2; %находим температуру нагревателя
dTc=(Tc-T0);
Tsc=(T0+Tc)/2;
tem=oprtem(nsic,Tsc,dTc,n);
qsc=oprConvSiCVs(A,sigma,tem,eps,F,n,1)*sign(dTc);
fc=RasdeltaT(h1,h,(Tsc+T0)/2,dTc,d,0,qs,qsc);
dTa=(Ta-T0);
Tsa=(Ta+T0)/2;
tem=oprtem(nsic,Tsa,dTa,n);
qsa=oprConvSiCVs(A,sigma,tem,eps,F,n,1)*sign(dTa);
fa=RasdeltaT(h1,h,(Tsa+T0)/2,dTa,d,0,qs,qsa);
dTb=(Tb-T0);
Tsb=(Tb+T0)/2;
tem=oprtem(nsic,Tsb,dTb,n);
qsb=oprConvSiCVs(A,sigma,tem,eps,F,n,1)*sign(dTb);
fb=RasdeltaT(h1,h,(Tsb+T0)/2,dTb,d,0,qs,qsb);
if (((fa*fc)<0) && ((fb*fc)>0))
    Tb=Tc; 
end
if (((fa*fc)>0) && ((fb*fc)<0))
    Ta=Tc; 
end
ra=abs(Ta-Tb);
k=k+1;
end
Tc1=Tc;
vm=oprConvSiC2(dT0, T0, nsic, h1, n, F, A, eps);
Tc2=vm(1);
vm=oprConvSiC3(dT0, T0, nsic, h1, n, F, A, eps);
Tc3=vm(1);
if (vyb==0)
Tc=mean([Tc1,Tc2,Tc3]);
Tc=max([Tc1,Tc2,Tc3])
else
end
dTc=Tc-T0;
Rts=RasdeltaT(h1,h,(Tc+T0)/2,dTc,d,1,qs,qsc);
ptpv=dTc/Rts;
Rts=RasdeltaT2(h1,(Tc+T0)/2,dTc,1,qs,qsc);
ptpn=dTc/Rts;
tem=oprtem(nsic,Tc,dTc,n);
Qist=oprConvSiCVs(A,sigma,tem,eps,F,n,3) %тепловой поток от нагревателя
ptpv=ptpv+oprConvSiCVs(A,sigma,tem,eps,F,n,1)*sign(dTc);
ptpn=ptpn+oprConvSiCVs(A,sigma,tem,eps,F,n,2)*sign(dTc);
Qobn=A(n-1)*ptpn;
Qobv=A(n)*ptpv;
Qob=Qobv+Qobn;
vykh=[Qob,Tc,dTc,ptpn,ptpv,(Qist+Qob)/2,Qist,(Tc+T0)/2,Fp,(Tc-dTc)];
end

function [ vykh ] = oprConvSiC2(dT0, T0, nsic, h1, n, F, A, eps)
dT0=2*dT0;
ktp=0.226439;
te0=273.15;
Tg=1e3+te0;
Tkh=226+te0;
hob=3e1*1e-3;
qs=ktp*abs(Tg-Tkh)/hob;
Fp=A(n-1);
ep=1e-3;
k=0;
kit=1e2;
ra=1e2;
Ta=T0-dT0+ep;
Tb=T0+dT0; %температура SiC
sigma=5.67e-8;
while ((ra>ep) && (k<kit))
Tc=(Ta+Tb)/2; %находим температуру нагревателя
dTc=(Tc-T0);
Tsc=(T0+Tc)/2;
tem=oprtem(nsic,Tsc,dTc,n);
qsc=oprConvSiCVs(A,sigma,tem,eps,F,n,2)*sign(dTc);
fc=RasdeltaT2(h1,(Tsc+T0)/2,dTc,0,qs,qsc);
dTa=(Ta-T0);
Tsa=(Ta+T0)/2;
tem=oprtem(nsic,Tsa,dTa,n);
qsa=oprConvSiCVs(A,sigma,tem,eps,F,n,2)*sign(dTa);
fa=RasdeltaT2(h1,(Tsa+T0)/2,dTa,0,qs,qsa);
dTb=(Tb-T0);
Tsb=(Tb+T0)/2;
tem=oprtem(nsic,Tsb,dTb,n);
qsb=oprConvSiCVs(A,sigma,tem,eps,F,n,2)*sign(dTb);
fb=RasdeltaT2(h1,(Tsb+T0)/2,dTb,0,qs,qsb);
if (((fa*fc)<0) && ((fb*fc)>0))
    Tb=Tc; 
end
if (((fa*fc)>0) && ((fb*fc)<0))
    Ta=Tc; 
end
ra=abs(Ta-Tb);
k=k+1;
end
Rts=RasdeltaT2(h1,(Tc+T0)/2,dTc,1,qs,qsc);
vykh=[Tc,dTc,(Tc+T0)/2,Fp,(Tc-dTc)];
end

function po = oprConvSiCVs(A,sigma,tem,eps,F,n,v)
nsic=n-2;
for k=1:n
    ho(k)=0;
    id(k)=1;
    pin(k)=sigma*(tem(k)^4);
end
pout=GRAYDIFF(1, A, eps, ho, F, id, pin);
if (v==1) %верхняя
qi=pout(n);
elseif (v==2) %нижняя
qi=pout(n-1);
end
s=0;
for k=1:nsic
    pout(k)=pout(k)*A(k);
    s=s+pout(k);
end
if (v==3)
    qi=s;
end
po=abs(qi);
end

function [ vykh ] = oprConvSiC3(dT0, T0, nsic, h1, n, F, A, eps)
qs=0;
ep=1e-3;
k=0;
kit=1e2;
ra=1e2;
Ta=T0-ep;
Tb=T0+3*dT0; %температура SiC
sigma=5.67e-8;
while ((ra>ep) && (k<kit))
Tc=(Ta+Tb)/2; %находим температуру нагревателя
dTc=(Tc-T0);
Tsc=(T0+Tc)/2;
tem=oprtem(nsic,Tsc,dTc,n);
qsc=oprConvSiCVs(A,sigma,tem,eps,F,n,1)*sign(dTc);
Rtsc=RasdeltaT2(h1,(Tsc+T0)/2,dTc,1,qs,qsc);
qc=dTc/Rtsc+qsc;
fc=oprRazTemp(qc,T0);
dTa=(Ta-T0);
Tsa=(Ta+T0)/2;
tem=oprtem(nsic,Tsa,dTa,n);
qsa=oprConvSiCVs(A,sigma,tem,eps,F,n,1)*sign(dTa);
Rtsa=RasdeltaT2(h1,(Tsa+T0)/2,dTa,1,qs,qsa);
qa=dTa/Rtsa+qsa;
fa=oprRazTemp(qa,T0);
dTb=(Tb-T0);
Tsb=(Tb+T0)/2;
tem=oprtem(nsic,Tsb,dTb,n);
qsb=oprConvSiCVs(A,sigma,tem,eps,F,n,1)*sign(dTb);
Rtsb=RasdeltaT2(h1,(Tsb+T0)/2,dTb,1,qs,qsb);
qb=dTb/Rtsb+qsb;
fb=oprRazTemp(qb,T0);
if (((fa*fc)<0) && ((fb*fc)>0))
    Tb=Tc; 
end
if (((fa*fc)>0) && ((fb*fc)<0))
    Ta=Tc; 
end
ra=abs(Ta-Tb);
k=k+1;
end
vykh=[Tc,dTc,qc];
end

function rT = oprRazTemp(qv,Tpo)
te0=273.15;
Tvnepo=7e1+te0;
hsh=3e1*1e-3; %толщина шупа
hsht1=55*1e-3; %толщина ШТ-1
kktpsh=polyfit([8,12,14]*1e2+te0,[16,23,28]*1e-2,1);
ktpsh=polyval(kktpsh,Tpo);
kktpsht1=polyfit([35,65]*1e1+te0,[5,6]*1e-1,1);
Tksh=Tpo-qv*hsh/ktpsh;
ktpsht1=polyval(kktpsht1,Tksh);
Tsht1=Tksh-qv*hsht1/ktpsht1;
rT=Tsht1-Tvnepo;
end

function vykh = RasdeltaT(ra1,ra2,T,dT,d,vy,qs,qsi)
sig=sign(dT);
dT=abs(dT);
g=9.8;
beta0=1/T;
PrVo=vychKoefVoz(1,T);
roVoz=vychKoefVoz(2,T);
muVoz=vychKoefVoz(3,T);
lambVoz=vychKoefVoz(4,T);
nuVo=muVoz/roVoz;
GrVo=beta0*dT*g*(d^3)/(nuVo^2);
RaVo=GrVo*PrVo;
if ((RaVo>1e-2) && (RaVo<1e2))
    NuVo=1.02*(RaVo^(0.15));
elseif ((RaVo>1e2) && (RaVo<1e4))
    NuVo=0.85*(RaVo^(0.19));
elseif ((RaVo>1e4) && (RaVo<1e7))
    NuVo=0.5*(RaVo^(0.25));
elseif ((RaVo>1e7) && (RaVo<1e10))
    NuVo=0.125*(RaVo^(0.33));
else NuVo=0;
end
alp=lambVoz*NuVo/d;
alp=abs(alp);
if (alp>0)
Rts=1/alp+(ra2+ra1)/lambVoz/2;
Rts=abs(Rts);
else Rts=0;
end
if (Rts>0)
qr=dT*sig/Rts;
qr=qr+qsi;
else qr=0;
end
if (vy==0)
vykh=(qr-qs);
else
vykh=Rts;
end
end

function vykh = RasdeltaT2(ra1,T,dT,vy,qs,qsi)
sig=sign(dT);
dT=abs(dT);
lambVoz=vychKoefVoz(4,T);
Rts=ra1/lambVoz;
Rts=abs(Rts);
if (Rts>0)
qtpvo=dT*sig/Rts;
qr=qtpvo+qsi;
else qr=0;
end
if (vy==0)
vykh=(qr-qs);
else
vykh=Rts;
end
end

function w = RasKonvVnesh(h,h1,ab,hv,htsic)
hsh=30*1e-3;
hsht1=55*1e-3;
hmet=10*1e-3;
c=(h+h1+hmet+hsh+hsht1);
a=c+hv;
b=abs(a-c);
te0=273.15;
Ts=7e1;
Ts=Ts+te0;
Tb=23;
Tb=Tb+te0;
T=(Ts+Tb)/2;
dT=abs(Tb-Ts);
g=9.8;
beta=1/T;
PrVo=vychKoefVoz(1,T);
roVoz=vychKoefVoz(2,T);
muVoz=vychKoefVoz(3,T);
lambVoz=vychKoefVoz(4,T);
nuVo=muVoz/roVoz;
alpb=RasAlp1(PrVo,dT,g,beta,b,nuVo,lambVoz);
alpc=RasAlp1(PrVo,dT,g,beta,c,nuVo,lambVoz);
htnsht1=140*1e-3;
htak=10*1e-3;
c1=(2*(hsh+hsht1)+htsic+2*hmet); %верхний
b1=(2*(hsh+hsht1)+ab+2*hmet); %верхний
a1=(2*htnsht1+2*hmet+2*htak+ab); %нижний %dl=150*1e-3; b=114*1e-3; 
Qb=alpb*b*4*a1; %вертикальный с торцов нижних
Qc=alpc*c*b1*4; %вертикальный с торцов верхних
Q1=RasAlp2(PrVo,dT,g,beta,a1,c1,nuVo,lambVoz)*abs((a1^2)-(b1*c1)); %горизонтальный нижний
Q2=RasAlp2(PrVo,dT,g,beta,a1,b1,nuVo,lambVoz)*(b1*c1); %горизонтальный верхний
w=(Qb+Qc+Q1+Q2);
end

function w = RasAlp1(PrVo,dT,g,beta,a,nuVo,lambVoz)
F=(0.437/PrVo)^(9/16);
F=(1+F)^(-16/9);
GrVoa=g*beta*dT*(a^3)/(nuVo^2);
NuVa=GrVoa*PrVo*F;
NuVa=0.67*(NuVa^(1/4));
w=NuVa*lambVoz*dT/a;
end

function w = RasAlp2(PrVo,dT,g,beta,a1,c1,nuVo,lambVoz)
p=2*(c1+a1);
F=a1*c1;
l=F/p;
GrV=g*beta*(l^3)*dT/(nuVo^2);
RaV=GrV*PrVo;
NuV=(0.322/PrVo)^(11/20);
if (RaV>1e5)
NuV=0.15*(RaV^(1/3))/(1+NuV)^(20/33);
else
NuV=0.766*(RaV^(1/5))/(1+NuV)^(4/11);    
end
w=NuV*lambVoz*dT/l;
end

function w = RasQNagr(qv,Tpo,Tsvo,a1,h1,h2,h,htp)
a2=a1; %a2=15e1*1e-3;
te0=273.15;
Tvnpo=7e1+te0;
kktpver=polyfit((arrTem1VVF1()+arrTem2VVF1())/2+te0,arrKTP_VVF1(),2);
%---вермикулит
Tkver=1e3+te0; Tnver=226+te0; Tsver=(Tkver+Tnver)/2;
lamver=polyval(kktpver,Tsver) %0.226439
ptpn=lamver*abs(Tkver-Tnver)/htp;
Ts=(Tnver+Tkver)/2;
Tn=23+te0; 
dT=Ts-Tn;
ro=25e1;
vv=(a1^2)*h;
kcpver=[-2e-4,0.2796,295.47];
cpver=(polyval(kcpver,Ts)+polyval(kcpver,Tn))/2;
mver=(ro*vv);
Qver=cpver*mver*dT
%---воздух---%h1 - до образца снизу--%h2 - до образца сверху
Vvo=(a1^2)*(h1+h2);
cpvo=vychKoefVoz(5,Tsvo);
rovo=vychKoefVoz(2,Tsvo);
mvo=(Vvo*rovo);
Qvo=cpvo*mvo*abs(Ts-Tn)
%---шуп и легковес
hsh=3e1*1e-3; %толщина шупа
hvsht1=55*1e-3; %толщина ШТ-1
kktpsh=polyfit([8,12,14]*1e2+te0,[16,23,28]*1e-2,1);
ktpsh=polyval(kktpsh,Tpo);
kktpsht1=polyfit([35,65]*1e1+te0,[5,6]*1e-1,1);
Tksh=Tpo-qv*hsh/ktpsh;
if (Tksh<Tvnpo)
    Tksh=Tvnpo;
end
ktpsht1=polyval(kktpsht1,Tksh);
cpsh=(84e1+96e1)/2; 
rosh=4e2;
Ssh=(2*hsh+a1)*(2*hsh+a2);
Vsh=(h1+h2)*(Ssh-a1*a2)+Ssh*hsh;
msh=(Vsh*rosh);
Qsh=msh*cpsh*abs((Tpo+Tksh)/2-Tn)
Tsht1=Tksh-qv*hvsht1/ktpsht1;
if (Tsht1<Tvnpo)
    Tsht1=Tvnpo;
end
a=1.033; b=0.175*1e-3; d=-0.301/1e-5;
Tsht1sr=(Tsht1+Tksh)/2;
cpk=a+Tsht1sr*b+d/(Tsht1sr^2);
cpn=a+Tn*b+d/(Tn^2);
kcpsht1=[a,b,d];
cpsht1=(cpk+cpn)*1e3/2; %уд. тепл. шамота
Ssht1=(2*(hsh+hvsht1)+a1)*(2*(hsh+hvsht1)+a2);
Vsht1=(h1+h2)*(Ssht1-Ssh)+Ssht1*hvsht1;
rosht1=1e3;
msht1=(Vsht1*rosht1);
Qshtv1=msht1*cpsht1*abs(Tsht1-Tvnpo);
d1=abs(h-htp)/2;
hnsht1=14e1*1e-3;
Qshtn1=RasQNagVs(Tnver,Tkver,Tkver,h,kktpsht1,d1,a1,kcpsht1,hnsht1,rosht1,ptpn,htp,kktpver,Tvnpo)
Qti=Qsh+Qshtv1+Qvo+Qver+Qshtn1;
w=Qti;
end

function [ vf ] = oprviewf(n,d,ne,a,r,b,h)
c=(a-n*d)/n;
c2=c/2;
for k=1:n
    x(k,1)=c2+(c+d)*(k-1)+r;
    x(k,2)=a-x(k,1);
end
l=a;
arg=[l,r,h,0];
hx=1e-3;
for j=1:ne
    for k=1:ne
        F(j,k)=0;
    end
end
for k=1:n
    F(k,n+1)=RasF21(x(k,1),x(k,2),hx,arg,l,r,a);
    F(k,n+2)=1-F(k,n+1);
end
for k=1:n
    F(k,k+1)=viewfactor38([d,b,c]);
end
vf=F;
end

function [ plo ] = oprPlo(a,b,h,h1,r,n,l)
A(n+1)=a*b;
A(n+2)=a*b+2*(b+a)*(h+h1);
for k=1:n
A(k)=2*pi*r*l;
end
plo=A;
end

function [ e ] = opreps(T,n,nv)
for k=1:n
eps(k)=1;
end
%eps(nv)=epsver(T);
eps(nv)=0.901173740694166;
e=eps;
end

function [ t ] = oprtem(n,T,dT,ne)
for k=1:ne
tem(k)=T;
end
for k=1:n
tem(k)=tem(k)+dT;
end
t=tem;
end

function w = RasQNagVs(Ta0,Tb0,Tc1,h,kktp,d1,h1,kcp,hti,rosht1,ptpn,htp,kktpver,Tvnesh)
k=0;
ra=1e2;
ep=1e-3;
kit=2e2;
d1a=-abs(h-htp);
d1b=abs(h-htp)+2*d1;
Ta=Ta0;
kabh=(Tb0-Ta0)/htp;
while ((ra>ep) && (k<kit))
d1c=(d1a+d1b)/2; 
Tba=Tb0+d1a*kabh;
lama=polyval(kktp,(Ta+Tba)/2);
lama=lama-polyval(kktpver,(Ta+Tba)/2);
lama=1/(2*lama);
dTa=(Tba-Ta);
d2a=htp+d1a;
koa=log(d2a/abs(d1a));
koa=h1/koa;
Qra=ptpn*h1*pi*abs(d1a)*2;
Tbb=Tb0+d1b*kabh;
lamb=polyval(kktp,(Ta+Tbb)/2);
lamb=lamb-polyval(kktpver,(Ta+Tbb)/2);
lamb=1/(2*lamb);
dTb=(Tbb-Ta);
d2b=htp+d1b;
kob=log(d2b/abs(d1b));
kob=h1/kob;
Qrb=ptpn*h1*pi*abs(d1b)*2;
Tbc=Tb0+d1c*kabh;
lamc=polyval(kktp,(Ta+Tbc)/2);
lamc=lamc-polyval(kktpver,(Ta+Tbc)/2);
lamc=1/(2*lamc);
dTc=(Tbc-Ta);
d2c=htp+d1c;
koc=log(d2c/abs(d1c));
koc=h1/koc;
Qrc=ptpn*h1*pi*abs(d1c)*2;
fc=pi*dTc*koc/lamc;
fc=(fc-Qrc)/4;
fa=pi*dTa*koa/lama;
fa=(fa-Qra)/4;
fb=pi*dTb*kob/lamb;
fb=(fb-Qrb)/4;
if (((fa*fc)<0) && ((fb*fc)>0))
    d1b=d1c; 
end
if (((fa*fc)>0) && ((fb*fc)<0))
    d1a=d1c; 
end
k=k+1;
ra=abs(d1a-d1b);
end
d1=d1c;
ht=1e-3; 
dT=1e2; 
d2=d1; 
te0=273.15;
a=kcp(1); b=kcp(2); d=kcp(3);
Tn=23+te0;
Qo=0;
Qr=ptpn*h1*pi*d1/2;
while (d2<hti)
d2=d2+ht; 
Ta=Ta0-dT;
Tb=Tb0+dT;
k=0;
ra=1e2;
ko=log(d2/d1);
ko=h1/ko;
while ((ra>ep) && (k<kit))
Tc=(Ta+Tb)/2; 
lamc=polyval(kktp,(Tc+Tc1)/2);
lamc=1/(2*lamc);
fc=pi*(Tc1-Tc)*ko/lamc;
fc=fc-Qr;
lama=polyval(kktp,(Ta+Tc1)/2);
lama=1/(2*lama);
fa=pi*(Tc1-Ta)*ko/lama;
fa=fa-Qr;
lamb=polyval(kktp,(Ta+Tc1)/2);
lamb=1/(2*lamb);
fb=pi*(Tc1-Tb)*ko/lamb;
fb=fb-Qr;
if (((fa*fc)<0) && ((fb*fc)>0))
    Tb=Tc; 
end
if (((fa*fc)>0) && ((fb*fc)<0))
    Ta=Tc; 
end
ra=abs(Ta-Tb);
k=k+1;
end
if (Tc<Tn)
    Tc=Tn;
end
cp=a+Tc*b+d/(Tc^2); cp=1e3*cp;
if (d2<h)
dm=pi*d2/2;
else phi=asin(h/d2);
dm=d2*phi;
end
dm=dm*h1*ht*rosht1;
dQsht1=cp*dm*(Tc-Tn);
Qo=Qo+dQsht1;
end
lam=polyval(kktp,(Ta0+Tb0)/2);
lam=lam-polyval(kktpver,(Ta0+Tb0)/2);
qiz=lam*abs(Ta0-Tb0)/htp;
lam=polyval(kktp,(Ta0+Tb0)/2);
Tva=Ta0-qiz*hti/lam;
Tvb=Tb0-qiz*hti/lam;
if (Tva<Tvnesh)
    Tva=Tvnesh;
end
if (Tvb<Tvnesh)
    Tvb=Tvnesh;
end
Tsa=(Tva+Tvb)/2;
Tsb=(Ta0+Tb0)/2;
cpa=a+Tsa*b+d/(Tsa^2); cpa=1e3*cpa;
cpb=a+Tsb*b+d/(Tsb^2); cpb=1e3*cpb;
cp=(cpa+cpb)/2;
mti=(((hti+h1)^2-h1^2)*h*rosht1);
dT=abs(Tsa-Tsb);
Qo=mti*dT*cp;
w=Qo;
end