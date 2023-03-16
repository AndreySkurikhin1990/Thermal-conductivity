function dul = DulnevKoefTeplPorZerMat()
m2=0.6635; %2-0.7 mm
v=1; d=15e-6; N=20; te0=273.15; %zpo=0.17;  por=por+zpo; m22=(0.3+0.45)/2;
dan=vyborTempTepl(v,d,N);
la_e=dan(3);
ts=dan(2);
eps=0.899851588244671;
lavo=opredTeploprovVozd(ts-te0)
dul=DulnevPoris207(m2,T,eps,lavo,lam1)
end

function lam4 = DulnevPoris207(m2,T,eps,lavo,lam1)
r=(2+0.7)/2/2; 
d=2*r;
Nk=sqrt(m2^2-10*m2+9);
Nk=(m2+3+Nk)/2/m2;
y2=(3.3e-3)*((1-m2)^(-2/9));
%pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2);
%y2=y2*(pud+9.8*ro1*(1-m2)*hsl)^(1/3);
%y1=1e-20; y2=y1;
eta=1e-4; %eta=1e-5...1e-2
y1=10e-4; 
r1=y1*r
y2=y1/sqrt(eta); %y1=(10...50)*1e-4; Hm=1...2
y3=2*sqrt(Nk-1)/Nk;
r3=y3*r
y4=y3/((1-m2)^(1/3));
%hsh=2e-3/2;
hsh=(0+2e-3)/2;
po=RazPor(r,r1,hsh,r3);
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

function rp = RazPor(r,r1,hsh,r3)
v1=pi*(r1^2)*hsh;
ha=r-sqrt(r^2-r3^2);
v2=pi*(ha-hsh)*(r3^2);
va=pi*(r-ha/3)*(ha^2);
vb=pi*(hsh^2)*(r-hsh/3);
vshse=abs(va-vb);
ha=ha-hsh;
rp=sqrt((v1+v2+vshse)/pi/ha);
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

function [ tem ] = RasTempZer(Temp,d,N)
c=Temp(1); h=30e-3; M=0;
M(1,1)=h^2; M(1,2)=h; M(2,1)=(h/2)^2; M(2,2)=h/2; b=0;
b(1)=Temp(3)-c; b(2)=Temp(2)-c;
y=0; y=inv(M)*b'; a=y(1); b=y(2);
xn=h/2-N*d; xk=h/2+N*d;
Tn=a*(xn^2)+b*xn+c; Tk=a*(xk^2)+b*xk+c;
tem = [Tn   Tk  a   b   c];
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
la = [q tsr lam ko(1) ko(2) ko(3) da(1) da(2) da(3) da(4) da(5)];
end