function [ ktp ] = DulnevZernkvi(por, te, srk, na, laefm, lavo, wsio, walo, wmgo, srch, n, tena, dtosc, dm, kuscv, tkuscv, stchv, x0)
up=1e0; koscve=1e0; uo=urovPod(por);
for k=1:n
ts=te(k); %koscve=opredKTPTKTochSha(kuscv, tkuscv, ts, dm);
stchve=opredKTPTKTochSha(stchv, te, ts, n);
epsil=stchve*koscve; 
stch(k)=epsil; 
la_voz=opredKTPTKTochSha(lavo, te, ts, n);
la_e=opredKTPTKTochSha(laefm, te, ts, n);
lam=opredDulnLam1kvi(por, ts, epsil, la_voz, la_e,x0); 
lam=urovOtsechen(lam, la_e, uo);
lam=urovPodder(lam, la_e, up);
srk(k)=lam;
end
srk = proverkakvi(por,srk,laefm,lavo,up,uo);
ktp=srk;
end

function ol1 = opredDulnLam1kvi(po,T,eps,lavo,lae,r)
ladb=1e5; lada=1e-5; ra=1e2; ep=1e-7; h=0; la=0; kit=1e3;
while ((ra>ep) && (h<kit))
ladc=(lada+ladb)/2;    
fa=DulnevSloForkvi(po,T,eps,lavo,lada,r);
fa=fa-lae;
fb=DulnevSloForkvi(po,T,eps,lavo,ladb,r);
fb=fb-lae;
fc=DulnevSloForkvi(po,T,eps,lavo,ladc,r);
fc=fc-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
ol1=ladc;
end

function lam4 = DulnevSloForkvi(m2,T,eps,lavo,lam1,r)
d=2e0*r;
Nk=sqrt(m2^2-10e0*m2+9e0);
Nk=(m2+3e0+Nk)/2/m2;
%y2=3.3e-3*(1-m2)^(-2/9); %pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2); %y2=y2*(pud+9.8*ro1*(1-m2)*hsl)^(1/3); %y1=1e-20; y2=y1;
eta=1e-4; 
y1=10e-4; y2=y1/sqrt(eta); %y1=(10+50)*1e-4/2;
y3=2*sqrt(Nk-1)/Nk;
y4=y3/((1-m2)^(1/3));
hsh=2e-3/2;
%hsh=1e-10; %hshsr=hsh/r;
epf=1e-2;
fbn=opredFiBolN1(y2,y2/y3);
fbnn=opredFiBolNN1(lavo/lam1,m2);
if (fbn>1) 
    fb=fbnn; 
if (fbnn>1)
    fb=0;
end
elseif (fbnn>1) 
    fb=fbn; 
elseif (abs(fbn-fbnn)<epf) 
    fb=(fbn+fbnn)/2; 
elseif (fbn>fbnn) 
    fb=fbnn; 
else fb=fbn;
end
%fbs=opredFiBol(y2/y3,y2); fbn=opredFiBolN(y2,y2/y3); %if (abs(fbn-fbs)>1e-1) fb=(fbs+fbn)/2; else fb=fbs; end
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

function r = ProvAdek(m)
if (m>1e0) 
    m=1e0; 
end
if (m<0) 
    m=0;
end
r=m;
end

function ktpo = opredKTPTKTochSha(ktptks, te, temp, n)
f = 1; p = 0; ep=1e-4; ktp = 0;
if ((temp>te(1)) && (temp<te(n)))
for k = 1:n
if ((te(k) >= temp) && (f>0))
        p = k; f = 0; break;
end
end
elseif (temp>=te(n))
    p=n; f=0;
elseif (temp<=te(1))
    p=2; f=0;
end
if ((f==0) && (p>1))
    x2=te(p);
    x1=te(p-1);
	dt = x2 - x1;
if (abs(dt) > ep)
    y2=ktptks(p);
    y1=ktptks(p - 1);
    b=y1;
			ko = (y2 - y1) / dt;
            if (p==n)
                b=y2;
            end
			ktp = b + ko*(temp - x1);
else 
    ktp=0;
end
end
ktpo=ktp;
end

function t = opredFiBolN1(x, y)
y1=[0,47e-3,113e-3,195e-3,3e-1,413e-3,563e-3,725e-3,832e-3,932e-3,1e0];
y01=[0,261e-3,411e-3,526e-3,618e-3,7e-1,774e-3,837e-3,9e-1,958e-3,1e0];
n=length(y1);
fib(1)=0; 
for k=2:n 
    fib(k)=fib(k-1)+1e-1; 
end
f1=1; f2=1; q1=0; q2=0;
for k=1:n
if ((f1>0) && (y<y01(k))) 
        q1=k; f1=0; break;
end
if ((f2>0) && (y<y1(k))) 
        q2=k; f2=0; break;
end
end
if (q1==0)
    q1=2; 
end
if (q2==0)
q2=2;
end
fb01=fib(q1-1)+(fib(q1)-fib(q1-1))*(y-y01(q1-1))/(y01(q1)-y01(q1-1));
fb1=fib(q2-1)+(fib(q2)-fib(q2-1))*(y-y1(q2-1))/(y1(q2)-y1(q2-1));
fibo=fb1+(fb01-fb1)*(x-1e0)/(1e-1-1e0);
t=fibo;
end

function t = opredFiBolNN1(nu, m2)
cb=1e2; ca=-1e1; ra=1e2; ep=1e-6;
h=0; kit=100;
while ((ra>ep) && (h<kit))
cc=(ca+cb)/2;
fa=2*(ca^3)-3*(ca^2)+1-m2;
fb=2*(cb^3)-3*(cb^2)+1-m2;
fc=2*(cc^3)-3*(cc^2)+1-m2;
if ((fc*fb>0) && (fa*fc<0)) 
    cb=cc;
end
if ((fc*fa>0) && (fb*fc<0)) 
    ca=cc;
end
ra=abs(fa-fb); h=h+1;
end
fi1=abs(cc/(2-cc)); 
fi1=sqrt(fi1);
t=(nu^2)*(7.5-11*nu+4.5*(nu^2))*(1-fi1)+fi1;
end

function  fb = opredFiBol(bb,rr)
nb=1e2;
sn=0;
for k=0:nb
    n=2*k+1;
    I1xy=besseli(1,n*pi*bb*rr);
    I1x=besseli(1,n*pi*rr);
    K1xy=besselk(1,n*pi*bb*rr);
    K1x=besselk(1,n*pi*rr);
   sn=sn+I1xy*(I1x*K1xy-K1x*I1xy)/(n^2)/I1x; 
end
fb=1-16*sn/(pi^2);
end

function fibo = opredFiBolN(x,y)
y1=[0,0.047,0.113,0.195,0.3,0.413,0.563,0.725,0.832,0.932,1];
y01=[0,0.261,0.411,0.526,0.618,0.7,0.774,0.837,0.9,0.958,1];
fib=0:1e-1:1e0; f1=1; f2=1; q1=1; q2=1;
for k=2:length(y1)
    if ((f1>0) && (y<y01(k)))
        q1=k; f1=0;
    end
    if ((f2>0) && (y<y1(k)))
        q2=k; f2=0;
    end
end
if (q1==1)
    q1=2;
end
if (q2==1)
    q2=2;
end
fb01=fib(q1-1)+(fib(q1)-fib(q1-1))*(y-y01(q1-1))/(y01(q1)-y01(q1-1));
fb1=fib(q2-1)+(fib(q2)-fib(q2-1))*(y-y1(q2-1))/(y1(q2)-y1(q2-1));
fibo=fb1+(fb01-fb1)*(x-1e0)/(1e-1-1e0);
end

function fibo = opredFiBolNN(nu,m2)
cb=1e2; ca=-1e1; ra=1e2; ep=1e-9; h=0; la=0; kit=1e3; m2=m2';
while ((ra>ep) && (h<kit))
cc=(ca+cb)/2;
fa=2*(ca^3)-3*(ca^2)+1-m2;
fb=2*(cb^3)-3*(cb^2)+1-m2;
fc=2*(cc^3)-3*(cc^2)+1-m2;
la = VybCAB(fa,fb,fc,ca,cb,cc);
ca=la(1); cb=la(2); cc=la(3);
ra=abs(fa-fb); h=h+1;
end
cc=cc';
fi1=sqrt(cc/(2-cc)); fi1=abs(fi1);
fibo=(nu^2)*(7.5-11*nu+4.5*(nu^2))*(1-fi1)+fi1;
end

function [ vy ] = VybCAB(fa,fb,fc,xa,xb,xc)
if (((fa*fc)<0) && ((fb*fc)>0))
    xb=xc; 
end
if (((fa*fc)>0) && ((fb*fc)<0) )
    xa=xc;
end
vy = [xa,xb,xc];
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
a=opredKoefAkkomodLandau(T);
%a=0.9;
cz=8.42e-3;
Ty=113;
la0=cz/H0/(1+Ty/T);
Kn=la0*H0/H/de;
B=4*gam/(gam+1)*(2-a)/a*Kn/Pr;
la=lamg/(1+B);
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
ro=opredPlotnVozd(T);
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

function u = urovPod(po)
porex=95e-2; pormakvi=53e-2; 
porist=28e-2; porg=5e-1; 
if (po>6e-1) 
     m2=porex;
elseif (po>36e-2)
    m2=pormakvi;
else
    m2=po;
end
 if (po<porist)
     uo=opredUrovPodderM03(po); 
elseif (po<porg) 
    uo=1e0/(1e0-po); 
else
     uo=1/(1-m2);
 end
u=(uo^3);
end

function fl = urovOtsechen(pro,ref,urot)
if (pro<0)
    f=0;
elseif (pro>(urot*ref))
    f=0;
else f=1;
end
if (f>0)
fl=pro;
else
    fl=0;
end
end

function t = opredUrovPodderM03(por)
n=3; 
pam(3)=-8.93227577791347; 
pam(2)=8.89444783404518; 
pam(1)=1.06833435021354;
p=0; s=0; 
for k=1:n 
    s=s+pam(k)*(por^p); 
    p=p+1;
end
t=1/(1-s*por);
end

function fl = urovPodder(pro,ref,urpo)
f=1; 
if (pro<0) 
    f=0; 
elseif (pro<(urpo*ref)) 
    f=0; 
else
    f=1;
end
if (f>0) 
    fl=pro;
else
    fl=0;
end
end

function [ pr ] = proverkakvi(po,ta,laefm,lavo,up,uo)
for k=1:length(ta)
	ta(k)=urovOtsechen(ta(k), laefm(k), uo); % не может быть выше Ё “ѕ*10 с учетом пористости
	ta(k)=urovPodder(ta(k), lavo(k), up); %не может быть ниже  “ѕ воздуха
end
pr=ta;
end