function [ t ] = VasFrayObsh(laefm,lavo, po, n)
nit=1000; h=0; L=1e0; 
maxlam=1e5; lamb=maxlam; minlam=1e-9; lama=minlam; 
tocras=1e-9; ep=tocras; ra=1e0; 
x=Poiskkx((1e0-po)/4e0,tocras,nit); 
uo=urovPod(po); up=1;
db=L*x; hlm=x/(0.5-x); 
hlb=hlm/(1e0+hlm);
mu=1e-2; E=8.428e6;
Pn=125.0*9.8*30e-3; kk=(1.5+2)/2; kb=(2.2+2.9)/2;
km=(4+5)/2; nu=2*(1-mu^2)/E; kc=(35e-2+45e-2)/2e0;
rp=(725e-3)*((nu*Pn*L/2e0)^(1/3)); Q=(74e-2/(1-po))^(1/3);
hsh=L*km*1e-3; sfk=pi*(rp^2);
for k=1:n
	laef=laefm(k); lavoz=lavo(k); h=0;
while ((ra>ep) && (h<nit))
lamc=(lama+lamb)/2e0;
lamz=lavoz; nugza=lamz/lama; nugzb=lamz/lamb; nugzc=lamz/lamc;
if (Pn<3e5) 
    lampna=lama*(Pn^(2/3))*kc/75/Q; 
    lampnb=lamb*(Pn^(2/3))*kc/75/Q;
    lampnc=lamc*(Pn^(2/3))*kc/75/Q;
else
    lampna=lama*((Pn/E)^(4/9))*kb/Q;
    lampnb=lamb*((Pn/E)^(4/9))*kb/Q;
    lampnc=lamc*((Pn/E)^(4/9))*kb/Q;
end
lamksha=2*sqrt(2)*sfk*lama/L/Q/hsh/kk; 
lamkshb=2*sqrt(2)*sfk*lamb/L/Q/hsh/kk; 
lamkshc=2*sqrt(2)*sfk*lamc/L/Q/hsh/kk;
lamka=lamksha+lampna;
lamkb=lamkshb+lampnb;
lamkc=lamkshc+lampnc;
Aa=(hlb^2)*1e3*nugza/4/kk/km+lamka/lama;
Ab=(hlb^2)*1e3*nugzb/4/kk/km+lamkb/lamb;
Ac=(hlb^2)*1e3*nugzc/4/kk/km+lamkc/lamc;
nuga=nugza; nugb=nugzb; nugc=nugzc;
fa=1/(hlb^(-2)+Aa)+nuga*((1-hlb)^2)+2/(1+hlm+1/(nuga*hlb))-laef/lama;
fb=1/(hlb^(-2)+Ab)+nugb*((1-hlb)^2)+2/(1+hlm+1/(nugb*hlb))-laef-lamb;
fc=1/(hlb^(-2)+Ac)+nugc*((1-hlb)^2)+2/(1+hlm+1/(nugc*hlb))-laef/lamc;
if ((fc*fb>0) && (fa*fc<0)) 
    lamb=lamc; 
elseif ((fc*fa>0) && (fb*fc<0)) 
    lama=lamc; 
end
h=h+1; 
ra=abs(lama-lamb); 
end
lto(k)=lamc;
end
lto=proverkakvi(po,lto,laefm,lavo,up,uo);
t=lto;
end

function t = Poiskkx(kx,ep,nit)
xa=0; xb=5e-1; ra=abs(xa-xb); h=0;
while ((ra>ep) && (h<nit))
xc=(xa+xb)/2;
fa=4*(xa^3)-3*(xa^2)+kx; 
fb=4*(xb^3)-3*(xb^2)+kx; 
fc=4*(xc^3)-3*(xc^2)+kx; 
if ((fc*fb>0) && (fa*fc<0)) 
    xb=xc; 
elseif ((fc*fa>0) && (fb*fc<0)) 
    xa=xc; 
end
ra=abs(xa-xb); h=h+1; 
end
t=xc; 
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
u=uo;
end

function [ pr ] = proverkakvi(po,ta,laefm,lavo,up,uo)
for k=1:length(ta)
	ta(k)=urovOtsechen(ta(k), laefm(k), uo); % не может быть выше Ё “ѕ*10 с учетом пористости
	ta(k)=urovPodder(ta(k), lavo(k), up); %не может быть ниже  “ѕ воздуха
end
pr=ta;
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

function fl = urovOtsechen(pro,ref,urot)
if (pro<0)
    f=0;
elseif (pro>(urot*ref))
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