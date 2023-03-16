function t = tmp58()
format longe;
di=1e-4;
di=2*erfcinv(di);
n=4;
mini(1)=0; maxi(1)=6e1; maxi(2)=1e2; maxi(3)=15e1; maxi(4)=2e2;
for k=1:n
    maxi(k)=1e-3*maxi(k); %в миллиметрах
end
for k=2:n
    mini(k)=maxi(k-1);
end
for k=1:n
obsr(k)=0; obpr(k)=0; nor(k)=0;
mx(k)=(maxi(k)+mini(k))/2; %МО
sig(k)=abs((mini(k)-maxi(k))/di); %СКО
y1(k)=(mini(k)-mx(k))/sig(k);
y2(k)=(maxi(k)-mx(k))/sig(k);
obsr(k)=SrObShary(mx(k),sig(k),y1(k),y2(k));
nor(k)=provFuncNorm(mini(k),maxi(k),sig(k),mx(k));
nx2(k)=rasnxx2(mx(k)/1e3,sig(k)/1e3,mini(k)/1e3,maxi(k)/1e3); %средняя площадь, м
src(k)=SrRazCh(mx(k),sig(k),mini(k),maxi(k)); %средний размер частиц, мм
%obpr(k)=SrObPPP(y1(k),y2(k));
end
sig=sig';
src=(src')*1e-3;
nx2=nx2';
b=1.8923;
lasr=(1.33+27)*(1e-6)/2;
for nom=1:n
    n=0; tol=0; vv=0;
    switch nom
        case 1
n=rasKonVer60(3);
tol=rasKonVer60(2);
vv=rasKonVer60(1);
disp('Фракция менее 60 мкм');
        case 2
n=rasKonVer60_100(3);
tol=rasKonVer60_100(2);
vv=rasKonVer60_100(1);          
disp('Фракция от 60 до 100 мкм');
        case 3
n=rasKonVer100_150(3);
tol=rasKonVer100_150(2);
vv=rasKonVer100_150(1);
disp('Фракция от 100 до 150 мкм');
        case 4
n=rasKonVer150_200(3);
tol=rasKonVer150_200(2);
vv=rasKonVer150_200(1);
disp('Фракция от 150 до 200 мкм');
    end
    ko=0; n0=0; nc=0;
for k=1:n
    vt=pi*((13e-3/2)^2)*tol(k); %объем таблетки в мкм3
    nc(k)=vv(k)/(obsr(nom)*((1e-3)^3)); %число частиц в таблетке
    n0(k)=nc(k)/vt; %число частиц в 1 м3
    ko(k)=n0(k)*(pi*nx2(nom)/2+b*lasr*src(nom)/2);
end
ko=ko';
n0=n0';
end
obsr=(1e3^3)*obsr'; %в микронах
nor=nor';
%obpr=(1e3^3)*obpr'
t=0;
end

function t = SrObShary(mx,sig,y1,y2)
ko=2*sqrt(2*pi)/(2^3); 
t2=-exp(-(y2^2)/2); t2=t2*((y2^2)+2);
t1=-exp(-(y1^2)/2); t1=t1*((y1^2)+2);
fv=(ko/3)*(sig^3)*(t2-t1);
t2=-exp(-(y2^2)/2); t1=-exp(-(y1^2)/2);
fv=fv+ko*sig*(mx^2)*(t2-t1);
z1=y1/sqrt(2); z2=y2/sqrt(2); ko2=2*sqrt(2);
t2=exp(-(z2^2)); t2=-z2*t2/2;
t1=exp(-(z1^2)); t1=-z1*t1/2;
f2=sqrt(pi)*erf(z2)/2; f1=sqrt(pi)*erf(z1)/2;
fv=fv+ko*mx*(sig^2)*ko2*(t2-t1+(f2-f1)/2);
fv=fv+(mx^3)*ko*(f2-f1)*sqrt(2)/3;
t=fv;
end

function t = SrObPPP(y1,y2)
t2=-exp(-(y2^2)/2)*sqrt(pi)/2; t1=-exp(-(y1^2)/2)*sqrt(pi)/2;
t=((t2-t1)^3);
end

function t = provFuncNorm(amin,amax,sig,mx)
z1=(amin-mx)/sig/sqrt(2); z2=(amax-mx)/sig/sqrt(2);
t=(erf(z2)-erf(z1))/2;
end

function [ t ] = rasKonVer60(kl)
mkbr=[250.239,249.740,249.223,250.395,250.336,249.55]; 
mv=[0.283,0.464,0.812,0.22,0.547,0.777]; 
tol=[0.73,0.72,0.72,0.72,0.71,0.7];
rokbr=2.75; rov=0.49; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %в см3
vv(k)=mv(k)/(1e3*rov); %в см3
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %в мкм
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
if (kl==1)
    t=vv/((1e3)^3)/((1e3)^3); %в м3
elseif (kl==2)
    t=tol*1e-3; %в м
elseif (kl==3)
    t=n;
end
end

function [ t ] = rasKonVer60_100(kl)
mkbr=[250,249.629,249.294,249.706,249.510,249.307,250.328,249.604,249.206]; 
mv=[0.255,0.539,0.809,0.295,0.517,0.756,0.36,0.534,0.843]; 
tol=[0.72,0.71,0.7,0.7,0.73,0.72,0.74,0.7,0.76]; 
rokbr=2.75; rov=0.52; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
if (kl==1)
    t=vv/((1e3)^3)/((1e3)^3); %в м3
elseif (kl==2)
    t=tol*1e-3; %в м
elseif (kl==3)
    t=n;
end
end

function [ t ] = rasKonVer100_150(kl)
mkbr=[249.913,249.607,249.218,249.929,249.695,249.306,250.405,249.625,249.348];
mv=[0.315,0.473,0.709,0.293,0.528,0.83,0.27,0.493,0.764]; 
tol=[0.74,0.74,0.72,0.72,0.71,0.7,0.78,0.73,0.76]; 
rokbr=2.75; rov=0.53; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
if (kl==1)
    t=vv/((1e3)^3)/((1e3)^3); %в м3
elseif (kl==2)
    t=tol*1e-3; %в м
elseif (kl==3)
    t=n;
end
end

function [ t ] = rasKonVer150_200(kl)
mkbr=[250.882,249.590,249.213,250.299,249.441,249.365];
mv=[0.320,0.533,0.849,0.223,0.502,0.797]; 
tol=[0.76,0.72,0.69,0.73,0.73,0.73]; 
rokbr=2.75; rov=0.56; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vv(k)=vv(k)*((1e4)^3); 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
end
if (kl==1)
    t=vv/((1e3)^3)/((1e3)^3); %в м3
elseif (kl==2)
    t=tol*1e-3; %в м
elseif (kl==3)
    t=n;
end
end

function t = rasnxx2(mx,sig,amin,amax)
z1=(amin-mx)/sqrt(2)/sig;
z2=(amax-mx)/sqrt(2)/sig;
t1=-(z1/2); t1=t1*exp(-(z1^2));
t2=-(z2/2); t2=t2*exp(-(z2^2));
f1=sqrt(pi)*erf(z1)/2;
f2=sqrt(pi)*erf(z2)/2;
nx2=(t2-t1+f2-f1)*2*(sig^2)/sqrt(pi);
t1=-exp(-(z1^2))/2;
t2=-exp(-(z2^2))/2;
nx2=nx2+2*sqrt(2)*sig*mx*(t2-t1)/sqrt(pi);
nx2=nx2+(mx^2)*(erf(z2)-erf(z1))/2;
t=nx2;
end

function t = SrRazCh(mx,sig,amin,amax)
z1=(amin-mx)/sqrt(2)/sig;
z2=(amax-mx)/sqrt(2)/sig;
t1=-exp(-(z1^2))/2; t2=-exp(-(z2^2))/2;
xs=sig*sqrt(2/pi)*(t2-t1);
xs=xs+mx*(erf(z2)-erf(z1))/2;
t=xs;
end