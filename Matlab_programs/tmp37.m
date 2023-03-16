function [ tm ] = tmp37()
format long g;
%определяет коэффициент рассеяния для таблетки KBr+вермикулит
t=2*PoisKorn();
%nom=1 - для цилиндров с эллипсами в основании, 2 - для прямоугольных параллелепипедов
for nom=1:2
    for no=1:4
switch (no)
    case 1
ti_60=RasFra60(t,nom);
ti_60=ti_60'
    case 2
ti_60_100=RasFra60_100(t,nom);
ti_60_100=ti_60_100'
    case 3
ti_100_150=RasFra100_150(t,nom);
ti_100_150=ti_100_150'
    case 4
ti_150_200=RasFra150_200(t,nom);
ti_150_200=ti_150_200'
end
    end
end
tm=[0];
end
function [ ras60 ] = RasFra60(di,ide)
mkbr=[250.239 249.740 249.223 250.395 250.336 249.55]; 
mv=[0.283 0.464 0.812 0.22 0.547 0.777]; 
tol=[0.73 0.72 0.72 0.72 0.71 0.7];
rokbr=2.75; rov=0.49; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %в см3
vv(k)=mv(k)/(1e3*rov); %в см3
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %в мкм
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
mo=6e1/2; %МО
si=6e1/di; %СКО
disp('Фракция менее 60 мкм');
Tr=PloshSeche(length(vv),vv,mo,si,0,60,ide);
ras60=alphaRas(tol,length(tol),Tr);
end
function p = prov(hk,mi,ma)
hp=hk;
if (hk<mi)
    hp=mi;
end
    if (hk>ma)
        hp=ma;
    end
    p=hp;
end
function [ al ] = alphaRas(tol,n,Tr)
alp=0;
for k=1:n
    alp(k)=-log(Tr(k))/tol(k);
end
al=alp;
end
function pp = PloshSeche(n,vv,mo,si,mi,ma,no)
Nr=1e3; secn=0;
for k=1:n
    secn(k)=0;
end
for g=1:Nr
for j=1:n
sv=vv(j); suv=0; 
k=1; nk=1e9; sech=0;
while (suv<sv)
hk=mo+1.*si*randn();
hk=prov(hk,mi,ma);
ak=mo+1.*si*randn();
ak=prov(ak,mi,ma);
bk=mo+1.*si*randn();
bk=prov(bk,mi,ma);
%rc(k)=13e3*rand()/2; %13 мм - диаметр таблетки
%phic(k)=2*pi*rand(); %угол phi центра цилиндра
%zc(k)=tol(j)*rand(); %положение на Oz
tetak=3*pi/2+pi*rand(); %угол наклона плоскости основания к оси Z
apk=ak*abs(cos(tetak));
nuk=pi/2-tetak; %угол наклона плоскости, перпендикулярной основанию
hpk=hk*cos(nuk);
if (no==1)
sech=sech+(pi*apk/2+hpk)*bk/2;
suv=suv+pi*hk*ak*bk/4;
else 
    sech=sech+(apk+hpk)*bk;
    suv=suv+hk*ak*bk;
end
k=k+1;
if (k>nk)
    break; 
end
end
secn(j)=secn(j)+sech;
end
end
secn=secn/Nr/1e6;
pt=pi*13^2/4;
pp=(pt-secn)/pt; %пропускание за счет рассеяния
end
function [ ras60_100 ] = RasFra60_100(di,ide)
mkbr=[250 249.629 249.294 249.706 249.510 249.307 250.328 249.604 249.206]; 
mv=[0.255 0.539 0.809 0.295 0.517 0.756 0.36 0.534 0.843]; 
tol=[0.72 0.71 0.7 0.7 0.73 0.72 0.74 0.7 0.76]; 
rokbr=2.75; rov=0.52; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
mo=(1e2+6e1)/2; %МО
si=(1e2-6e1)/di; %СКО
disp('Фракция 60-100 мкм');
Tr=PloshSeche(length(vv),vv,mo,si,6e1,1e2,ide);
ras60_100=alphaRas(tol,length(tol),Tr);
end
function [ ras100_150 ] = RasFra100_150(di,ide)
mkbr=[249.913 249.607 249.218 249.929 249.695 249.306 250.405 249.625 249.348];
mv=[0.315 0.473 0.709 0.293 0.528 0.83 0.27 0.493 0.764]; 
tol=[0.74 0.74 0.72 0.72 0.71 0.7 0.78 0.73 0.76]; 
rokbr=2.75; rov=0.53; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
mo=(1e2+15e1)/2; %МО
si=abs(15e1-1e2)/di; %СКО
disp('Фракция 100-150 мкм');
Tr=PloshSeche(length(vv),vv,mo,si,1e2,15e1,ide);
ras100_150=alphaRas(tol,length(tol),Tr);
end
function [ ras150_200 ] = RasFra150_200(di,ide)
mkbr=[250.882 249.590 249.213 250.299 249.441 249.365];
mv=[0.320 0.533 0.849 0.223 0.502 0.797]; 
tol=[0.76 0.72 0.69 0.73 0.73 0.73]; 
rokbr=2.75; rov=0.56; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
mo=(2e2+15e1)/2; %МО
si=(2e2-15e1)/di; %СКО
disp('Фракция 150-200 мкм');
Tr=PloshSeche(length(vv),vv,mo,si,15e1,2e2,ide);
ras150_200=alphaRas(tol,length(tol),Tr);
end
function pk = PoisKorn()
a=1e-5; b=1e5; ep=1e-6; Nit=1e3;
k=0; ffv=1e-4; ra=abs(a-b); %ffv = 0,01 % - такая доля частиц не попадает в диапазон фракции
while (ra>ep)
    c=(a+b)/2;
    fa=erfc(a)-ffv;
    fb=erfc(b)-ffv;
    fc=erfc(c)-ffv;
    if ((fc*fb)>0) 
        if ((fa*fc)<0) 
            b=c; 
        end
    end
    if ((fc*fa)>0) 
        if ((fb*fc)<0) 
            a=c; 
        end
    end
    k=k+1;
    if (k>Nit) 
        break; 
    end
    ra=abs(a-b);
end
pk=c;
end