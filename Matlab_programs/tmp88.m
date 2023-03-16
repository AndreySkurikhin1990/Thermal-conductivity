%функция находит показатель преломления по заданным коэффициентам отражения
function t = tmp88()
format long g;
k=0;
dl=vyborVesch(k,0); %0 - kremnezyom, 1 - shamot
sc=vyborVesch(k,1);
mdk=vyborVesch(k,2);
for k=1:length(sc)
    k
%np(k)=opredPokPrel(sc(k));
%reflOp(k)=Refl_Shamot(np(k)); %опорное значение
reflOp(k)=1-sc(k);
np(k)=opredPokPrel2(reflOp(k)); %Поиск n2
end
reflOp
np
t=0;
end

function [ vy ] = vyborVesch(v,m)
if (v==0) %Silika - kremnezyom
dv=[1,1.45,1.79,2,2.51,2.88,3,3.33,3.47,3.66,3.96,4,4.29,4.63,5,5.45,6,6.41,7,7.67,8,8.83,9,9.64,10,11,12,13,14,15];
stch=[0.44,0.4,0.347,0.333,0.311,0.384,0.367,0.4,0.5,0.6,0.7,0.721,0.8,0.9,0.965,0.988,0.97,0.958,0.974,0.903,0.971,0.929,0.907,0.936,0.953,0.968,0.974,0.805,0.966,0.814];
sod=[93,7,0,0,0,0]; %SiO2, Al2O3, MgO, Fe2O3, FeO, CaO
elseif (v==1) %Shamot
dv=[1,1.24,2,2.71,2.86,3,3.16,3.77,4,4.09,4.47,4.85,5.61,6,6.19,7,7.83,8,9,10,11,12,13,14,14.59,15];    
stch=[0.468,0.472,0.467,0.5,0.517,0.513,0.5,0.6,0.68,0.7,0.8,0.9,0.968,0.973,0.978,0.975,0.991,0.986,0.926,0.955,0.968,0.957,0.955,0.952,0.9,0.845];
sod=[56, 39.6, 0.7, 2.3, 0.7, 0.7];
end
if (m==0)
    vy=dv;
elseif (m==1)
    vy=stch;
elseif (m==2)
    vy=sod;
end
end

function pr = poiskOtrSposFrenel(n2)
koot=0; phin=0; phik=pi/2; hphi=1e-4; n1=1;
%Усреднение коэффициента отражения по углам
Nphi=(phik-phin)/hphi;
for k=1:Nphi
    phi=phin+(phik-phin)*k/Nphi;
    phis=PoiskPhiShtr(phi,n1,n2);
    ko=PoiskReflPerv(phi,phis,n1,n2);
    koot=koot+ko*hphi;
end
pr=koot/abs(phik-phin);
end

%в плоскости ksis, etas
function pr10 = PraCha10(phi,phis,n1,n2)
pr10=n1*sin(phi)-n2*sin(phis);
end

function pu = provUgla(phi)
if (phi<(-pi))
    phi=phi+pi;
end
if (phi<(-pi/2))
    phi=-phi-pi;
end
if (phi>pi)
    phi=phi-pi;
end
if (phi>(pi/2))
    phi=pi-phi;
end
if (phi<0)
    phi=abs(phi);
end
pu=phi;
end

function [ pesdp ] = PoiskPhiShtr(phieta,n1,n2)
a=0; b=pi/2; ep=1e-6; Nit=1e4; k=0; ra=abs(a-b); %находим угол преломления
while ((ra>ep) && (k<Nit))
    c=(a+b)/2;
    fa=PraCha10(phieta,a,n1,n2);
    fb=PraCha10(phieta,b,n1,n2);
    fc=PraCha10(phieta,c,n1,n2);
    ro=VybCAB(fa,fb,fc,a,b,c); 
    a=ro(1); b=ro(2); c=ro(3); 
    k=k+1; ra=abs(a-b);
end
pesdp=provUgla(c);
end

%поиск изменения угла после преломления
function rop = PoiskReflPerv(phieta,phietas,n1,n2)
rpa=(n2*cos(phieta)-n1*cos(phietas))/(n2*cos(phieta)+n1*cos(phietas)); %коэффициент отражения параллельный плоскости падения
rpe=(n1*cos(phieta)-n2*cos(phietas))/(n1*cos(phieta)+n2*cos(phietas)); %коэффициент отражения перпендикулярный к плоскости падения
rop=(abs(rpa^2)+abs(rpe^2))/2;
end

function pp = opredPokPrel(eps)
nb=1e3; na=1+1e-6; ra=1e2; ep=1e-6; h=0; kit=1e4;
while ((ra>ep) && (h<kit))
nc=(na+nb)/2;
fa=emdiel(na)-eps;
fb=emdiel(nb)-eps;
fc=emdiel(nc)-eps;
la = VybCAB(fa,fb,fc,na,nb,nc);
na=la(1); nb=la(2); nc=la(3);
ra=abs(fa-fb); h=h+1;
end
pp=nc;
end

function pp = opredPokPrel2(roop)
ep=1e-5; nb=1e2; na=1+ep; ra=1e2; h=0; kit=1e4;
while ((ra>ep) && (h<kit))
nc=(na+nb)/2;
fa=poiskOtrSposFrenel(na)-roop;
fb=poiskOtrSposFrenel(nb)-roop;
fc=poiskOtrSposFrenel(nc)-roop;
la = VybCAB(fa,fb,fc,na,nb,nc);
na=la(1); nb=la(2); nc=la(3);
ra=abs(fa-fb); h=h+1;
end
pp=nc;
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

function kootr = Refl_Shamot(np)
    np=abs(np);
    ronu=(5e-1)+((np-1)*(3*np+1))/(6*(np+1)^2);
    ronu=ronu-(2*(np^3)*(np^2+2*np-1))/((np^2+1)*(np^4-1)); 
    podln=(np-1)/(np+1);
    podln=abs(podln);
    ronu=ronu+(8*(np^4)*((np^4)+1)*log(np))/((np^2+1)*((np^4-1)^2));
    ronu=ronu+(np^2)*((np^2-1)^2)*log(podln)/((np^2+1)^3);
    ronu=abs(ronu);
    if (ronu>1)
        ronu=1;
    end
kootr=ronu;
end

function [ dlv ] = dliny_voln()
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end