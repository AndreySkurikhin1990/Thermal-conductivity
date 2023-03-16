function [ alsr ] = SredGrafMod()
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=privedkEdiPropus(TrkoVer5311()); 
Tal2=privedkEdiPropus(TrkoVer5312());
Tal3=privedkEdiPropus(TrkoVer5313()); 
Tal4=privedkEdiPropus(TrkoVer53101()); 
Tal5=privedkEdiPropus(TrkoVer53102()); 
Tal6=privedkEdiPropus(TrkoVer53103()); 
mkbr=0; mkbr(1)=250.239; 
mv=0; mv(1)=0.283; 
tol=0; tol(1)=0.73; 
rokbr=2.75; 
rov=0; rov(1)=0.49;
%---
mkbr(2)=249.740; 
mv(2)=0.464; 
tol(2)=0.72; 
rov(2)=rov(1); 
%---
mkbr(3)=249.223; 
mv(3)=0.812; 
tol(3)=0.72; 
mkbr(4)=250.395; 
mv(4)=0.22; 
tol(4)=0.72; 
rov(3)=0.49; 
rov(4)=rov(3); 
%---
mkbr(5)=250.336; 
mv(5)=0.547; 
tol(5)=0.71; 
rov(5)=0.49; 
%---
mkbr(6)=249.55; 
mv(6)=0.777; 
tol(6)=0.7; 
rov(6)=rov(5); 
%---
md=0; od=0; xv=0; vkbr=0; vv=0; 
for k=1:6
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov(k)); 
    md(k)=mv(k)/(mv(k)+mkbr(k));
    od(k)=(vv(k)/(vv(k)+vkbr(k)));
    xv(k)=od(k)*tol(k)*1e3;
end
%---
nTal=length(Tal);
kopo1=0;
kopo2=0;
kopo3=0;
kopo4=0;
kopo5=0;
kopo6=0;
for k=1:nTal 
    kopo1(k)=-log(Tal(k))/xv(1);
    kopo2(k)=-log(Tal2(k))/xv(2);
    kopo3(k)=-log(Tal3(k))/xv(3);
    kopo4(k)=-log(Tal4(k))/xv(4);
    kopo5(k)=-log(Tal5(k))/xv(5);
    kopo6(k)=-log(Tal6(k))/xv(6); 
end
kopo7=0;
kopo8=0;
kopo9=0;
%dvmakopo=obrazDvuMas(kopo1,kopo2,kopo3,kopo4,kopo5,kopo6,kopo7,kopo8,kopo9,0);
%alsre=postGrafAl(dvmakopo,1e2*od,3,2,length(kopo1));
alsre=vyvodZnachMas(mkbr,mv,tol,xv,md,od,rov,vkbr,vv);
%kop1=utochAlpha(Tal,kopo1,6,xv(1));
%kop2=utochAplha(Tal2,kop2,6,xv2);
%kop3=utochAplha(Tal3,kop3,6,xv3);
%kop4=utochAplha(Tal4,kop4,6,xv4);
%kop5=utochAplha(Tal5,kop5,6,xv5);
%kop6=utochAplha(Tal5,kop6,6,xv6);
alsr1=0;
for k=1:nTal 
    alsr1(k)=kopo1(k)+kopo2(k)+kopo3(k)+kopo4(k);
    alsr1(k)=alsr1(k)+kopo5(k)+kopo6(k);
    alsr1(k)=alsr1(k)/6;
end
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=privedkEdiPropus(TrkoVer5421()); 
Tal2=privedkEdiPropus(TrkoVer5422()); 
mkbr=0; mkbr(1)=250; 
mv=0; mv(1)=0.255; 
tol=0; tol(1)=0.72; 
rokbr=2.75; 
rov=0; rov(1)=0.52; 
%---
mkbr(2)=249.629; 
mv(2)=0.539; 
tol(2)=0.71; 
rov(2)=rov(1); 
%---
Tal3=privedkEdiPropus(TrkoVer5423()); 
Tal4=privedkEdiPropus(TrkoVer5431());
mkbr(3)=249.294; 
mv(3)=0.809; 
tol(3)=0.7;  
%---
mkbr(4)=249.706; 
mv(4)=0.295; 
tol(4)=0.7; 
rov(3)=0.52; 
rov(4)=rov(3);
%---
Tal5=privedkEdiPropus(TrkoVer5432()); 
Tal6=privedkEdiPropus(TrkoVer5433());
mkbr(5)=249.510; 
mv(5)=0.517; 
tol(5)=0.73; 
rokbr=2.75; 
rov(5)=0.52; 
%---
mkbr(6)=249.307; 
mv(6)=0.756; 
tol(6)=0.72; 
rov(6)=rov(5);
%---
Tal7=privedkEdiPropus(TrkoVer5441()); 
Tal8=privedkEdiPropus(TrkoVer5442());
Tal9=privedkEdiPropus(TrkoVer5443());
mkbr(7)=250.328; 
mv(7)=0.36; 
tol(7)=0.74; 
rov(7)=0.52;
%---
mkbr(8)=249.604; 
mv(8)=0.534; 
tol(8)=0.7; 
rov(8)=rov(7);
%---
mkbr(9)=249.206; 
mv(9)=0.843; 
tol(9)=0.76; 
rov(9)=0.52;
%---
md=0; od=0; xv=0; vkbr=0; vv=0; 
for k=1:9
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov(k)); 
    md(k)=mv(k)/(mv(k)+mkbr(k));
    od(k)=(vv(k)/(vv(k)+vkbr(k)));
    xv(k)=od(k)*tol(k)*1e3;
end
for k=1:nTal
    kopo1(k)=-log(Tal(k))/xv(1);
    kopo2(k)=-log(Tal2(k))/xv(2);
    kopo3(k)=-log(Tal3(k))/xv(3);
    kopo4(k)=-log(Tal4(k))/xv(4);
    kopo5(k)=-log(Tal5(k))/xv(5);
    kopo6(k)=-log(Tal6(k))/xv(6);
    kopo7(k)=-log(Tal7(k))/xv(7);
    kopo8(k)=-log(Tal8(k))/xv(8);
    kopo9(k)=-log(Tal9(k))/xv(9); 
end;
alsr2=0;
for k=1:nTal 
    alsr2(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k);
    alsr2(k)=alsr2(k)+Tal6(k)+Tal7(k)+Tal8(k)+Tal9(k);
    alsr2(k)=alsr2(k)/9;
end
%dvmakopo=obrazDvuMas(kopo1,kopo2,kopo3,kopo4,kopo5,kopo6,kopo7,kopo8,kopo9,1);
%alsre=postGrafAl(dvmakopo,1e2*od,5,2,length(kopo1));
alsre=vyvodZnachMas(mkbr,mv,tol,xv,md,od,rov,vkbr,vv);
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56;
Tal=privedkEdiPropus(TrkoVer5551()); 
Tal2=privedkEdiPropus(TrkoVer5552()); 
mkbr=0; mkbr(1)=249.913; 
mv=0; mv(1)=0.315; 
tol=0; tol(1)=0.74; 
rokbr=2.75; 
rov=0; rov(1)=0.53; 
%---
mkbr(2)=249.607; 
mv(2)=0.473; 
tol(2)=0.74; 
rov(2)=rov(1); 
%---
Tal3=privedkEdiPropus(TrkoVer5553()); 
Tal4=privedkEdiPropus(TrkoVer5561());
mkbr(3)=249.218; 
mv(3)=0.709; 
tol(3)=0.72; 
rov(3)=0.53; 
%---
mkbr(4)=249.929; 
mv(4)=0.293; 
tol(4)=0.72; 
rov(4)=rov(3); 
%---
Tal5=privedkEdiPropus(TrkoVer5562()); 
Tal6=privedkEdiPropus(TrkoVer5563());
mkbr(5)=249.695; 
mv(5)=0.528; 
tol(5)=0.71; 
rokbr=2.75; 
rov(5)=0.53; 
%---
mkbr(6)=249.306; 
mv(6)=0.83; 
tol(6)=0.7; 
rov(6)=rov(5); 
%---
Tal7=privedkEdiPropus(TrkoVer5571()); 
Tal8=privedkEdiPropus(TrkoVer5572());
mkbr(7)=250.405; 
mv(7)=0.27; 
tol(7)=0.78; 
rov(7)=0.53;
%---
mkbr(8)=249.625; 
mv(8)=0.493; 
tol(8)=0.73; 
rov(8)=rov(7);
%---
Tal9=privedkEdiPropus(TrkoVer5573());
mkbr(9)=249.348; 
mv(9)=0.764; 
tol(9)=0.76; 
rov(9)=0.53;
%---
md=0; od=0; xv=0; vkbr=0; vv=0; 
for k=1:9
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov(k)); 
    md(k)=mv(k)/(mv(k)+mkbr(k));
    od(k)=(vv(k)/(vv(k)+vkbr(k)));
    xv(k)=od(k)*tol(k)*1e3;
end
for k=1:nTal
    kopo1(k)=-log(Tal(k))/xv(1);
    kopo2(k)=-log(Tal2(k))/xv(2);
    kopo3(k)=-log(Tal3(k))/xv(3);
    kopo4(k)=-log(Tal4(k))/xv(4);
    kopo5(k)=-log(Tal5(k))/xv(5);
    kopo6(k)=-log(Tal6(k))/xv(6);
    kopo7(k)=-log(Tal7(k))/xv(7);
    kopo8(k)=-log(Tal8(k))/xv(8);
    kopo9(k)=-log(Tal9(k))/xv(9); 
end
alsr3=0;
for k=1:nTal 
    alsr3(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k);
    alsr3(k)=alsr3(k)+Tal6(k)+Tal7(k)+Tal8(k)+Tal9(k);
    alsr3(k)=alsr3(k)/9;
end
%dvmakopo=obrazDvuMas(kopo1,kopo2,kopo3,kopo4,kopo5,kopo6,kopo7,kopo8,kopo9,1);
%alsre=postGrafAl(dvmakopo,1e2*od,5,2,length(kopo1));
alsre=vyvodZnachMas(mkbr,mv,tol,xv,md,od,rov,vkbr,vv);
kopo7=0;
kopo8=0;
kopo9=0;
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=0; Tal=privedkEdiPropus(TrkoVer5681()); 
Tal2=0; Tal2=privedkEdiPropus(TrkoVer5682());
mkbr=0; mkbr(1)=250.882; 
mv=0; mv(1)=0.320; 
tol=0; tol(1)=0.76;
rokbr=2.75; 
rov=0; rov(1)=0.56; 
mkbr(2)=249.590; 
mv(2)=0.533; 
tol(2)=0.72; 
rov(2)=rov(1); 
%---
Tal3=0; Tal3=privedkEdiPropus(TrkoVer5683()); 
Tal4=0; Tal4=privedkEdiPropus(TrkoVer5691());
mkbr(3)=249.213; 
mv(3)=0.849; 
tol(3)=0.69; 
mkbr(4)=250.299; 
mv(4)=0.223; 
tol(4)=0.73; 
rov(3)=0.56; 
rov(4)=rov(3); 
%---
Tal5=privedkEdiPropus(TrkoVer5692()); 
Tal6=privedkEdiPropus(TrkoVer5693());
mkbr(5)=249.441; 
mv(5)=0.502; 
tol(5)=0.73; 
rokbr=2.75; 
rov(5)=0.56; 
mkbr(6)=249.365; 
mv(6)=0.797; 
tol(6)=0.73; 
rov(6)=rov(5); 
%---
md=0; od=0; xv=0; vkbr=0; vv=0; 
for k=1:6
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov(k)); 
    md(k)=mv(k)/(mv(k)+mkbr(k));
    od(k)=(vv(k)/(vv(k)+vkbr(k)));
    xv(k)=od(k)*tol(k)*1e3;
end
for k=1:nTal
    kopo1(k)=-log(Tal(k))/xv(1);
    kopo2(k)=-log(Tal2(k))/xv(2);
    kopo3(k)=-log(Tal3(k))/xv(3);
    kopo4(k)=-log(Tal4(k))/xv(4);
    kopo5(k)=-log(Tal5(k))/xv(5);
    kopo6(k)=-log(Tal6(k))/xv(6); 
end
alsr4=0;
for k=1:nTal 
    alsr4(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k)+Tal6(k);
    alsr4(k)=alsr4(k)/6;
end
%dvmakopo=obrazDvuMas(kopo1,kopo2,kopo3,kopo4,kopo5,kopo6,kopo7,kopo8,kopo9,0);
%alsre=postGrafAl(dvmakopo,1e2*od,3,2,length(kopo1));
alsre=vyvodZnachMas(mkbr,mv,tol,xv,md,od,rov,vkbr,vv);
alsre=0;
for k=1:nTal
    alsre(k)=(alsr1(k)+alsr2(k)+alsr3(k)+alsr4(k))/4;
end
alsr=0;
end

function [ prked ] = privedkEdiPropus(ar)
arr=ar;
p=length(arr);
for k=1:p
    if (arr(k)>1)
        arr(k)=1;
    end
end
prked=arr;
end

function [ PokPrel ] = Kramers_Kron(alsr,Sp)
format long g;
c0=299792458;
enu=1e-3;
leSp=length(Sp);
nu=0; nus=0; hi=0;
for k=1:leSp
    hi(k)=Sp(k)*alsr(k)/4/pi;
    nu(k)=2*pi*c0/Sp(k);
end
nus=izmMasChast(nu,enu);
nnus=PokazPrelomAl(nu,nus,alsr,c0)-1+1.563;
nnu=PokazPrelomPribl(nu,nus,nnus);
PokPrel=nnu;
end

function [ oprn ] = opredpokprel(nu,nus,hi,epsSt,leSp)
n=0;
for k=1:leSp
f=0;
for j=1:leSp
    if (nu(j)==nus(k))
        f(j)=0;
    else
        f(j)=(hi(j)*nu(j))/(nu(j)^2-nus(k)^2);
    end
end
n(k)=sqrt(epsSt)+(2/pi)*integpo2ma(nu,f);
end
oprn=n;
end

function inte = integpo2ma(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ nus ] = izmMasChast(nu,enu)
nust=0;
lenu=length(nu);
for k=1:lenu-1
    nust(k)=enu*(nu(k+1)-nu(k))+nu(k);
end
nust(lenu)=2*nust(lenu-1)-nust(lenu-2);
nus=nust;
end

function [ pop ] = PokazPrelomAl(nu,nus,ka,vl)
np=0;
p=length(nu);
for k=1:p
    fn=0; q=1;
    for j=1:p-1
        dkpodo=(ka(j+1)-ka(j))/(nu(j+1)-nu(j));
        podln=((nu(j)+nus(k))/(nu(j)-nus(k)));
        podln=abs(podln);
        fn(q)=dkpodo*log(podln);
        q=q+1;
    end
    fn(p)=fn(p-1);
    np(k)=1+(vl/pi)*integpo2ma(nu,fn)/2/nu(k);
end
pop=np;
end

function [ ns ]  = PokazPrelomPribl(nu,nus,na)
nnu=0;
p=length(nu);
for k=1:p-1
    de=(-nus(k)+nu(k))*(na(k+1)-na(k))/(nu(k+1)-nu(k));
    nnu(k)=de+na(k);
end
nnu(p)=(-nus(p)+nu(p))*(na(p)-na(p-1))/(nu(p)-nu(p-1))+na(p);
ns = nnu;
end

function [ sred ] = SredGraSt()
Tal=privedkEdiPropus(Trko2()); Tal2=privedkEdiPropus(Trko4());
mkbr=232.188; mv=1.041; tol=0.67; rokbr=2.75; rov=0.112; 
mkbr2=228.831; mv2=0.979; tol2=0.64; rov2=0.092; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); 
vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
xv=(vv/(vv+vkbr))*tol*1e3; 
xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
nTal=length(Tal);
for k=1:nTal 
    Tal(k)=-log(Tal(k))/xv; 
    Tal2(k)=-log(Tal2(k))/xv2;
end
alsr1=0;
for k=1:nTal 
    alsr1(k)=(Tal(k)+Tal2(k))/2;
end
sred=alsr1;
end

function [ obma ] = skleiv2Mas(ar1,ar2)
q=1;
obm=0;
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
dl1=0; dl1=dlvoi(); 
dl2=0; dl2=dlvoVer53101();
kon=dl2(1);
p=length(dl1);
q=1;
for k=1:p
    if (dl1(k)>kon)
        q=q+1;
    else
        break;
    end
end
t1=(ar1(q)-ar1(q-1))*(dl2(1)-dl1(q))/(dl1(q)-dl1(q-1));
t1=ar1(q-1)+t1;
t2=ar2(1);
t2=(t2-t1);
arr=0;
for k=1:q+1
    if (k<(q+1))
    arr(k)=ar1(k)+t2;
    end
end
dlv = arr;
end

function [ oa ] = utochAlpha(T,alp,n,d)
Nh=1e3; eps=1e-2; 
ap=0; ak=0; To=0;
To=oprTrans(T);
ap=1e6*SredGraSt();
alsr=soglMas(ap,1e6*alp);
alsr=skleiv2Mas(alsr,1e6*alp);
ap=0; ap=alsr;
dl=RasshDiapDlinVoln();
    pp=Kramers_Kron(ap,dl);
    ak=oprAlp(ap,pp,To,dl,d);
    ep=oprMax(RaznMasi(ak,ap))/oprMax(ap);
    m1=ap; m2=ak;
    h=1; 
    while (ep>eps)
        ap=0; ap=ak;
    if (h>Nh)
        break;
    end
    To=oprTransN(ap,d);
    pp=Kramers_Kron(ap,dl);
    ak=oprAlp(ap,pp,To,dl,d);
    ep=oprMax(RaznMasi(ak,ap))/oprMax(ap);
    h=h+1;
    end
    m3=ak;
    h=h*1
oa=ak;
end

function R = oprKoeOtr(n,al,la)
hi=(al*la)/(4*pi);
R=((n-1)^2+hi^2)/((n+1)^2+hi^2);
end

function [ al ] = oprAlp(apr,npp,Tr,dl,d)
n=oprMin([length(npp)   length(Tr) length(dl)]);
alp=0;
for k=1:n
    R=oprKoeOtr(npp(k),apr(k),dl(k));
    alp(k)=RootsFind(apr(k),d,Tr(k),R);
end
al=alp;
end

function kor = RootsFind(as,d,t,r)
aa=0; ab=as+1e7; nh=1e4;
eps=1; e=1e-7; h=0;
while (eps>e)
ac=(aa+ab)/2;    
fc=(1-r^2)*exp(-ac*d)/(1-(r*exp(-ac*d))^2)-t;
fa=(1-r^2)*exp(-aa*d)/(1-(r*exp(-aa*d))^2)-t;
fb=(1-r^2)*exp(-ab*d)/(1-(r*exp(-ab*d))^2)-t;
if ((fa*fc)<0)
    if ((fb*fc)>0)
    ab=ac; 
    end
end
if ((fa*fc)>0) 
    if ((fb*fc)<0) 
    aa=ac; 
    end
end
eps=abs(fa-fb);
h=h+1;
if (h>nh) 
    break; 
end
end
kor=ac;
end

function [ pg ] = postrGraf(la1,la2,al1,al2)
pl=plot(la1,al1,'-r',la2,al2,'-k');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'}); 
ylabel({'Коэффициент поглощения, м^(-1)'}); 
title({'График зависимости коэффициента поглощения от длины волны'});
pg=[0];
end

function m = oprMax(ar)
n=length(ar);
ma=ar(1);
for k=1:n
    if (ar(k)>ma)
        ma=ar(k);
    end
end
m=ma;
end

function m = oprMin(ar)
n=length(ar);
mi=ar(1);
for k=1:n
    if (ar(k)<mi)
        mi=ar(k);
    end
end
m=mi;
end

function [ ar ] = RaznMasi(ar1,ar2)
n=oprMin([length(ar1)   length(ar2)]);
ar3=0;
for k=1:n
    ar3(k)=ar1(k)-ar2(k);
    ar3(k)=abs(ar3(k));
end
ar=ar3;
end

function [ tr ] = oprTrans(T)
al=SredGraSt();
mkbr=232.188; mv=1.041; tol=0.67; rokbr=2.75; rov=0.112; 
mkbr2=228.831; mv2=0.979; tol2=0.64; rov2=0.092; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); 
vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
xv=(vv/(vv+vkbr))*tol*1e3; 
xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
p=length(al);
Tt=0; xsr=(xv+xv2)/2;
for k=1:p
    Tt(k)=exp(-al(k)*xsr);
end
Tt=soglMas(Tt,T);
To=skleiv2Mas(Tt,T);
dl1=RasshDiapDlinVoln();
dl2=dlivoln();
dl1=utochDlin(dl1,To,1);
To=utochDlin(dl1,To,2);
h=postrGraf(dl1,dl2,To,T); 
tr=To;
end

function [ tr ] = oprTransN(al,d)
n=length(al);
t=0;
for k=1:n
    t(k)=exp(-al(k)*d);
end
tr=t;
end

function [ m ] = utochDlin(a1,a2,p)
n=oprMin([length(a1)   length(a2)]);
ar=0;
if (p==1)
for k=1:n
    ar(k)=a1(k);
end
else 
    for k=1:n
    ar(k)=a2(k);
    end
end
m=ar;
end

function [ pg ] = postGrafAl(dvmakp,to,n,k,u)
if (k==1)
    p=1;
end
if (n==3)
if (k>1)
    p=4;
    end
end
if (n>3)
if (k>1)
    p=5;
    end
end
la=1e6*dlivoln(); 
tol=0;
if (n==3)
    al1=perepVOdnMas(dvmakp,p,u);
    al2=perepVOdnMas(dvmakp,p+1,u);
    al3=perepVOdnMas(dvmakp,p+2,u);
pl=plot(la,al1,'-b',la,al2,'-k',la,al3,'-m');
for l=1:n
    tol(l)=to(p+l-1);
end
sto1=sprintf('%1.3f',tol(1)); 
sto1=strcat(sto1,' %');
%sto1=strcat(sto1,' мкм');
sto2=sprintf('%1.3f',tol(2));
sto2=strcat(sto2,' %');
%sto2=strcat(sto2,' мкм');
sto3=sprintf('%1.3f',tol(3));
sto3=strcat(sto3,' %');
%sto3=strcat(sto3,' мкм');
legend(sto1,sto2,sto3,'location','best');
end
if (n==4)
    if (k==1)
    al1=perepVOdnMas(dvmakp,p,u);
    al2=perepVOdnMas(dvmakp,p+1,u);
    al3=perepVOdnMas(dvmakp,p+2,u);
    al4=perepVOdnMas(dvmakp,p+3,u);
pl=plot(la,al1,'-k',la,al2,'-b',la,al3,'-m',la,al4,'-r');
for l=1:n
    tol(l)=to(p+l-1);
end
sto1=sprintf('%1.3f',tol(1)); 
sto1=strcat(sto1,' %');
%sto1=strcat(sto1,' мкм');
sto2=sprintf('%1.3f',tol(2)); 
sto2=strcat(sto2,' %');
%sto2=strcat(sto2,' мкм');
sto3=sprintf('%1.3f',tol(3)); 
sto3=strcat(sto3,' %');
%sto3=strcat(sto3,' мкм');
sto4=sprintf('%1.3f',tol(4)); 
sto4=strcat(sto4,' %');
%sto4=strcat(sto4,' мкм');
legend(sto1,sto2,sto3,sto4,'location','best');
    end
end
if (n==5)
    if (k>1)
    al1=perepVOdnMas(dvmakp,p,u);
    al2=perepVOdnMas(dvmakp,p+1,u);
    al3=perepVOdnMas(dvmakp,p+2,u);
    al4=perepVOdnMas(dvmakp,p+3,u);
    al5=perepVOdnMas(dvmakp,p+4,u);
pl=plot(la,al1,'-k',la,al2,'-b',la,al3,'-m',la,al4,'-r',la,al5,'-g');
for l=1:n
    tol(l)=to(p+l-1);
end
sto1=sprintf('%1.3f',tol(1)); 
sto1=strcat(sto1,' %');
%sto1=strcat(sto1,' мкм');
sto2=sprintf('%1.3f',tol(2)); 
sto2=strcat(sto2,' %');
%sto2=strcat(sto2,' мкм');
sto3=sprintf('%1.3f',tol(3)); 
sto3=strcat(sto3,' %');
%sto3=strcat(sto3,' мкм');
sto4=sprintf('%1.3f',tol(4)); 
sto4=strcat(sto4,' %');
%sto4=strcat(sto4,' мкм');
sto5=sprintf('%1.3f',tol(5)); 
sto5=strcat(sto5,' %');
%sto5=strcat(sto5,' мкм');
legend(sto1,sto2,sto3,sto4,sto5,'location','best');
    end
end
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'}); 
ylabel({'Коэффициент поглощения, мкм(-1)'}); 
title({'График зависимости коэффициента поглощения от длины волны'});
pg=[0];
end

function [ mas ] = obrazDvuMas(al1,al2,al3,al4,al5,al6,al7,al8,al9,n)
p=length(al1);
dm=0;
for j=1:9
    for k=1:p
        dm(j,k)=0;
        if (n==0)
        al7(k)=n;
        al8(k)=n;
        al9(k)=n;
        end
    end
end
dm=RecordMasToDM(al1,dm,p,1);
dm=RecordMasToDM(al2,dm,p,2);
dm=RecordMasToDM(al3,dm,p,3);
dm=RecordMasToDM(al4,dm,p,4);
dm=RecordMasToDM(al5,dm,p,5);
dm=RecordMasToDM(al6,dm,p,6);
dm=RecordMasToDM(al7,dm,p,7);
dm=RecordMasToDM(al8,dm,p,8);
dm=RecordMasToDM(al9,dm,p,9);
mas=dm;
end

function [ m ] = RecordMasToDM(al,dm,p,j)
for k=1:p
    dm(j,k)=al(k);
end
m=dm;
end

function [ m ] = perepVOdnMas(dm,j,p)
for k=1:p
    al(k)=dm(j,k);
end
m=al;
end

function p = vyvodZnachMas(ar1,ar2,ar3,ar4,ar5,ar6,ar7,ar8,ar9)
for k=1:length(ar1)
s=sprintf('ar1 [ %d ] = %1.3f\n',k,ar1(k))
end
for k=1:length(ar2)
s=sprintf('ar2 [ %d ] = %1.3f\n',k,ar2(k))
end
for k=1:length(ar3)
s=sprintf('ar3 [ %d ] = %1.2f\n',k,ar3(k))
end
for k=1:length(ar4)
s=sprintf('ar4 [ %d ] = %1.3f\n',k,ar4(k))
end
for k=1:length(ar5)
s=sprintf('ar5 [ %d ] = %1.4f\n',k,1e2*ar5(k))
end
for k=1:length(ar6)
s=sprintf('ar6 [ %d ] = %1.4f\n',k,1e2*ar6(k))
end
for k=1:length(ar7)
s=sprintf('ar7 [ %d ] = %1.2f\n',k,ar7(k))
end
for k=1:length(ar8)
s=sprintf('ar8 [ %d ] = %1.4f\n',k,ar8(k))
end
for k=1:length(ar9)
s=sprintf('ar9 [ %d ] = %1.3f\n',k,1e3*ar9(k))
end
p=0;
end