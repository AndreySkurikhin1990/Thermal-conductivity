function [ PokPrel ] = Kramers_n()
format long g;
c0=299792458;
enu=1e-3;
%alsr=1e6*SredGrafKoAbsVer();
%alsr=skleiv2Mas(soglMas(1e6*SredGraf(),1e6*SredGraSt()),1e6*SredGraf());
alsr=RasMasKoAbs(); 
%n0=ZapisFileOptio(alsr);
Sp=RasshDiapDlinVoln();
leSp=min([length(alsr),length(Sp)]);
%n0=postrGraf(alsr,alph,Sp,dlinvolny());
%Sp=dlivoln();
%alsr=1e6*SredGraf();
nu=0; nus=0; hi=0;
for k=1:leSp
    hi(k)=Sp(k)*alsr(k)/4/pi;
    nu(k)=2*pi*c0/Sp(k);
end
nus=izmMasChast(nu,enu);
nnus=PokazPrelomAl(nu,nus,alsr,c0)-1+1.543;
nnu=PokazPrelomPribl(nu,nus,nnus);
%n0=postrGraf(nnu,nnus,nu,nus);
%hia=PokazPoglPribl(nus,alsr,nu);
%ep=opredEpsilSt(nu/2/pi,nus/2/pi,hi,hia)
%ns=PokazPrelomProv2(nu,nus,hi,hia,nnu,nnus);
%------
%n0=opredn0(Sp,alsr);
n0=ZapisFile(nnu);
PokPrel=nnu;
end

function t = ZapisFileOptio(massi)
fid = fopen('Koefficient_pogloscheniya.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFile(massi)
fid = fopen('Pokazatel_prelomleniya.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
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

function epsi = opredEpsilSt(nu,nus,hinu,hia)
leSp=length(nu);
eps=1e-4;
epst=0;
for k=1:leSp
euta=1+1e-6;
eutb=1e1; 
Nit=2e2;
p=hia(k)
t=1;
while (abs(eutb-euta)>eps)
eutc=(euta+eutb)/2;
ffa=0; ffb=0; ffc=0;
na=opredpokprel(nu,nus,hinu,euta,leSp);
nb=opredpokprel(nu,nus,hinu,eutb,leSp);
nc=opredpokprel(nu,nus,hinu,eutc,leSp);
     for j=1:leSp
         ffc(j)=(nc(j)-sqrt(eutc))/(nu(j)^2-nus(k)^2);
         ffa(j)=(na(j)-sqrt(euta))/(nu(j)^2-nus(k)^2);
         ffb(j)=(nb(j)-sqrt(eutb))/(nu(j)^2-nus(k)^2);
     end  
ma=-(2*nus(k)/pi)*integpo2ma(nu,ffa)
ma=ProvAdekv(ma);
ma=ma-p;
mb=-(2*nus(k)/pi)*integpo2ma(nu,ffb)
mb=ProvAdekv(mb);
mb=mb-p;
mc=-(2*nus(k)/pi)*integpo2ma(nu,ffc)
mc=ProvAdekv(mc);
mc=mc-p;
if (ma*mc>0)        
euta=eutc;    
else
if (mc*mb>0)
eutb=eutc;
else disp(eutc);
end
end
t=t+1;
if (t>Nit)
    disp(t); break; end
end
epst(k)=eutc;
disp(k);
disp(eutc);
end
end

function [ npopr ] = PokazPrelomProv(nu,nus,hi,hia,nnu,na)
p=length(hi);
for j=1:p
nmk(j)=2*na(j)*hia(j);
end
nmkp=0;
for k=1:p
    ff=0;
    t=(na(k)^2)-(hia(k)^2);
    for j=1:p
        ff(j)=(nnu(j)^2-hi(j)^2-t)/(nu(j)^2-nus(k)^2);
    end
    nmkp(k)=-(2*nus(k)/pi)*integpo2ma(nu,ff);
end
npopr=postrGraf(nmk,nmkp,nu,nus);
end

function n0 = opredn0(dl,ka)
n0=1+(1/2/pi^2)*integpo2ma(dl,ka);
end

function [ dl ] = izmMasDlVol(lam)
om=0;
c0=299792458;
lenu=length(lam);
for k=1:lenu
    om(k)=2*pi*c0/lam(k);
end
dl=om;
end

function pg = postrGraf(nnu,nnus,nu,nus)
dlnu=izmMasDlVol(nu);
dlnus=izmMasDlVol(nus); 
%dlnu=nu;
%dlnus=nus;
p=length(dlnu); 
dlnuo=0; dlnuso=0;
j=1; nnuo=0; na=0; 
for k=1:p  
    if (dlnu(k)<20e-6)    
      nnuo(j)=nnu(k); 
      dlnuo(j)=dlnu(k); 
      j=j+1; 
    end
end
p1=j-1;
p=length(dlnus);
j=1;
for k=1:p  
    if (dlnus(k)<20e-6)    
        na(j)=nnus(k);
        dlnuso(j)=dlnus(k); 
        j=j+1;
    end
end
p2=j-1;
pl=plot(dlnuo(2:p1-1)*1e6,nnuo(2:p1-1),'-b',dlnuso(2:p2-1)*1e6,na(2:p2-1),'-k');
%pl=plot(dlnuo*1e6,nnuo,'-b');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'}); 
ylabel({'Показатель преломления'}); 
title({'График зависимости показателя преломления от длины волны'});
pg=0;
end

function [ npopr ] = PokazPrelomProv2(nu,nus,hi,hia,nnu,na)
p=length(hi);
nmk=0;
for j=1:p
nmk(j)=(na(j)^2)-(hia(j)^2);
end
nmkp=0;
for k=1:p
    ff=0;
    tt=2*na(k)*hia(k)*nus(k);
    for j=1:p
        t=(nu(j)^2)-(nus(k)^2);
        ff(j)=(2*nnu(j)*hi(j)*nu(j)-tt)/t;
    end
    nmkp(k)=1+(2/pi)*integpo2ma(nu,ff);
end
npopr=postrGraf(nmk,nmkp,nu,nus);
end

function [ hia ] = PokazPoglPribl(nus,alsr,nu)
dls=0; 
c0=299792458;
p=length(nus);
als=0;
his=0;
for k=1:p-1
    de=(-nus(k)+nu(k))*(alsr(k+1)-alsr(k))/(nu(k+1)-nu(k));
    als(k)=alsr(k)+de;
end
als(p)=(-nus(p)+nu(p))*(alsr(p)-alsr(p-1))/(nu(p)-nu(p-1))+alsr(p);
for k=1:p
dls(k)=2*pi*c0/nus(k);
his(k)=dls(k)*als(k)/4/pi;
end
hia=his;
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
    alsr1(k)=Tal(k)+Tal2(k);
end
sred=alsr1/2;
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

function [ dlv ] = dlinvolny()
dl=0; 
dl=dlvoi();
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end

function [ dlv ] = konMas(ar)
dl=0; 
dl=dlvoi(); 
dlvv=0;
dlvv=dlvoVer53101();
kon=dlvv(1);
ar1=0;
p=length(ar); 
q=1;
for k=1:p  
    if (dl(k)>kon)
        ar1(q)=ar(k);
        q=q+1;
    end
end
dlv = ar1;
end

function [ dlv ] = soglMas(ar1,ar2)
dl=0; 
dl=dlvoi(); 
dlvv=0;
dlvv=dlvoVer53101();
kon=dlvv(1);
ar3=0;
p=length(dl); 
q=1;
arr=0;
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
    else
%    arr(k)=ar1(j);
%    j=j+1;
    end
end
dlv = arr;
end