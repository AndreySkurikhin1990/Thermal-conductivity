%определяет средний КП по Росселанду и по Планку во всем и в заданном интервале длин волн
function tm = tmp30()
format long g;
alsr=1e6*SredGraf();
alfs = 0; alfs = RasMasKoAbs30(alsr);
    pp=0; pp=Kramers_n30(alfs);
    dl=0; dl=RasshDiapDlinVoln30();
T = 3e2:1e3:4e2;
al=0; alp=0;
for k=1:length(T)
al(k) = 1e6 / sredRosSieg30(T(k),alfs,dl,pp);
alp(k) = 1e6 / sredRosSieg(T(k));
end
al=1*al'
al=0; 
for k=1:length(T)
al(k)=AlphaSredDVVSr30(dl,pp,alfs,T(k));
end
al=1*al'
dv = SuzhDiapDlVo(1e6*dl);
popr=SuzhDiapDanPoPr(1e6*dl,pp);
kop=SuzhDiapDanKoPo(1e6*dl,alfs);
al=0; 
for k=1:length(T)
al(k) = 1e6 / sredRosSieg30(T(k),kop,dv,popr);
end
al=1*al'
al=0;
for k=1:length(T)
al(k)=AlphaSredDVVSr30(dv,popr,kop,T(k));
end
al=1*al'
tm = alp';
%tm = KoPogTabl26(1e6*dlivoln(), Tal);
end

function [ kp ] = KoPogTabl30(wl, alsr)
dl = 15e-1:5e-1:15e0;
nd=length(dl);
j = 1;
kopo = 0;
n = length(wl);
kopo = 0;
for j=1:nd
    for k=1:n
        if (wl(k) > dl(j))
        kopo(j) = alsr(k);
        %disp(wl(k));
        break;
        end
    end
end
kp = kopo';
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

function [ rs ] = sredRosSieg30(tem,alfs,dl,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4));
np=nsreddvvsr30(dl,npp,tem);
Cons=(pi/2)*(C1*C2)/sig;
%temp=2e2:1e2:1e3; temp=temp+te0; 
Ibc=0; Ibz=0; chasc=0; chacz=0; 
dlsvprfo1=0; dlsvprfo2=0; me=0; 
ote=0; chi=0; zna=0; c1m=0; c2m=0;
for k=1:length(npp)
c1m(k)=C1/(npp(k)^2);
c2m(k)=C2/npp(k);
Cons(k)=(pi/2)*(c1m(k)*c2m(k))/sig;
dlv(k)=dl(k)/npp(k);
%dlv(k)=dl(k);
%me=exp(C2/(dlv(k)*tem))-1;Ib(k)=2*pi*C1/(dlv(k)^5)/me;ote=(sig/Ib(k))^(1/4);
eb=(np^2)*sig*(tem^4);
ote=(sig*(np^2)/eb)^(1/4);
chi=exp(c2m(k)*ote/dlv(k));
zna=(chi-1)^2;
Ibz(k)=Cons(k)*(ote^5)*chi/zna/(dlv(k)^6); 
Ibc(k)=Ibz(k)/alfs(k);
end
chasc=integpo2ma30(dlv,Ibc);
chasz=integpo2ma30(dlv,Ibz);
dlsvprfo2=chasc/chasz;
rs=1e6*dlsvprfo2;
end

function [ ndv ] = RasshDiapDlinVoln30()
Sp=skleiv2Mas30(dlinvolndoMas30(),dlinvoln30());
ndv=Sp;
end

function [ dlv ] = dlinvoln30()
dl=0; 
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end

function [ dlv ] = dlinvolndoMas30()
dl=0; 
dl=dlvoi();
wl=dlvoVer53101();
kon=wl(1);
dli=0;
p=length(dl); 
q=1;
for k=1:p  
    if (dl(k)>kon)
    dli(q)=1e-2/dl(k);
    q=q+1;
    end
end
dlv = dli;
end

function [ obma ] = skleiv2Mas30(ar1,ar2)
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

function [ rasmaskopo ] = RasMasKoAbs30(alph)
alsr=skleiv2Mas30(soglMas30(alph,SredGraSt30()),alph);
rasmaskopo=alsr;
end

function [ dlv ] = soglMas30(ar1,ar2)
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
    end
end
dlv = arr;
end

function [ sred ] = SredGraSt30()
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
sred=1e6*alsr1;
end

function ns = AlphaSredDVVSr30(dv,npp,alph,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
vl=c0;
for k=1:length(npp)
    vl=vl/npp(k);
    c1=PP*c0^2;
    c2=PP*c0/PB;
    %lambda=dv(k);
    lambda=dv(k)/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
    %Ib(k)=Ib(k)/(npp(k)^2);    
Ibn(k)=alph(k)*Ib(k);
vl=c0;
end
nc=integpo2ma30(dv,Ibn);
ns=abs(real(nc/integpo2ma30(dv,Ib)));
end

function [ dv ] = SuzhDiapDlVo(dlv)
nac=3;
kon=10;
le=length(dlv);
q=1;
dl=0;
for k=1:le
    if (dlv(k)>nac)
        if (dlv(k)<kon)
            dl(q)=dlv(k);
            q=q+1;
        end
    end
end
dv = dl;
end

function [ dv ] = SuzhDiapDanPoPr(dlv,npp)
nac=3;
kon=10;
le=length(dlv);
q=1;
np=0;
for k=1:le
    if (dlv(k)>nac)
        if (dlv(k)<kon)
            np(q)=npp(k);
            q=q+1;
        end
    end
end
dv = np;
end

function [ dv ] = SuzhDiapDanKoPo(dlv,alp)
nac=3;
kon=10;
le=length(dlv);
q=1;
alph=0;
for k=1:le
    if (dlv(k)>nac)
        if (dlv(k)<kon)
            alph(q)=alp(k);
            q=q+1;
        end
    end
end
dv = alph;
end

function [ PokPrel ] = Kramers_n30(alp)
c0=299792458;
enu=1e-3;
alsr=skleiv2Mas30(soglMas30(alp,SredGraSt30()),alp);
Sp=RasshDiapDlinVoln30();
leSp=length(Sp);
nu=0; nus=0; hi=0;
for k=1:leSp
    hi(k)=Sp(k)*alsr(k)/4/pi;
    nu(k)=2*pi*c0/Sp(k);
end
nus=izmMasChast30(nu,enu);
nnus=PokazPrelomAl30(nu,nus,alsr,c0)-1+1.563;
nnu=PokazPrelomPribl30(nu,nus,nnus);
PokPrel=nnu;
end

function inte = integpo2ma30(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ nus ] = izmMasChast30(nu,enu)
nust=0;
lenu=length(nu);
for k=1:lenu-1
    nust(k)=enu*(nu(k+1)-nu(k))+nu(k);
end
nust(lenu)=2*nust(lenu-1)-nust(lenu-2);
nus=nust;
end

function [ pop ] = PokazPrelomAl30(nu,nus,ka,vl)
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
    np(k)=1+(vl/pi)*integpo2ma30(nu,fn)/2/nu(k);
end
pop=np;
end

function [ ns ]  = PokazPrelomPribl30(nu,nus,na)
nnu=0;
p=length(nu);
for k=1:p-1
    de=(-nus(k)+nu(k))*(na(k+1)-na(k))/(nu(k+1)-nu(k));
    nnu(k)=de+na(k);
end
nnu(p)=(-nus(p)+nu(p))*(na(p)-na(p-1))/(nu(p)-nu(p-1))+na(p);
ns = nnu;
end

function ns = nsreddvvsr30(dv,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
vl=c0;
for k=1:length(npp)
    vl=vl/npp(k);
    c1=PP*c0^2;
    c2=PP*c0/PB;
    %lambda=dv(k);
    lambda=dv(k)/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
    %Ib(k)=Ib(k)/(npp(k)^2);    
Ibn(k)=npp(k)*Ib(k);
vl=c0;
end
nc=integpo2ma30(dv,Ibn);
ns=abs(real(nc/integpo2ma30(dv,Ib)));
end