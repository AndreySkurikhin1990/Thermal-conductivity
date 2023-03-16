function tm = tmp26()
format long g;
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=privedkEdiPropus(TrkoVer5311());
Tal2=privedkEdiPropus(TrkoVer5312());
Tal3=privedkEdiPropus(TrkoVer5313()); 
Tal4=privedkEdiPropus(TrkoVer53101()); 
Tal5=privedkEdiPropus(TrkoVer53102()); 
Tal6=privedkEdiPropus(TrkoVer53103()); 
mkbr=250.239; mv=0.283; tol=0.73; rokbr=2.75; rov=0.49; 
mkbr2=249.740; mv2=0.464; tol2=0.72; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); 
vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
mkbr3=249.223; mv3=0.812; tol3=0.72; 
mkbr4=250.395; mv4=0.22; tol4=0.72; rov3=0.49; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); 
vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
mkbr5=250.336; mv5=0.547; tol5=0.71; rokbr=2.75; rov5=0.49; 
mkbr6=249.55; mv6=0.777; tol6=0.7; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); 
vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
xv=(vv/(vv+vkbr))*tol*1e3;
xv2=(vv2/(vv2+vkbr2))*tol2*1e3;
xv3=(vv3/(vv3+vkbr3))*tol3*1e3;
xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3;
xv6=(vv6/(vv6+vkbr6))*tol6*1e3;
nTal=length(Tal);
for k=1:nTal 
    Tal(k)=-log(Tal(k))/xv; 
    Tal2(k)=-log(Tal2(k))/xv2;
    Tal3(k)=-log(Tal3(k))/xv3;
    Tal4(k)=-log(Tal4(k))/xv4;
    Tal5(k)=-log(Tal5(k))/xv5;
    Tal6(k)=-log(Tal6(k))/xv6; 
end
alsr1=0;
for k=1:nTal 
    alsr1(k)=(Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k)+Tal6(k))/6;
end
alfs1 = 0; alfs1 = RasMasKoAbs26(1e6*Tal);
alfs2 = 0; alfs2 = RasMasKoAbs26(1e6*Tal2);
alfs3 = 0; alfs3 = RasMasKoAbs26(1e6*Tal3);
alfs4 = 0; alfs4 = RasMasKoAbs26(1e6*Tal4);
alfs5 = 0; alfs5 = RasMasKoAbs26(1e6*Tal5);
alfs6 = 0; alfs6 = RasMasKoAbs26(1e6*Tal6);
    pp1=0; pp1=Kramers_n26(alfs1);
    pp2=0; pp2=Kramers_n26(alfs2);
    pp3=0; pp3=Kramers_n26(alfs3);
    pp4=0; pp4=Kramers_n26(alfs4);
    pp5=0; pp5=Kramers_n26(alfs5);
    pp6=0; pp6=Kramers_n26(alfs6);
    dl=0; dl=RasshDiapDlinVoln26();
    pp=0; pp=Kramers_n26(alsr1);
    %pp=Kramers_n();
T = 3e2:1e3:4e2;
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0; alp=0;
for k=1:length(T)
al1(k) = 1e6 / sredRosSieg26(T(k),alfs1,dl,pp1);
al2(k) = 1e6 / sredRosSieg26(T(k),alfs2,dl,pp2);
al3(k) = 1e6 / sredRosSieg26(T(k),alfs3,dl,pp3);
al4(k) = 1e6 / sredRosSieg26(T(k),alfs4,dl,pp4);
al5(k) = 1e6 / sredRosSieg26(T(k),alfs5,dl,pp5);
al6(k) = 1e6 / sredRosSieg26(T(k),alfs6,dl,pp6);
alp(k) = 1e6 / sredRosSieg(T(k));
end
al1=1*al1'
al2=1*al2'
al3=1*al3'
al4=1*al4'
al5=1*al5'
al6=1*al6'
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0;
for k=1:length(T)
al1(k)=AlphaSredDVVSr26(dl,pp1,alfs1,T(k));
al2(k)=AlphaSredDVVSr26(dl,pp2,alfs2,T(k));
al3(k)=AlphaSredDVVSr26(dl,pp3,alfs3,T(k));
al4(k)=AlphaSredDVVSr26(dl,pp4,alfs4,T(k));
al5(k)=AlphaSredDVVSr26(dl,pp5,alfs5,T(k));
al6(k)=AlphaSredDVVSr26(dl,pp6,alfs6,T(k));
end
al1=1*al1'
al2=1*al2'
al3=1*al3'
al4=1*al4'
al5=1*al5'
al6=1*al6'
dv = SuzhDiapDlVo(1e6*dl);
popr1=SuzhDiapDanPoPr(1e6*dl,pp1);
kop1=SuzhDiapDanKoPo(1e6*dl,alfs1);
popr2=SuzhDiapDanPoPr(1e6*dl,pp2);
kop2=SuzhDiapDanKoPo(1e6*dl,alfs2);
popr3=SuzhDiapDanPoPr(1e6*dl,pp3);
kop3=SuzhDiapDanKoPo(1e6*dl,alfs3);
popr4=SuzhDiapDanPoPr(1e6*dl,pp4);
kop4=SuzhDiapDanKoPo(1e6*dl,alfs4);
popr5=SuzhDiapDanPoPr(1e6*dl,pp5);
kop5=SuzhDiapDanKoPo(1e6*dl,alfs5);
popr6=SuzhDiapDanPoPr(1e6*dl,pp6);
kop6=SuzhDiapDanKoPo(1e6*dl,alfs6);
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0;
for k=1:length(T)
al1(k) = 1e6 / sredRosSieg26(T(k),kop1,dv,popr1);
al2(k) = 1e6 / sredRosSieg26(T(k),kop2,dv,popr2);
al3(k) = 1e6 / sredRosSieg26(T(k),kop3,dv,popr3);
al4(k) = 1e6 / sredRosSieg26(T(k),kop4,dv,popr4);
al5(k) = 1e6 / sredRosSieg26(T(k),kop5,dv,popr5);
al6(k) = 1e6 / sredRosSieg26(T(k),kop6,dv,popr6);
end
al1=1*al1'
al2=1*al2'
al3=1*al3'
al4=1*al4'
al5=1*al5'
al6=1*al6'
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0;
for k=1:length(T)
al1(k)=AlphaSredDVVSr26(dv,popr1,kop1,T(k));
al2(k)=AlphaSredDVVSr26(dv,popr2,kop2,T(k));
al3(k)=AlphaSredDVVSr26(dv,popr3,kop3,T(k));
al4(k)=AlphaSredDVVSr26(dv,popr4,kop4,T(k));
al5(k)=AlphaSredDVVSr26(dv,popr5,kop5,T(k));
al6(k)=AlphaSredDVVSr26(dv,popr6,kop6,T(k));
end
al1=1*al1'
al2=1*al2'
al3=1*al3'
al4=1*al4'
al5=1*al5'
al6=1*al6'
tm = alp';
%tm = KoPogTabl26(1e6*dlivoln(), Tal);
end

function [ kp ] = KoPogTabl26(wl, alsr)
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

function [ rs ] = sredRosSieg26(tem,alfs,dl,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4));
np=nsreddvvsr26(dl,npp,tem);
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
chasc=trapz(dlv,Ibc);
chasz=trapz(dlv,Ibz);
dlsvprfo2=chasc/chasz;
rs=1e6*dlsvprfo2;
end

function [ ndv ] = RasshDiapDlinVoln26()
Sp=skleiv2Mas26(dlinvolndoMas26(),dlinvoln26());
ndv=Sp;
end

function [ dlv ] = dlinvoln26()
dl=0; 
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end

function [ dlv ] = dlinvolndoMas26()
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

function [ obma ] = skleiv2Mas26(ar1,ar2)
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

function [ rasmaskopo ] = RasMasKoAbs26(alph)
alsr=skleiv2Mas26(soglMas26(alph,1e6*SredGraSt26()),alph);
rasmaskopo=alsr;
end

function [ dlv ] = soglMas26(ar1,ar2)
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

function [ sred ] = SredGraSt26()
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

function ns = AlphaSredDVVSr26(dv,npp,alph,tem)
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
nc=trapz(dv,Ibn);
ns=abs(real(nc/trapz(dv,Ib)));
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

function [ PokPrel ] = Kramers_n26(alp)
format long g;
c0=299792458;
enu=1e-3;
alsr=skleiv2Mas26(soglMas26(alp,1e6*SredGraSt26()),alp);
Sp=RasshDiapDlinVoln26();
leSp=length(Sp);
nu=0; nus=0; hi=0;
for k=1:leSp
    hi(k)=Sp(k)*alsr(k)/4/pi;
    nu(k)=2*pi*c0/Sp(k);
end
nus=izmMasChast26(nu,enu);
nnus=PokazPrelomAl26(nu,nus,alsr,c0)-1+1.563;
nnu=PokazPrelomPribl26(nu,nus,nnus);
PokPrel=nnu;
end

function inte = integpo2ma26(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ nus ] = izmMasChast26(nu,enu)
nust=0;
lenu=length(nu);
for k=1:lenu-1
    nust(k)=enu*(nu(k+1)-nu(k))+nu(k);
end
nust(lenu)=2*nust(lenu-1)-nust(lenu-2);
nus=nust;
end

function [ pop ] = PokazPrelomAl26(nu,nus,ka,vl)
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
    np(k)=1+(vl/pi)*integpo2ma26(nu,fn)/2/nu(k);
end
pop=np;
end

function [ ns ]  = PokazPrelomPribl26(nu,nus,na)
nnu=0;
p=length(nu);
for k=1:p-1
    de=(-nus(k)+nu(k))*(na(k+1)-na(k))/(nu(k+1)-nu(k));
    nnu(k)=de+na(k);
end
nnu(p)=(-nus(p)+nu(p))*(na(p)-na(p-1))/(nu(p)-nu(p-1))+na(p);
ns = nnu;
end

function ns = nsreddvvsr26(dv,npp,tem)
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
nc=trapz(dv,Ibn);
ns=abs(real(nc/trapz(dv,Ib)));
end