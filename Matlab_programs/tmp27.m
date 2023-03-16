function tm = tmp27()
format long g;
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; Tal8=0; Tal9=0;
Tal=privedkEdiPropus(TrkoVer5421()); Tal2=privedkEdiPropus(TrkoVer5422()); 
mkbr=250; mv=0.255; tol=0.72; rokbr=2.75; rov=0.52; 
mkbr2=249.629; mv2=0.539; tol2=0.71; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
Tal3=privedkEdiPropus(TrkoVer5423()); Tal4=privedkEdiPropus(TrkoVer5431());
mkbr3=249.294; mv3=0.809; tol3=0.7;  
mkbr4=249.706; mv4=0.295; tol4=0.7; rov3=0.52; rov4=rov3;
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
Tal5=privedkEdiPropus(TrkoVer5432()); Tal6=privedkEdiPropus(TrkoVer5433());
mkbr5=249.510; mv5=0.517; tol5=0.73; rokbr=2.75; rov5=0.52; 
mkbr6=249.307; mv6=0.756; tol6=0.72; rov6=rov5;
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
Tal7=privedkEdiPropus(TrkoVer5441()); Tal8=privedkEdiPropus(TrkoVer5442());
Tal9=privedkEdiPropus(TrkoVer5443());
mkbr7=250.328; mv7=0.36; tol7=0.74; rov7=0.52;
mkbr8=249.604; mv8=0.534; tol8=0.7; rov8=rov7;
vkbr7=mkbr7/(1e3*rokbr); vv7=mv7/(1e3*rov7); 
vkbr8=mkbr8/(1e3*rokbr); vv8=mv4/(1e3*rov8);
mkbr9=249.206; mv9=0.843; tol9=0.76; rov9=0.52;
vkbr9=mkbr9/(1e3*rokbr); vv9=mv9/(1e3*rov9); 
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
xv7=(vv7/(vv7+vkbr7))*tol7*1e3; xv8=(vv8/(vv8+vkbr8))*tol8*1e3;
xv9=(vv9/(vv9+vkbr9))*tol9*1e3;
nTal=length(Tal);
Tal=abs(Tal);
Tal2=abs(Tal2);
Tal3=abs(Tal3);
Tal4=abs(Tal4);
Tal5=abs(Tal5);
Tal6=abs(Tal6);
Tal7=abs(Tal7);
Tal8=abs(Tal8);
Tal9=abs(Tal9);
for k=1:nTal
    Tal(k)=-log(Tal(k))/xv;
    Tal2(k)=-log(Tal2(k))/xv2;
    Tal3(k)=-log(Tal3(k))/xv3;
    Tal4(k)=-log(Tal4(k))/xv4;
    Tal5(k)=-log(Tal5(k))/xv5;
    Tal6(k)=-log(Tal6(k))/xv6;
    Tal7(k)=-log(Tal7(k))/xv7;
    Tal8(k)=-log(Tal8(k))/xv8;
    Tal9(k)=-log(Tal9(k))/xv9; 
end;
alsr2=0;
for k=1:nTal 
    alsr2(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k);
    alsr2(k)=alsr2(k)+Tal6(k)+Tal7(k)+Tal8(k)+Tal9(k);
end
alsr2=alsr2/9;
alsr2=arrachec(alsr2);
Tal=arrachec(Tal);
Tal2=arrachec(Tal2);
Tal3=arrachec(Tal3);
Tal4=arrachec(Tal4);
Tal5=arrachec(Tal5);
Tal6=arrachec(Tal6);
Tal7=arrachec(Tal7);
Tal8=arrachec(Tal8);
Tal9=arrachec(Tal9);
alfs1 = 0; alfs1 = RasMasKoAbs27(1e6*Tal);
alfs2 = 0; alfs2 = RasMasKoAbs27(1e6*Tal2);
alfs3 = 0; alfs3 = RasMasKoAbs27(1e6*Tal3);
alfs4 = 0; alfs4 = RasMasKoAbs27(1e6*Tal4);
alfs5 = 0; alfs5 = RasMasKoAbs27(1e6*Tal5);
alfs6 = 0; alfs6 = RasMasKoAbs27(1e6*Tal6);
alfs7 = 0; alfs7 = RasMasKoAbs27(1e6*Tal7);
alfs8 = 0; alfs8 = RasMasKoAbs27(1e6*Tal8);
alfs9 = 0; alfs9 = RasMasKoAbs27(1e6*Tal9);
alfs1=arrachec(alfs1);
alfs2=arrachec(alfs2);
alfs3=arrachec(alfs3);
alfs4=arrachec(alfs4);
alfs5=arrachec(alfs5);
alfs6=arrachec(alfs6);
alfs7=arrachec(alfs7);
alfs8=arrachec(alfs8);
alfs9=arrachec(alfs9);
    pp1=0; pp1=Kramers_n27(alfs1);
    pp2=0; pp2=Kramers_n27(alfs2);
    pp3=0; pp3=Kramers_n27(alfs3);
    pp4=0; pp4=Kramers_n27(alfs4);
    pp5=0; pp5=Kramers_n27(alfs5);
    pp6=0; pp6=Kramers_n27(alfs6);
    pp7=0; pp7=Kramers_n27(alfs7);
    pp8=0; pp8=Kramers_n27(alfs8);
    pp9=0; pp9=Kramers_n27(alfs9);
    dl=0; dl=RasshDiapDlinVoln27();
    pp=0; pp=Kramers_n27(alsr2);
    pp1=arrachec(pp1);
    pp2=arrachec(pp2);
    pp3=arrachec(pp3);
    pp4=arrachec(pp4);
    pp5=arrachec(pp5);
    pp6=arrachec(pp6);
    pp7=arrachec(pp7);
    pp8=arrachec(pp8);
    pp9=arrachec(pp9);
    pp=arrachec(pp);
    dl=arrachec(dl);
    %wl = 1e6*RasshDiapDlinVoln26();
T = 3e2:1e3:4e2;
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0; al7=0; al8=0; al9=0;
for k=1:length(T)
al1(k) = 1e6 / sredRosSieg27(T(k),alfs1,dl,pp1);
al2(k) = 1e6 / sredRosSieg27(T(k),alfs2,dl,pp2);
al3(k) = 1e6 / sredRosSieg27(T(k),alfs3,dl,pp3);
al4(k) = 1e6 / sredRosSieg27(T(k),alfs4,dl,pp4);
al5(k) = 1e6 / sredRosSieg27(T(k),alfs5,dl,pp5);
al6(k) = 1e6 / sredRosSieg27(T(k),alfs6,dl,pp6);
al7(k) = 1e6 / sredRosSieg27(T(k),alfs7,dl,pp7);
al8(k) = 1e6 / sredRosSieg27(T(k),alfs8,dl,pp8);
al9(k) = 1e6 / sredRosSieg27(T(k),alfs9,dl,pp9);
alp(k) = 1e6 / sredRosSieg(T(k));
end
al1=1*al1'
al2=1*al2'
al3=1*al3'
al4=1*al4'
al5=1*al5'
al6=1*al6'
al7=1*al7'
al8=1*al8'
al9=1*al9'
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0; al7=0; al8=0; al9=0;
for k=1:length(T)
al1(k)=AlphaSredDVVSr27(dl,pp1,alfs1,T(k));
al2(k)=AlphaSredDVVSr27(dl,pp2,alfs2,T(k));
al3(k)=AlphaSredDVVSr27(dl,pp3,alfs3,T(k));
al4(k)=AlphaSredDVVSr27(dl,pp4,alfs4,T(k));
al5(k)=AlphaSredDVVSr27(dl,pp5,alfs5,T(k));
al6(k)=AlphaSredDVVSr27(dl,pp6,alfs6,T(k));
al7(k)=AlphaSredDVVSr27(dl,pp7,alfs7,T(k));
al8(k)=AlphaSredDVVSr27(dl,pp8,alfs8,T(k));
al9(k)=AlphaSredDVVSr27(dl,pp9,alfs9,T(k));
end
al1=1*al1'
al2=1*al2'
al3=1*al3'
al4=1*al4'
al5=1*al5'
al6=1*al6'
al7=1*al7'
al8=1*al8'
al9=1*al9'
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0; al7=0; al8=0; al9=0;
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
popr7=SuzhDiapDanPoPr(1e6*dl,pp7);
kop7=SuzhDiapDanKoPo(1e6*dl,alfs7);
popr8=SuzhDiapDanPoPr(1e6*dl,pp8);
kop8=SuzhDiapDanKoPo(1e6*dl,alfs8);
popr9=SuzhDiapDanPoPr(1e6*dl,pp9);
kop9=SuzhDiapDanKoPo(1e6*dl,alfs9);
for k=1:length(T)
al1(k) = 1e6 / sredRosSieg27(T(k),kop1,dv,popr1);
al2(k) = 1e6 / sredRosSieg27(T(k),kop2,dv,popr2);
al3(k) = 1e6 / sredRosSieg27(T(k),kop3,dv,popr3);
al4(k) = 1e6 / sredRosSieg27(T(k),kop4,dv,popr4);
al5(k) = 1e6 / sredRosSieg27(T(k),kop5,dv,popr5);
al6(k) = 1e6 / sredRosSieg27(T(k),kop6,dv,popr6);
al7(k) = 1e6 / sredRosSieg27(T(k),kop7,dv,popr7);
al8(k) = 1e6 / sredRosSieg27(T(k),kop8,dv,popr8);
al9(k) = 1e6 / sredRosSieg27(T(k),kop9,dv,popr9);
end
al1=1*al1'
al2=1*al2'
al3=1*al3'
al4=1*al4'
al5=1*al5'
al6=1*al6'
al7=1*al7'
al8=1*al8'
al9=1*al9'
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0; al7=0; al8=0; al9=0;
for k=1:length(T)
al1(k)=AlphaSredDVVSr27(dv,popr1,kop1,T(k));
al2(k)=AlphaSredDVVSr27(dv,popr2,kop2,T(k));
al3(k)=AlphaSredDVVSr27(dv,popr3,kop3,T(k));
al4(k)=AlphaSredDVVSr27(dv,popr4,kop4,T(k));
al5(k)=AlphaSredDVVSr27(dv,popr5,kop5,T(k));
al6(k)=AlphaSredDVVSr27(dv,popr6,kop6,T(k));
al7(k)=AlphaSredDVVSr27(dv,popr7,kop7,T(k));
al8(k)=AlphaSredDVVSr27(dv,popr8,kop8,T(k));
al9(k)=AlphaSredDVVSr27(dv,popr9,kop9,T(k));
end
al1=1*al1'
al2=1*al2'
al3=1*al3'
al4=1*al4'
al5=1*al5'
al6=1*al6'
al7=1*al7'
al8=1*al8'
al9=1*al9'
al1=0; al2=0; al3=0; al4=0; al5=0; al6=0; al7=0; al8=0; al9=0;
tm = alp';
%tm = KoPogTabl27(1e6*dlivoln(), Tal);
end

function [ kp ] = KoPogTabl27(wl, alsr)
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

function [ rs ] = sredRosSieg27(tem,alfs,dl,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4));
np=nsreddvvsr27(dl,npp,tem);
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
Ibc=arrachec(Ibc);
Ibz=arrachec(Ibz);
chasc=trapz(dlv,Ibc);
chasz=trapz(dlv,Ibz);
dlsvprfo2=chasc/chasz;
rs=1e6*dlsvprfo2;
end

function [ ndv ] = RasshDiapDlinVoln27()
Sp=skleiv2Mas27(dlinvolndoMas27(),dlinvoln27());
ndv=Sp;
end

function [ dlv ] = dlinvoln27()
dl=0; 
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end

function [ dlv ] = dlinvolndoMas27()
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

function [ obma ] = skleiv2Mas27(ar1,ar2)
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

function [ rasmaskopo ] = RasMasKoAbs27(alph)
alsr=skleiv2Mas27(soglMas27(alph,1e6*SredGraSt27()),alph);
alsr=arrachec(alsr);
rasmaskopo=alsr;
end

function [ dlv ] = soglMas27(ar1,ar2)
dl=0; 
dl=dlvoi(); 
dlvv=0;
dlvv=dlvoVer53101();
kon=dlvv(1);
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

function [ sred ] = SredGraSt27()
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
alsr1=arrachec(alsr1);
sred=alsr1;
end

function ns = AlphaSredDVVSr27(dv,npp,alph,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
vl=c0;
Ib=0; Ibn=0;
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
Ibn=arrachec(Ibn);
Ib=arrachec(Ib);
nc=trapz(dv,Ibn);
nz=trapz(dv,Ib);
nz=abs(real(nc/nz));
ns=nz;
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

function [ PokPrel ] = Kramers_n27(alp)
c0=299792458;
enu=1e-3;
alsr=skleiv2Mas27(soglMas27(alp,1e6*SredGraSt27()),alp);
alsr=arrachec(alsr);
Sp=RasshDiapDlinVoln27();
leSp=length(Sp);
nu=0; nus=0; hi=0;
for k=1:leSp
    hi(k)=Sp(k)*alsr(k)/4/pi;
    nu(k)=2*pi*c0/Sp(k);
end
nus=izmMasChast27(nu,enu);
nnus=PokazPrelomAl27(nu,nus,alsr,c0)-1+1.563;
nnu=PokazPrelomPribl27(nu,nus,nnus);
nnu=arrachec(nnu);
PokPrel=nnu;
end

function inte = integpo2ma27(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ nus ] = izmMasChast27(nu,enu)
nust=0;
lenu=length(nu);
for k=1:lenu-1
    nust(k)=enu*(nu(k+1)-nu(k))+nu(k);
end
nust(lenu)=2*nust(lenu-1)-nust(lenu-2);
nus=nust;
end

function [ pop ] = PokazPrelomAl27(nu,nus,ka,vl)
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
    np(k)=1+(vl/pi)*integpo2ma27(nu,fn)/2/nu(k);
end
pop=np;
end

function [ ns ]  = PokazPrelomPribl27(nu,nus,na)
nnu=0;
p=length(nu);
for k=1:p-1
    de=(-nus(k)+nu(k))*(na(k+1)-na(k))/(nu(k+1)-nu(k));
    nnu(k)=de+na(k);
end
nnu(p)=(-nus(p)+nu(p))*(na(p)-na(p-1))/(nu(p)-nu(p-1))+na(p);
ns = nnu;
end

function ns = nsreddvvsr27(dv,npp,tem)
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
Ibn=arrachec(Ibn);
Ib=arrachec(Ib);
nc=trapz(dv,Ibn);
nz=trapz(dv,Ib);
nz=abs(real(nc/nz));
ns=nz;
end

function [ arache ] = arrachec(ar)
for k=1:length(ar)
    if (isnan(ar(k))>0)
        ar(k)=0;
    end
    if (isinf(ar(k))>0)
        ar(k)=0;
    end    
end
arache=ar;
end