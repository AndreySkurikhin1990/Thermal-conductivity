function t = tmp5()
%определение средних показателя преломления и коэффициента отражения шамота и вермикулита
format long g;
%alfs=0; alfs=SredGraf(); 
%knuSha=preobMas();
%knuSha=preobMasSha_n();
%npp=0; npp=Kramers_n(); 
knuSha=1e6*SredGrafSha();
dl=dlina_vol();
%ronu=0; nsh=0; nsha=0; nv=0; hi=0; hiSha=0; np=0;  nsham=0;
nsh=Kramers_n_Sham(knuSha,dl); p=length(dl);
for k=1:p    
    %hiSha(k)=dl(k)*knuSha(k)/(4*pi); 
    %hi(k)=dl(k)*alfs(k)/(4*pi);
    %nsha(k)=((nsh(k)^2+hiSha(k)^2)/(npp(k)^2+hi(k)^2))^0.5;
    %nsha(k)=sqrt(nsh(k)^2+hiSha(k)^2);
    nsha(k)=abs(nsh(k));
    ronu(k)=0.5+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); 
    ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));
    ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log((nsha(k)-1)/(nsha(k)+1))/((nsha(k)^2+1)^3);
end
nsha=nsha';
ronu=ronu';
tem=1e2:1e2:12e2; tem=tem+273.15;
for k=1:length(tem)
ro(k) = usrednen(tem(k),ronu,nsh,dl);
end
ro=abs(ro')
t=0;
end

function [ d ] = dlina_vol()
dv=(dlvoSham1()+dlvoSham2())/2;
for k=1:length(dv)
    dv(k)=1e-2/dv(k);
end
d=dv;
end

function uv = usrednen(T,usv,npp,dv)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*c0^2;
c2=PP*c0/PB;
iz0=0; iz1=0; dl=0;
c1t=c1; c2t=c2;
for j=1:length(dv)
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
iz0(j)=c1t/((la^5)*(exp(c2t/(la*T))-1));
iz1(j)=usv(j)*iz0(j); %усредняемая величина
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
uv=chi/zna; %усредненная величина
end

function t = tmp5_1()
p=length(dl); 
su=0; 
for j=1:p-1 
    su=su+(dl(j+1)-dl(j))*(nsha(j)+nsha(j+1))/2; 
end
nsham=real(su/(dl(p)-dl(1))); 
nsham=real(nsham)
su=0; 
for j=1:p-1 
    su=su+(dl(j+1)-dl(j))*(ronu(j)+ronu(j+1))/2; 
end
Refl=su/(dl(p)-dl(1)); 
Refl=real(Refl)
su=0; 
for j=1:p-1 
    su=su+(dl(j+1)-dl(j))*(npp(j)+npp(j+1))/2; 
end
np=su/(dl(p)-dl(1)); 
np=real(np)
nsham=trapz(dl,nsha)/(dl(p)-dl(1)); 
nsham=real(nsham)
Refl=trapz(dl,ronu)/(dl(p)-dl(1)); 
Refl=real(Refl)
%np=trapz(dl,npp)/(dl(p)-dl(1)); np=real(np)
t=0;
end

%определяет показатель преломления шамота
function [ PokPrel ] = Kramers_n_Sham(alsr,Sp)
format long g;
c0=299792458;
%SpS=dlvoSham1();
%alsrS=1e6*SredGrafSha();
%leSp=length(SpS)
%for k=1:leSp
	%SpS(k)=1e-2/SpS(k);
    %nuS(k)=2*pi*c0/SpS(k);
%end
%Sp=dlivoln();
%Sp=dliny_voln();
leSp=length(Sp);
%alsr=preobMasSha_n();
for k=1:leSp
    nu(k)=2*pi*c0/Sp(k);
    hi(k)=Sp(k)*alsr(k)/4/pi;
end
enu=1e-3;
ma=izmMasChast(nu,enu); %определение массива nus на enu от шага nu
na=PokazPrelomAl(nu,ma,alsr,c0); %главная часть - определение ПП по КК
nnu=PokazPrelomPribl(nu,ma,na); %определение второго массива ПП
%ep0=opredEpsilSt(nu,ma,hi,alsr)
%PokPrel=PokazPrelomProver(nu,ma,hi,nnu,na,alsr);
%PokPrel=PoiskPokPre(Sp,nnu+0.6);
PokPrel=nnu-1+1.56;
end

function [ dlv ] = dliny_voln()
%vl=299792458; 
format long g;
dl=0; dl=dlvoVer53101(); 
%dl=RasshDiapDlinVoln();
p=length(dl); 
%fid = fopen('Dlina_volny.txt','w');
for k=1:p  
    dl(k)=1e-2/dl(k); %fr(k)=vl/dl(k); 
    %fprintf(fid,'%0.20f\n',dl(k));
end
%fclose(fid);
dlv = dl;
end

function pg = postGraf(lb,ns,la,pp)
%c0=299792458;
%leSp=length(la);
%for k=1:leSp
	%lb(k)=2*pi*c0/lb(k);
	%la(k)=2*pi*c0/la(k);
%end
%pl=plot(lb*1e6,ns,'-b');
pl=plot(lb*1e6,ns,'-b',la*1e6,pp,'-k');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'});
%ylabel({'Показатель преломления'}); 
%title({'График зависимости ПП от длины волны'});
%xlabel({'Показатель преломления'}); 
ylabel({'Степень черноты'}); 
title({'График зависимости СЧ от ДВ'});
pg=0;
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
%for k=1:lenu
%nust(k)=(nu(k)+nu(k+1))/2;
%end
%nust(lenu)=(nust(lenu-1)+nu(lenu))/2;
%nus=nust;
end

function [ pop ] = PokazPrelomAl(om,mao,ka,vl)
np=0;
p=length(om);
for k=1:p
    q=1;
    for j=1:p-1
        dkpodo=(ka(j+1)-ka(j))/(om(j+1)-om(j));
        podln=(om(j)+mao(k))/(om(j)-mao(k));
        podln=abs(podln);
        fn(q)=dkpodo*log(podln);
        q=q+1;
    end
    fn(p)=fn(p-1);
    np(k)=1+(vl/pi)*integpo2ma(om,fn)/2/mao(k);
end
pop=np;
end

function pn = PoiskPokPre(dv,ep)
%Routschka
sc = [0.468 0.472 0.4667 0.5 0.5171 0.5132 0.5 0.6 0.6811 0.7 ...
0.8 0.9 0.9675 0.9727 0.9779 0.9753 0.9909 0.9857 0.9260 ...
0.9553 0.9675 0.9571 0.9545 0.9519 0.9 0.8447];
la=[1 1.24 2 2.71 2.86 3 3.16 3.77 4 4.09 4.47 4.85 5.61 ...
6 6.19 7 7.83 8 9 10 11 12 13 14 14.59 15]*1e-6;
lela=length(la);
c0=299792458;
enu=1e-3;
for k=1:lela
    om(k)=2*pi*c0/la(k);
end
%Шамот
%Zhuravlyov, Blokh 293 K Average
sc1 = [0.225 0.19 0.205 0.4 0.54 0.425 0.43 0.625 0.825 ...
    0.875 0.94 0.915];
la1 = [1 1.5 2 2.5 3 3.5 4 4.5 5 6 7 8]*1e-6;
%Zhuravlyov, Blokh 800-1500 K
sc2 = [0.26 0.225 0.24 0.37 0.525 0.455 0.5 0.675 0.8 ...
    0.85 0.91 0.915];
%Zhuravlyov, Blokh 293 K max
sc4 = [0.4 0.3 0.31 0.5 0.75 0.45 0.48 0.7 0.9 0.95 0.98 0.98];
%ma=izmMasChast(om,enu);
%na=PokazPrelomAl(om,ma,alsr,c0); 
%nom=PokazPrelomPribl(om,ma,na);
na=poikor(sc); %na=na'
%nb=poikor(sc1); nb=1*nb'
%nc=poikor(sc2); nc=1*nc'
%ne=poikor(sc4); ne=ne'
%q=1;
%for k=1:length(na)
    %if (na(k)>1)
    %naa(q)=na(k);
    %q=q+1;
    %else
        %naa(q)=1;
        %q=q+1;
    %end
%end
%es=epsiAver(la,sc,300,naa);
for k=1:length(ep)
eps(k)=epsilam(ep(k));
end
%eps=1*eps'
pn=postGraf(la,sc,dv,eps);
end

function epsi = epsilam(n)
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(abs(n))/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log(abs(n-1)/abs(n+1))*((n^2-1)^2)/((n^2+1)^3);
epsi=eps;
end

function es = epsiAver(dv,eps,T,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*c0^2;
c2=PP*c0/PB;
iz0=0; iz1=0; dl=0;
c1t=c1; c2t=c2;
for j=1:length(dv)
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
iz0(j)=c1t/((la^5)*(exp(c2t/(la*T))-1));
iz1(j)=eps(j)*iz0(j);
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
es=chi/zna;
end

function [ ko ] = poikor(izsp)
Nit=2e3;
eps=1e-12;
for k=1:length(izsp)
na=1e-3;
nb=1e3; 
t=1; sc=izsp(k);
do=abs(nb-na);
while (do>eps)
nc=(na+nb)/2;
ffc=epsilam(nc)-sc;
ffa=epsilam(na)-sc;
ffb=epsilam(nb)-sc;
if (ffa>0)
        if (ffb>0)
    nb=1e3;
    na=nc;
    %ffa=1*ffa
    %ffb=1*ffb
    %ffc=1*ffc
        end
end
if (ffa<0)
        if (ffb<0)
    na=1e-3;
    nb=nc;
    %ffa=1*ffa
    %ffb=1*ffb
    %ffc=1*ffc
        end
end
if (ffc*ffb>0)
    if (ffc*ffa<0)
nb=nc;    
    end
end
if (ffc*ffa>0)
    if (ffc*ffb<0)
na=nc;
    end
end
t=t+1;
if (t>Nit)
    break; 
end
do=abs(nb-na);
end
n(k)=nc;
end
ko=n;
end

%--------------
%--------------
%--------------

function [ npopr ] = PokazPrelomProver(nu,nus,hi,nnu,na,als)
hia=PokazPoglPribl(nus,als,nu);
p=length(hi);
for j=1:p
nmk(j)=(na(j)^2)-(hia(j)^2);
end
for k=1:p
    ff=0;
    tt=2*na(k)*hia(k)*nus(k);
    for j=1:p
        t=(nu(j)^2)-(nus(k)^2);
        ff(j)=(2*nnu(j)*hi(j)*nu(j)-tt)/t;
    end
    nmkp(k)=1+(2/pi)*integpo2ma(nu,ff);
end
npopr=postGraf(nu,nus,nmk,nmkp);
end

function n0 = opredn0(dl,ka)
n0=1+(1/2/pi^2)*integpo2ma(dl,ka);
end

function epsi = opredEpsilSt(nu,nus,hinu,alsr)
hia=PokazPoglPribl(nus,alsr,nu);
dl=dlinyvoln();
n0=opredn0(dl,alsr)
leSp=length(nu);
eps=1e-4;
epst=0;
for k=1:leSp
euta=1e-3;
eutb=1e3; 
Nit=2e2;
p=hia(k);
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
ma=-(2*nus(k)/pi)*integpo2ma(nu,ffa);
ma=ProvAdekv(ma);
ma=ma-p;
mb=-(2*nus(k)/pi)*integpo2ma(nu,ffb);
mb=ProvAdekv(mb);
mb=mb-p;
mc=-(2*nus(k)/pi)*integpo2ma(nu,ffc);
mc=ProvAdekv(mc);
mc=mc-p;
if (ma*mc>0)        
euta=eutc;    
end
if (mc*mb>0)
eutb=eutc;
end
t=t+1;
if (t>Nit)
    break; 
end
end
epst(k)=eutc;
end
epsi=eutc;
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

function [ ns ]  = PokazPrelomPribl(nu,nus,na)
nnu=0;
p=length(nu);
for k=1:p-1
    de=(nus(k)-nu(k))*(na(k+1)-na(k))/(nu(k+1)-nu(k));
    nnu(k)=de+na(k);
end
nnu(p)=(nus(p)-nu(p))*(na(p)-na(p-1))/(nu(p)-nu(p-1))+na(p);
ns = nnu;
end

function [ hia ] = PokazPoglPribl(nus,alsr,nu)
c0=299792458;
p=length(nus);
als=0;
his=0;
for k=1:p-1
    de=(nus(k)-nu(k))*(alsr(k+1)-alsr(k))/(nu(k+1)-nu(k));
    als(k)=alsr(k)+de;
end
als(p)=(nus(p)-nu(p))*(alsr(p)-alsr(p-1))/(nu(p)-nu(p-1))+alsr(p);
for k=1:p
dls(k)=2*pi*c0/nus(k);
his(k)=dls(k)*als(k)/4/pi;
end
hia=his;
end

function [ dl ] = dlinyvoln()
Sp=dlvoSham1();
leSp=length(Sp);
for k=1:leSp
	Sp(k)=1e-2/Sp(k);
end
dl=Sp;
end