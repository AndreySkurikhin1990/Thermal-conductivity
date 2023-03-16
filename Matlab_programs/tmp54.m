%Только для вермикулита: определение средних значений степени черноты, показателя преломления и его
%квадрата, отражательной способности, коэффициента поглощения
function t = tmp54()
te0=273.15; te=23; %te=2e1:1:2e3; 
te=te+te0; te=te';
npp=Kramers_n_uk();
%npp=Kramers_n();
%dv=RasshDiapDlinVoln();
dv=dlinyvoln();
%fileID = fopen('Koefficient_pogloscheniya_ver.txt','r'); formatSpec='%f';
%als=fscanf(fileID,formatSpec); fclose(fileID);
%n2s=SreZnaPokPrel2(te,npp,dv); n2s=n2s'; %показатель преломления в квадрате
for k=0:3
%alpm = vyborka(dv,k);
end
%alp=podmassiv(dv,als,0);
%dl=podmassiv(dv,dv,1);
%np=podmassiv(dv,npp,0);
%n2s=SreZnaPokPrel2(te,np,dl)
%kp=SreZnaKoefPogl(te,npp,dv,als); %коэффициент поглощения
%skor=sredRosSiegel(te(1), np, dl, alp, n2s)
%ns=SreZnaPokPrel(te,npp,dv); %показатель преломления
%posp=SreZnaPogSpos(te,npp,dv); %поглощательная способность или степень черноты
%npp=PoiskSredAlpha(npp,dv);
t=0;
end

%определение среднего значения степени черноты
function szkp = SreZnaPogSpos(te,npp,dv)
eps=epsillam(npp);
for k =1:length(te)
ab(k)=usrednen(te(k),eps,npp,dv); 
end
ab=real(ab')
%p=plot(te,ab,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Степень черноты (поглощательная способность)'}); 
%title({'График зависимости степени черноты от температуры'});
szkp=0;
end

%определение среднего значения коэффициента поглощения
function szkp = SreZnaKoefPogl(te,npp,dv,als)
%als=RasMasKoAbs();
%als=1e6*SredGraf();
for k =1:length(te)
kp(k)=usrednen(te(k),als,npp,dv); 
end
kp=real(kp')
leAb=1e6/kp
%t=0;
%t=ZapisFileOptio(kp);
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости коэффициента поглощения от температуры'});
szkp=ZapisFileOptioPokPogl(dv,als,npp);
end

%определение среднего значения квадрата показателя преломления
function [ szkp ] = SreZnaPokPrel2(te,npp,dv)
for k=1:length(npp)
    npp2(k)=(npp(k))^2;
end
for k =1:length(te)
npp2s(k)=usrednen(te(k),npp2,npp,dv); 
end
npp2s=real(npp2s');
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости квадрата показателя преломления от температуры'});
szkp=npp2s;
end

%определение среднего значения квадрата показателя преломления
function szkp = SreZnaPokPrel(te,npp,dv)
for k =1:length(te)
npps(k)=usrednen(te(k),npp,npp,dv); 
end
npps=real(npps');
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости показателя преломления от температуры'});
szkp=0;
end

function uv = usrednen(T,usv,npp,dv)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*(c0^2);
c2=PP*c0/PB;
c1t=c1; c2t=c2;
for k=1:length(dv)
    np=npp(k);
    c1t=c1t/(np^2); c2t=c2t/np;
    la=dv(k)/np; dl(k)=la;
izz(k)=c1t/(la^5)/(exp(c2t/(la*T))-1);
izc(k)=usv(k)*izz(k); %усредняемая величина
c1t=c1;
c2t=c2;
end
izz=izz';
chi=integpo2ma(dl,izc);
zna=integpo2ma(dl,izz);
uv=chi/zna; %усредненная величина
end

function inte = integpo2ma(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end
%определение поглощательной способности по формуле Данкла
function [ epsla ] = epsillam(pp)
ep=0;
for k=1:length(pp)
    n=abs(pp(k));
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(n)/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log((n-1)/(n+1))*((n^2-1)^2)/((n^2+1)^3);
ep(k)=abs(eps);
eps=0;
end
epsla=real(ep);
end

function [ rs ] = sredRosSiegel(tem,npp,dl, alfs, n2)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
eps=1e-6;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4));
for k=1:length(npp)
c1m=C1/(npp(k)^2);
c2m=C2/npp(k);
Consk=(pi/2)*(c1m*c2m)/(sig*n2);
la=dl(k)/npp(k);
dlv(k)=la;
chi=exp(c2m/la/tem);
zna=(chi-1)^2;
Ibz(k)=Consk*chi/zna/(la^6)/(tem^5); 
if (abs(alfs(k))>eps)
Ibc(k)=Ibz(k)/alfs(k);
end
end
chasc=trapz(dlv,Ibc);
chasz=trapz(dlv,Ibz);
dlsvprfo2=chasc/chasz;
dlsvprfo2=1/dlsvprfo2;
rs=dlsvprfo2';
end

function [ np ] = vyborka(Sp,ident)
if (ident==0)
fileID=fopen('Koefficient_pogloscheniya_ver_60.txt','r'); 
elseif (ident==1)
fileID=fopen('Koefficient_pogloscheniya_ver_100.txt','r'); 
elseif (ident==2)
fileID=fopen('Koefficient_pogloscheniya_ver_150.txt','r'); 
elseif (ident==3)
fileID=fopen('Koefficient_pogloscheniya_ver_200.txt','r'); 
end
formatSpec='%f';
nnu=fscanf(fileID,formatSpec); 
fclose(fileID);
n=length(Sp); shag=1e0; nac=2; kon=27;
mr=nac:shag:kon; ko=1e-6;
mr=mr*ko;
nonu=0;
for k=1:length(mr)
    la=mr(k); f=1; q=2;
    for j=2:n
        if ((Sp(j)>la) && (f>0))
            q=j; 
            f=0; 
            break;
        end
    end
    if (f<1)
        pr=(nnu(q)-nnu(q-1))/(Sp(q)-Sp(q-1));
        nonu(k)=pr*(la-Sp(q-1))+nnu(q-1);
    end
end
nonu=nonu'
if (ident==0)
t=ZapisFileOptioVybKoefPogl60(nonu);
elseif (ident==1)
t=ZapisFileOptioVybKoefPogl60_100(nonu);
elseif (ident==2)
t=ZapisFileOptioVybKoefPogl100_150(nonu);
elseif (ident==3)
t=ZapisFileOptioVybKoefPogl150_200(nonu);
mr=mr'
end
np=[0];
end

function pn = PoiskNomera(dlvo,dvr)
n=length(dlvo); f=1; q=1;
for k=1:n
    if ((dlvo(k)>dvr) && (f>0))
        f=0; 
        q=k;
        break;
    end
end
pn=q;
end

function [ pm ] = podmassiv(dlv,arr,iden)
na=3*1e-6; ko=10*1e-6;
ary=0;
    w1=PoiskNomera(dlv,na);
    ary(1)=arr(w1-1)+(arr(w1)-arr(w1-1))*(na-dlv(w1-1))/(dlv(w1)-dlv(w1-1));
    w2=PoiskNomera(dlv,ko); 
    ary(2)=arr(w1);
    q=w1+1;
    w3=w2-w1+2;
    for k=3:w3
        ary(k)=arr(q); 
        q=q+1;
    end
    ary(w3)=arr(w2-1)+(arr(w2)-arr(w2-1))*(ko-dlv(w2-1))/(dlv(w2)-dlv(w2-1));
    if (iden==1)
        ary(1)=na;
        ary(w3)=ko;
    end
pm=ary;
end

function [ dlv ] = dlinyvoln()
format long g;
dl=0; dl=dlvoVer53101(); 
p=length(dl); 
fid = fopen('Dliny_voln_ver.txt','w');
for k=1:p  
    fprintf(fid,'%0.20f\n',dl(k));
    dl(k)=1e-2/dl(k);
end
fclose(fid);
dlv = dl;
end

function t = ZapisFileOptio(massi)
fid = fopen('Koefficient_pogloscheniya_ver_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioVybKoefPogl60(massi)
fid = fopen('Koefficient_pogloscheniya_Vyborka60.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioVybKoefPogl60_100(massi)
fid = fopen('Koefficient_pogloscheniya_Vyborka100.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioVybKoefPogl100_150(massi)
fid = fopen('Koefficient_pogloscheniya_Vyborka150.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioVybKoefPogl150_200(massi)
fid = fopen('Koefficient_pogloscheniya_Vyborka200.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioPokPogl(dv,als,npp)
for k =1:length(als)
    pp(k)=dv(k)*als(k)/4/pi/npp(k);
end
massi=pp';
fid = fopen('Pokazatel_pogloscheniya_ver.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioFractions(massi)
fid = fopen('Koefficient_pogloscheniya_Fractions.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioRoss(massi)
fid = fopen('Koefficient_pogloscheniya_Rossel.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function [ alsr ] = PoiskSredAlpha(PokPrel,dv)
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
format long g; schot=1;
%менее 60 мкм
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; Tal8=0; Tal9=0;
Tal=privedkEdiPropus(TrkoVer5311()); Tal2=privedkEdiPropus(TrkoVer5312());
Tal3=privedkEdiPropus(TrkoVer5313()); Tal4=privedkEdiPropus(TrkoVer53101()); 
Tal5=privedkEdiPropus(TrkoVer53102()); Tal6=privedkEdiPropus(TrkoVer53103()); 
mkbr=250.239; mv=0.283; tol=0.73; rokbr=2.75; rov=0.49; 
mkbr2=249.740; mv2=0.464; tol2=0.72; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); 
vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
mkbr3=249.223; mv3=0.812; tol3=0.72; 
mkbr4=250.395; mv4=0.22; tol4=0.72; rov3=0.49; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); 
vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
mkbr5=250.366; mv5=0.547; tol5=0.71; rokbr=2.75; rov5=0.49; 
mkbr6=249.55; mv6=0.777; tol6=0.7; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); 
vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
xv=(vv/(vv+vkbr))*tol*1e3; 
xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3;
xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; 
xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
nTal=length(Tal); xvt=0; xvt=[xv,xv2,xv3,xv4,xv5,xv6]';
for k=1:nTal 
    Tal(k)=-log(Tal(k))/xv; 
    Tal2(k)=-log(Tal2(k))/xv2;
    Tal3(k)=-log(Tal3(k))/xv3;
    Tal4(k)=-log(Tal4(k))/xv4;
    Tal5(k)=-log(Tal5(k))/xv5;
    Tal6(k)=-log(Tal6(k))/xv6; 
end
Tsr=23+273.15; koe=1e6; mn2=0;
for k=1:length(PokPrel)
mn2(k)=PokPrel(k)^2;
end
n2=usrednen(Tsr,mn2,PokPrel,dv)
alsrSredVer(schot)=usrednen(Tsr,koe*Tal,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal2,PokPrel,dv);
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal2, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal3,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal3, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal4,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal4, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal5,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal5, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal6,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal6, n2); schot=schot+1;
alsr1=0;
for k=1:nTal 
    xv=Tal(k)+Tal2(k)+Tal3(k);
    xv2=Tal4(k)+Tal5(k)+Tal6(k);
    alsr1(k)=xv+xv2;
end
alsr1=alsr1/length(xvt);
k=ZapisFile60(koe*alsr1);
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
%фракция от 60 до 100 мкм
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
vkbr8=mkbr8/(1e3*rokbr); vv8=mv8/(1e3*rov8);
mkbr9=249.206; mv9=0.843; tol9=0.76; rov9=0.52;
vkbr9=mkbr9/(1e3*rokbr); vv9=mv9/(1e3*rov9); 
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
xv7=(vv7/(vv7+vkbr7))*tol7*1e3; xv8=(vv8/(vv8+vkbr8))*tol8*1e3;
xv9=(vv9/(vv9+vkbr9))*tol9*1e3; xvt=[xv,xv2,xv3,xv4,xv5,xv6,xv7,xv8,xv9]';
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
alsrSredVer(schot)=usrednen(Tsr,koe*Tal,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal2,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal2, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal3,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal3, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal4,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal4, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal5,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal5, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal6,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal6, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal7,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal7, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal8,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal8, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal9,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal9, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsr2=0; s1=0; s2=0;
for k=1:nTal 
    s1=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k);
    s2=Tal6(k)+Tal7(k)+Tal8(k)+Tal9(k);
    alsr2(k)=s1+s2; s1=0; s2=0;
end
alsr2=alsr2/length(xvt);
k=ZapisFile60_100(koe*alsr2);
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
%фракция от 100 до 150 мкм
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; Tal8=0; Tal9=0;
Tal=privedkEdiPropus(TrkoVer5551()); Tal2=privedkEdiPropus(TrkoVer5552()); 
mkbr=249.913; mv=0.315; tol=0.74; rokbr=2.75; rov=0.53; 
mkbr2=249.607; mv2=0.473; tol2=0.74; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
Tal3=privedkEdiPropus(TrkoVer5553()); Tal4=privedkEdiPropus(TrkoVer5561());
mkbr3=249.218; mv3=0.709; tol3=0.72; 
mkbr4=249.929; mv4=0.293; tol4=0.72; rov3=0.53; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
Tal5=privedkEdiPropus(TrkoVer5562()); Tal6=privedkEdiPropus(TrkoVer5563());
mkbr5=249.695; mv5=0.528; tol5=0.71; rokbr=2.75; rov5=0.53; 
mkbr6=249.306; mv6=0.83; tol6=0.7; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); 
vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
Tal7=privedkEdiPropus(TrkoVer5571()); Tal8=privedkEdiPropus(TrkoVer5572());
mkbr7=250.405; mv7=0.27; tol7=0.78; rov7=0.53;
mkbr8=249.625; mv8=0.493; tol8=0.73; rov8=rov7;
vkbr7=mkbr7/(1e3*rokbr); vv7=mv7/(1e3*rov7); 
vkbr8=mkbr8/(1e3*rokbr); vv8=mv8/(1e3*rov8);
Tal9=privedkEdiPropus(TrkoVer5573());
mkbr9=249.348; mv9=0.764; tol9=0.76; rov9=0.53;
vkbr9=mkbr9/(1e3*rokbr); vv9=mv9/(1e3*rov9);
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
xv7=(vv7/(vv7+vkbr7))*tol7*1e3; xv8=(vv8/(vv8+vkbr8))*tol8*1e3;
xv9=(vv9/(vv9+vkbr9))*tol9*1e3; xvt=0; xvt=[xv,xv2,xv3,xv4,xv5,xv6,xv7,xv8,xv9]';
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
end
alsrSredVer(schot)=usrednen(Tsr,koe*Tal,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal2,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal2, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal3,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal3, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal4,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal4, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal5,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal5, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal6,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal6, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal7,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal7, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal8,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal8, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal9,PokPrel,dv); 
spr=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal9, n2); 
alSrRos(schot)=spr; schot=schot+1;
alsr3=0;
for k=1:nTal
    xv=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k);
    xv2=Tal6(k)+Tal7(k)+Tal8(k)+Tal9(k);
    alsr3(k)=xv+xv2;
end
alsr3=alsr3/length(xvt);
k=ZapisFile100_150(koe*alsr3);
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
%фракция от 150 до 200 мкм
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; Tal8=0; Tal9=0;
Tal=privedkEdiPropus(TrkoVer5681()); Tal2=privedkEdiPropus(TrkoVer5682());
mkbr=250.882; mv=0.320; tol=0.76; rokbr=2.75; rov=0.56; 
mkbr2=249.590; mv2=0.533; tol2=0.72; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
Tal3=privedkEdiPropus(TrkoVer5683()); Tal4=privedkEdiPropus(TrkoVer5691());
mkbr3=249.213; mv3=0.849; tol3=0.69; 
mkbr4=250.299; mv4=0.223; tol4=0.73; rov3=0.56; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
Tal5=privedkEdiPropus(TrkoVer5692()); Tal6=privedkEdiPropus(TrkoVer5693());
mkbr5=249.441; mv5=0.502; tol5=0.73; rokbr=2.75; rov5=0.56; 
mkbr6=249.365; mv6=0.797; tol6=0.73; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3;
xvt=0; xvt=[xv,xv2,xv3,xv4,xv5,xv6]';
for k=1:nTal
    Tal(k)=-log(Tal(k))/xv;
    Tal2(k)=-log(Tal2(k))/xv2;
    Tal3(k)=-log(Tal3(k))/xv3;
    Tal4(k)=-log(Tal4(k))/xv4;
    Tal5(k)=-log(Tal5(k))/xv5;
    Tal6(k)=-log(Tal6(k))/xv6; 
end
alsrSredVer(schot)=usrednen(Tsr,koe*Tal,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal2,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal2, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal3,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal3, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal4,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal4, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal5,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal5, n2); schot=schot+1;
alsrSredVer(schot)=usrednen(Tsr,koe*Tal6,PokPrel,dv); 
alSrRos(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal6, n2); schot=schot+1;
alsr4=0;
for k=1:nTal 
    xv=Tal(k)+Tal2(k)+Tal3(k);
    xv2=Tal4(k)+Tal5(k)+Tal6(k);
    alsr4(k)=xv+xv2;
end
alsr4=alsr4/length(xvt);
k=ZapisFile150_200(koe*alsr4);
alsre=0;
for k=1:nTal
    alsre(k)=(alsr1(k)+alsr2(k)+alsr3(k)+alsr4(k));
end
dlinaalsrSredVer=length(alsrSredVer');
dlinaalSrRos=length(alSrRos');
alsrSredVer=alsrSredVer';
alSrRos=alSrRos'
xv=ZapisFileOptioFractions(alsrSredVer);
xv=ZapisFileOptioRoss(alSrRos);
alsre=alsre/4;
k=ZapisFile150_200(koe*alsre);
alsr=alsre;
end

function [ prked ] = privedkEdiPropus(ar)
arr=ar;
p=length(arr);
for k=1:p
    if (arr(k)>1)
        arr(k)=1;
    end
    if (arr(k)<0)
        arr(k)=0;
    end
end
prked=arr;
end

function t = ZapisFile(massi)
fid = fopen('Koefficient_pogloscheniya_ver.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFile60(massi)
fid = fopen('Koefficient_pogloscheniya_ver_60.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFile60_100(massi)
fid = fopen('Koefficient_pogloscheniya_ver_100.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFile100_150(massi)
fid = fopen('Koefficient_pogloscheniya_ver_150.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFile150_200(massi)
fid = fopen('Koefficient_pogloscheniya_ver_200.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end