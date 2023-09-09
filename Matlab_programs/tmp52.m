%Только для ИТОМ и КВИ: определение средних значений степени черноты, показателя преломления и его
%квадрата, отражательной способности, коэффициента поглощения
function t = tmp52()
%te0=273.15; te=2e1:1e0:2e3; te=te+te0;
%dv=(dlvoItom1()+dlvoItom2())/2;
%d = ZapisFileOptiodv(dv);
%npp=Kramers_n_Itom();
%for k=1:length(dv)
    %dv(k)=1e-2/dv(k);
%end
%posp=SreZnaPogSpos(te,npp,dv); %поглощательная способность
%kp=SreZnaKoefPogl440(te); %коэффициент поглощения
%kp=SreZnaKoefPogl620(te); %коэффициент поглощения
%kp=SreZnaKoefPogl860(te); %коэффициент поглощения
%kp=SreZnaKoefPogl1000(te); %коэффициент поглощения
%n2s=SreZnaPokPrel2(te,npp,dv); %показатель преломления в квадрате
%ns=SreZnaPokPrel(te,npp,dv); %показатель преломления
%kp=SreZnaKoefPoglkvi400(te);
%kp=SreZnaKoefPoglkvi500(te);
%kp=SreZnaKoefPoglkvi600(te);
%kp=SreZnaKoefPoglkvi700(te);
%kp=SreZnaKoefPoglkvi800(te);
%kp=SreZnaKoefPoglkvi900(te);
%kp=SreZnaKoefPoglkvi1000(te);
%kp=RasshDiapDlinVoln();
%kp=ZapisFileOptiodvkvi(kp);
%kp=opredSredKPFrackvi();
for kp=0:3
    kp=kp'
    e(kp+1)=SreZnaKoefPoglitom620_1000(kp);
end
kp=mean(kp);
kp=opredSredKPFracitom();
t=0;
end

%определение среднего значения степени черноты
function szkp = SreZnaPogSpos(te,npp,dv)
eps=epsillam(npp);
for k =1:length(te)
ab(k)=usrednen(te(k),eps,dv,npp);
end
ab=real(ab');
%p=plot(te,ab,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Степень черноты (поглощательная способность)'}); 
%title({'График зависимости степени черноты от температуры'});
szkp=0;
end

%определение среднего значения коэффициента поглощения
function szkp = SreZnaKoefPogl440(te)
fileID=fopen('Koefficient_pogloscheniya_itom.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_Itom();
dv=dlinyvolnitom();
for k =1:length(te)
    kp(k)=usrednen(te(k),als,dv,npp);
end
kp=kp';
t=ZapisFileOptio440(kp);
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости коэффициента поглощения от температуры'});
szkp=t;
end

%определение среднего значения коэффициента поглощения
function szkp = SreZnaKoefPogl620(te)
fileID=fopen('Koefficient_pogloscheniya_itom620.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_Itom620();
dv=RasshDiapDlinVoln();
for k =1:length(te)
    kp(k)=usrednen(te(k),als,dv,npp);
end
kp=kp';
t=ZapisFileOptio620(kp);
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости коэффициента поглощения от температуры'});
szkp=t;
end

function szkp = SreZnaKoefPogl860(te)
fileID=fopen('Koefficient_pogloscheniya_itom860.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_Itom860();
dv=RasshDiapDlinVoln();
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
kp=kp';
t=ZapisFileOptio860(kp);
szkp=t;
end

function szkp = SreZnaKoefPogl1000(te)
fileID=fopen('Koefficient_pogloscheniya_itom1000.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_Itom1000();
dv=RasshDiapDlinVoln();
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
kp=kp';
t=ZapisFileOptio1000(kp);
szkp=t;
end

%определение среднего значения квадрата показателя преломления
function szkp = SreZnaPokPrel2(te,npp,dv)
for k=1:length(npp)
    npp2(k)=npp(k)*npp(k);
end
for k =1:length(te)
npp2s(k)=usrednen(te(k),npp2,dv,npp); 
end
npp2s=real(npp2s')
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости квадрата показателя преломления от температуры'});
szkp=0;
end

%определение среднего значения квадрата показателя преломления
function szkp = SreZnaPokPrel(te,npp,dv)
for k =1:length(te)
npps(k)=usrednen(te(k),npp,dv,npp); 
end
npps=real(npps')
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости показателя преломления от температуры'});
szkp=0;
end

function knus = usrednen(tem,knuSha,dv,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*(c0^2);
c2=PP*c0/PB;
ct1=0; ct2=0; dl=0;
for k=1:length(npp)
    ct1=ct1/(npp(k)^2);
    ct2=ct2/npp(k);
    lambda=dv(k)/npp(k);
    dl(k)=lambda;
    Ib(k)=ct1/((lambda^5)*(exp(ct2/(lambda*tem))-1));
    Ib(k)=ProvAd(Ib(k));
    Ibc(k)=knuSha(k)*Ib(k);
ct1=c1;
ct2=c2;
end
nc=trapz(dl,Ibc);
nz=trapz(dl,Ib);
knus=nc/nz;
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

function t = ZapisFileOptio440(massi)
fid = fopen('Koefficient_pogloscheniya_itom440_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptio620(massi)
fid = fopen('Koefficient_pogloscheniya_itom620_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptio860(massi)
fid = fopen('Koefficient_pogloscheniya_itom860_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptio1000(massi)
fid = fopen('Koefficient_pogloscheniya_itom1000_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiodv(massi)
fid = fopen('Dliny_voln_itom.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function [ dlv ] = dlinyvolnitom()
dl=(dlvoItom1()+dlvoItom2())/2;
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
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

function nk = opredn2(PokPrel,Tsr,dv)
mn2=0;
for k=1:length(PokPrel)
mn2(k)=PokPrel(k)^2;
end
nk=usrednen(Tsr,mn2,dv,PokPrel);
end

function t = opredSredKPFrackvi()
schot=1; koe=1e0; dv=RasshDiapDlinVoln();
Tsr=23+273.15; 
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; 
%
fileID = fopen('Koefficient_pogloscheniya_kvi400.txt','r'); formatSpec='%f';
Tal=fscanf(fileID,formatSpec); fclose(fileID);
PokPrel=0; PokPrel=Kramers_n_kvi400();
n2=opredn2(PokPrel,Tsr,dv);
alsrSredkvi(schot)=usrednen(Tsr,koe*Tal,dv,PokPrel); 
alSrRoskvi(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_kvi500.txt','r'); formatSpec='%f';
Tal2=fscanf(fileID,formatSpec); fclose(fileID);
PokPrel=0; PokPrel=Kramers_n_kvi500();
n2=opredn2(PokPrel,Tsr,dv);
alsrSredkvi(schot)=usrednen(Tsr,koe*Tal2,dv,PokPrel);
alSrRoskvi(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal2, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_kvi600.txt','r'); formatSpec='%f';
Tal3=fscanf(fileID,formatSpec); fclose(fileID);
PokPrel=0; PokPrel=Kramers_n_kvi600();
n2=opredn2(PokPrel,Tsr,dv);
alsrSredkvi(schot)=usrednen(Tsr,koe*Tal3,dv,PokPrel); 
alSrRoskvi(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal3, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_kvi700.txt','r'); formatSpec='%f';
Tal4=fscanf(fileID,formatSpec); fclose(fileID);
PokPrel=0; PokPrel=Kramers_n_kvi700();
n2=opredn2(PokPrel,Tsr,dv);
alsrSredkvi(schot)=usrednen(Tsr,koe*Tal4,dv,PokPrel); 
alSrRoskvi(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal4, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_kvi800.txt','r'); formatSpec='%f';
Tal5=fscanf(fileID,formatSpec); fclose(fileID);
PokPrel=0; PokPrel=Kramers_n_kvi800();
n2=opredn2(PokPrel,Tsr,dv);
alsrSredkvi(schot)=usrednen(Tsr,koe*Tal5,dv,PokPrel); 
alSrRoskvi(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal5, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_kvi900.txt','r'); formatSpec='%f';
Tal6=fscanf(fileID,formatSpec); fclose(fileID);
PokPrel=0; PokPrel=Kramers_n_kvi900();
n2=opredn2(PokPrel,Tsr,dv);
alsrSredkvi(schot)=usrednen(Tsr,koe*Tal6,dv,PokPrel); 
alSrRoskvi(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal6, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_kvi1000.txt','r'); formatSpec='%f';
Tal7=fscanf(fileID,formatSpec); fclose(fileID);
PokPrel=0; PokPrel=Kramers_n_kvi1000();
n2=opredn2(PokPrel,Tsr,dv);
alsrSredkvi(schot)=usrednen(Tsr,koe*Tal7,dv,PokPrel); 
alSrRoskvi(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal7, n2); schot=schot+1;
%
alsrSredkvi=alsrSredkvi'
alSrRoskvi=alSrRoskvi'
xv=ZapisFileOptioFractionskvi(alsrSredkvi);
xv=ZapisFileOptioRosskvi(alSrRoskvi);
t=xv;
end

function t = ZapisFileOptioFractionskvi(massi)
fid = fopen('Koefficient_pogloscheniya_frac_kvi.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioRosskvi(massi)
fid = fopen('Koefficient_pogloscheniya_Rossel_kvi.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = opredSredKPFracitom()
schot=1; koe=1e0; koef=1e0;
dv=RasshDiapDlinVoln();
Sp=dlinyvolnitom();
xv=0; Tsr=22+273.15; 
Tal=0; Tal2=0; Tal3=0; Tal4=0; 
%
Tal=SredGrafItom();
%Tal=preobMasitom_n(alph,koef*Sp,koef*dv);
PokPrel=Kramers_n_Itom();
%PokPrel=preobMasitom_n(PokPrel,koef*Sp,koef*dv);
n2=opredn2(PokPrel,Tsr,Sp);
alsrSreditom(schot)=usrednen(Tsr,koe*Tal,Sp,PokPrel); 
alSrRositom(schot)=sredRosSiegel(Tsr,PokPrel,Sp, koe*Tal, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_itom620.txt','r'); formatSpec='%f';
Tal2=fscanf(fileID,formatSpec); fclose(fileID);
%Tal2=preobMasitom_n(Tal2,koef*Sp,koef*dv);
PokPrel=Kramers_n_Itom620();
%PokPrel=preobMasitom_n(PokPrel,koef*Sp,koef*dv);
n2=opredn2(PokPrel,Tsr,dv);
alsrSreditom(schot)=usrednen(Tsr,koe*Tal2,dv,PokPrel);
alSrRositom(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal2, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_itom860.txt','r'); formatSpec='%f';
Tal3=fscanf(fileID,formatSpec); fclose(fileID);
%Tal3=preobMasitom_n(Tal3,koef*Sp,koef*dv);
PokPrel=Kramers_n_Itom860();
%PokPrel=preobMasitom_n(PokPrel,koef*Sp,koef*dv);
n2=opredn2(PokPrel,Tsr,dv);
alsrSreditom(schot)=usrednen(Tsr,koe*Tal3,dv,PokPrel); 
alSrRositom(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal3, n2); schot=schot+1;
%
fileID = fopen('Koefficient_pogloscheniya_itom1000.txt','r'); formatSpec='%f';
Tal4=fscanf(fileID,formatSpec); fclose(fileID);
%Tal4=preobMasitom_n(Tal4,koef*Sp,koef*dv);
PokPrel=Kramers_n_Itom1000();
%PokPrel=preobMasitom_n(PokPrel,koef*Sp,koef*dv);
n2=opredn2(PokPrel,Tsr,dv);
alsrSreditom(schot)=usrednen(Tsr,koe*Tal4,dv,PokPrel); 
alSrRositom(schot)=sredRosSiegel(Tsr,PokPrel,dv, koe*Tal4, n2); schot=schot+1;
%
Tsr=Tsr'
alsrSreditom=alsrSreditom'
alSrRositom=alSrRositom'
alSritom=mean(alsrSreditom)
alSrRositom=mean(alSrRositom)
xv=ZapisFileOptioFractionsitom(alsrSreditom);
xv=ZapisFileOptioRossitom(alSrRositom);
t=xv;
end

function t = ZapisFileOptioFractionsitom(massi)
fid = fopen('Koefficient_pogloscheniya_frac_itom.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptioRossitom(massi)
fid = fopen('Koefficient_pogloscheniya_Rossel_itom.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

%-----------------------

%определение среднего значения степени черноты
function szkp = SreZnaPogSposkvi(te,npp,dv)
eps=epsillam(npp);
for k =1:length(te)
ab(k)=usrednen(te(k),eps,dv,npp); 
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
function szkp = SreZnaKoefPoglkvi400(te)
fileID=fopen('Koefficient_pogloscheniya_kvi400.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_kvi400();
%dv=dlinyvolnkvi();
dv=RasshDiapDlinVoln(); 
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
t=ZapisFileOptiokvi400(kp);
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости коэффициента поглощения от температуры'});
szkp=t;
end

%определение среднего значения коэффициента поглощения
function szkp = SreZnaKoefPoglkvi500(te)
fileID=fopen('Koefficient_pogloscheniya_kvi500.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_kvi500();
dv=RasshDiapDlinVoln();
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
t=ZapisFileOptiokvi500(kp);
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости коэффициента поглощения от температуры'});
szkp=t;
end

function szkp = SreZnaKoefPoglkvi600(te)
fileID=fopen('Koefficient_pogloscheniya_kvi600.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_kvi600();
dv=RasshDiapDlinVoln();
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
t=ZapisFileOptiokvi600(kp);
szkp=t;
end

function szkp = SreZnaKoefPoglkvi700(te)
fileID=fopen('Koefficient_pogloscheniya_kvi700.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_kvi700();
dv=RasshDiapDlinVoln();
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
t=ZapisFileOptiokvi700(kp);
szkp=t;
end

function szkp = SreZnaKoefPoglkvi800(te)
fileID=fopen('Koefficient_pogloscheniya_kvi800.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_kvi800();
dv=RasshDiapDlinVoln();
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
t=ZapisFileOptiokvi800(kp);
szkp=t;
end

function szkp = SreZnaKoefPoglkvi900(te)
fileID=fopen('Koefficient_pogloscheniya_kvi900.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_kvi900();
dv=RasshDiapDlinVoln();
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
t=ZapisFileOptiokvi900(kp);
szkp=t;
end

function szkp = SreZnaKoefPoglkvi1000(te)
fileID=fopen('Koefficient_pogloscheniya_kvi1000.txt','r'); 
formatSpec='%f';
als=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_n_kvi1000();
dv=RasshDiapDlinVoln();
for k =1:length(te)
kp(k)=usrednen(te(k),als,dv,npp); 
end
t=ZapisFileOptiokvi1000(kp);
szkp=t;
end

%определение среднего значения квадрата показателя преломления
function szkp = SreZnaPokPrel2kvi(te,npp,dv)
for k=1:length(npp)
    npp2(k)=npp(k)*npp(k);
end
for k =1:length(te)
npp2s(k)=usrednen(te(k),npp2,dv,npp); 
end
npp2s=real(npp2s')
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости квадрата показателя преломления от температуры'});
szkp=0;
end

%определение среднего значения квадрата показателя преломления
function szkp = SreZnaPokPrelkvi(te,npp,dv)
for k =1:length(te)
npps(k)=usrednen(te(k),npp,dv,npp); 
end
npps=real(npps')
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости показателя преломления от температуры'});
szkp=0;
end

function t = ZapisFileOptiokvi400(massi)
fid = fopen('Koefficient_pogloscheniya_kvi400_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiokvi500(massi)
fid = fopen('Koefficient_pogloscheniya_kvi500_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiokvi600(massi)
fid = fopen('Koefficient_pogloscheniya_kvi600_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiokvi700(massi)
fid = fopen('Koefficient_pogloscheniya_kvi700_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiokvi800(massi)
fid = fopen('Koefficient_pogloscheniya_kvi800_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiokvi900(massi)
fid = fopen('Koefficient_pogloscheniya_kvi900_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiokvi1000(massi)
fid = fopen('Koefficient_pogloscheniya_kvi1000_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiodvkvi(massi)
fid = fopen('Dliny_voln_kvi.txt','w');
p=length(massi);
for k=1:p
    massi(k)=1e-2/massi(k);
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function [ dlv ] = dlinyvolnkvi()
dl=(dlvoItom1()+dlvoItom2())/2;
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end

function [ pre ] = preobMasitom_n(alp2,dlv2,dlv1)
format longg; eps=1e-10;
lealp=length(alp2);
ledlv2=length(dlv2);
ledlv1=length(dlv1);
%dlv1=preobDlinVoln(ledlv1,ledlv2,dlv1,dlv2,1);
%dlv1=dlvoVer53101(); ledlv1=length(dlv1);
for p=1:ledlv1
    a=dlv1(p);
    fa=-1;
    a1=0;
    b1=0;
    f1=0;
    f2=0;
for k=2:ledlv2
        if ((dlv2(k-1)<=a) && (dlv2(k)>a) && (fa<0))
            a1=dlv2(k-1);
            b1=dlv2(k);
            f1=alp2(k-1);
            f2=alp2(k);
            fa=1;
        end
        if (fa>0)
            break;
        end
end
if (fa>0)
    if (abs(b1-a1)>eps)
    ko=(f2-f1)/(b1-a1);
    alp1(p)=f1+(a-a1)*ko;
    end
else
    alp1(p)=0;
end
end
a1=dlv2(1);            
b1=dlv2(2);
f1=alp2(1);
f2=alp2(2);
alp1(1)=alp2(1)+(f2-f1)*(dlv1(1)-a1)/(b1-a1);
n=ledlv1;
m=ledlv2;
a1=dlv2(m-1);            
b1=dlv2(m);
f1=alp2(m-1);
f2=alp2(m);
alp1(n)=f2+(f2-f1)*(dlv1(n)-b1)/(b1-a1);
alp1(n-1)=(alp1(n)+alp1(n-2))/2;
for k=1:n
    alp1(k)=ProvAd(alp1(k));
end
pre=alp1;
end

function [ w ] = preobDlinVoln(ledlv1,ledlv2,dlv1,dlv2,v)
for k=1:ledlv1
    dlv1(k)=1e-2/dlv1(k);
end
for k=1:ledlv2
    dlv2(k)=1e-2/dlv2(k);
end
if (v==1)
    w=dlv1;
else
    w=dlv2;
end
end

function mk = ProvAd(m)
ep=1e-30;
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
if (abs(m)<ep)  
    m=0;  
end; 
mk=real(abs(m));
end

function szkp = SreZnaKoefPoglitom440()
tena=2e2; teko=9e2; te0=273.15; dt=1e2;
tem=tena:dt:teko; tem=tem+te0; koe=1e0;
n=length(tem);
dv=RasshDiapDlinVoln();
Sp=dlinyvolnitom();
Tal=SredGrafItom();
%Tal=preobMasitom_n(alph,koef*Sp,koef*dv);
PokPrel=Kramers_n_Itom();
%PokPrel=preobMasitom_n(PokPrel,koef*Sp,koef*dv);
for k=1:n
Tsr=tem(k);
n2=opredn2(PokPrel,Tsr,Sp);
alsrSredVer(k)=usrednen(Tsr,koe*Tal,dv,PokPrel); 
alSrRos(k)=sredRosSiegel(Tsr,PokPrel,Sp, koe*Tal, n2);
end
disp('itom-440');
alsrSredVer=alsrSredVer'
alSrRos=alSrRos'
szkp=0;
end

function szkp = SreZnaKoefPoglitom620_1000(ide)
tena=2e2; teko=9e2; 
te0=273.15; dt=1e2; tem=tena:dt:teko; 
%tem=23;
tem=tem+te0; koe=1e0;
n=length(tem);
Sp=RasshDiapDlinVoln();
%Tal=preobMasitom_n(alph,koef*Sp,koef*dv);
if (ide==0)
PokPrel=Kramers_n_Itom();
fileID=fopen('Koefficient_pogloscheniya_itom440.txt','r'); 
elseif (ide==1)
PokPrel=Kramers_n_Itom620();
fileID=fopen('Koefficient_pogloscheniya_itom620.txt','r'); 
elseif (ide==2)
PokPrel=Kramers_n_Itom860();
fileID=fopen('Koefficient_pogloscheniya_itom860.txt','r'); 
elseif (ide==3)
PokPrel=Kramers_n_Itom1000();
fileID=fopen('Koefficient_pogloscheniya_itom1000.txt','r'); 
end
%PokPrel=preobMasitom_n(PokPrel,koef*Sp,koef*dv);
formatSpec='%f';
Tal=fscanf(fileID,formatSpec); 
fclose(fileID);
for k=1:n
Tsr=tem(k);
n2=opredn2(PokPrel,Tsr,Sp);
alsrSreditom(k)=usrednen(Tsr,koe*Tal,Sp,PokPrel); 
alSrRositom(k)=sredRosSiegel(Tsr,PokPrel,Sp, koe*Tal, n2);
end
tem=tem'
n2=n2'
alsrSreditom=alsrSreditom'
alSrRositom=alSrRositom'
szkp=0;
end