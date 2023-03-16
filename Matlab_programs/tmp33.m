%программа рассчитывает ЛКТП от температуры (или координаты). Исходные
%значения ЛПТП даны в функциях самого скрипта, аппроксимируются параболой
function tm = tmp33()
format long g;
%t=SkhodPosled800();
%t=SkhodPosled1000();
t=RasTem585();
tm=t;
end

function tm = RasTem585()
na=NachToch();
ko=KonechToch();
Nt=[3e1,6e1,12e1,24e1];
tol = 30e-3;
Tk = 392.9;
T0 = 858.15;
Tsred = 650.15;
te = [T0 Tsred Tk];
x = [0 tol/2 tol];
koe = polyfit(x,te,length(x)-1);
qq=1;
for k = 1:length(Nt)
    kt = Nt(k);
    ksi=0:tol/kt:tol;
    for j=1:length(ksi)
    tem(j) = polyval(koe,ksi(j));
    end
ppie=poluchMassivPTP(k,na,ko,ksi);
Tz=poluchMassivTemp(na,ko,ksi,tem);
xi=poluchMassivKoor(na,ko,ksi);
q=0; 
for j=1:length(ppie)
    dmq(qq,j)=ppie(j);
    dmt(qq,j)=Tz(j);
    dmxr(qq,j)=xi(j);
    q=q+1;
end
dlma(qq)=q;
qq=qq+1;
end
kt=PriblizhkBolChiToc(dmt,dmq,dlma,dmxr,Nt);
tm=0;
end

function [ tm ] = RasTem800(v)
na=NachToch();
ko=KonechToch();
Ntn=[3e1,6e1,12e1,24e1,48e1];
Nt=[3e1,6e1,24e1,48e1];
tol=30e-3;
q=1;
na=1;
for k = na:length(Nt)
    kt = Nt(k);
    ksi=0:tol/kt:tol;
ppie=poluchMassivPTP800(k,na,ko,ksi,kt);
dm(q)=length(ppie);
for j=1:dm(q)
    dmq(q,j)=ppie(j);
end
q=q+1;
end
if (v==0)
tm=dmq;
elseif (v==1)
    tm=dm;
    elseif (v==2)
    tm=Nt;
end
end

function [ tm ] = RasTem800_t()
na=NachToch();
ko=KonechToch();
Ntn=[3e1 6e1 12e1 24e1 48e1];
Nt=[3e1,6e1,24e1,48e1];
Tk = 473.15;
T0 = 1073.15;
Tsred = 821.15;
te = [T0 Tsred Tk];
tol=30e-3;
x = [0 tol/2 tol];
koe = polyfit(x,te,length(x)-1);
q=1;
na=1;
for k = na:length(Nt)
    kt = Nt(k);
    ksi=0:tol/kt:tol;
    for j=1:length(ksi)
    tem(j) = polyval(koe,ksi(j));
    end
Tz=poluchMassivTemp(na,ko,ksi,tem);
for j=1:length(Tz)
    dmt(q,j)=Tz(j);
end
q=q+1;
end
tm=dmt;
end

function [ tm ] = RasTem800_x()
na=NachToch();
ko=KonechToch();
tol=30e-3;
Ntn=[3e1 6e1 12e1 24e1 48e1];
Nt=[3e1,6e1,24e1,48e1];
q=1;
na=1;
for k = na:length(Nt)
    kt = Nt(k);
    ksi=0:tol/kt:tol;
xi=poluchMassivKoor(na,ko,ksi);
for j=1:length(xi)
    dmxr(q,j)=xi(j);
end
q=q+1;
end
tm=dmxr;
end

function [ tm ] = RasTem1000(v,Nt)
na=NachToch();
ko=KonechToch();
%Ntn=[3e1,45,6e1,1e2,12e1,15e1,2e2];
%Nt=[3e1,45,15e1,20e1];
tol=30e-3; 
q=1;
nac=1;
for k = nac:length(Nt)
    kt = Nt(k);
    ksi=0:tol/kt:tol;
ppie=poluchMassivPTP1000(k,na,ko,ksi,kt);
dm(q)=length(ppie);
for j=1:length(ppie)
    dmq(q,j)=ppie(j);
end
q=q+1;
end
if (v==0)
tm=dmq;
elseif (v==1)
    tm=dm;
end
end

function [ tm ] = RasTem1000_new()
na=NachToch1000();
ko=KonechToch1000();
Ntn=[3e1,45,6e1,1e2,12e1,15e1,2e2,24e1];
Nt=[3e1,45,15e1,20e1];
tol = 30e-3;
Tk = 546.15;
T0 = 1273.15;
Tsred = 976.15;
te = [T0,Tsred,Tk];
x = [0,tol/2,tol];
p=opredPorPribl();
koe = polyfit(x,te,p);
qq=1;
nac=1;
for k = nac:length(Nt)
    kt = Nt(k);
    ksi=0;
    ksi=0:tol/kt:tol;
    for j=1:length(ksi)
    tem(j) = polyval(koe,ksi(j));
    end
ppie=poluchMassivPTP1000(k,na,ko,ksi,kt);
Tz=poluchMassivTemp(na,ko,ksi,tem);
xi=poluchMassivKoor(na,ko,ksi);
q=0; 
for j=1:length(ppie)
    dmq(qq,j)=ppie(j);
    dmt(qq,j)=Tz(j);
    dmxr(qq,j)=xi(j);
    q=q+1;
end
dlma(qq)=q;
qq=qq+1;
end
kt=PriblizhkBolChiToc1000(dmt,dmq,dlma,dmxr,Nt);
tm=kt;
end

function [ tm ] = RasTem1000_t(Nt)
na=NachToch();
ko=KonechToch();
%Ntn=[3e1,45,6e1,1e2,12e1,15e1,2e2];
%Nt=[3e1,45,15e1,20e1];
Tk = 546.15;
T0 = 1273.15;
Tsred = 976.15;
te = [T0 Tsred Tk];
tol=30e-3;
x = [0 tol/2 tol];
koe = polyfit(x,te,length(x)-1);
q=1;
nac=1;
for k = nac:length(Nt)
    kt = Nt(k);
    ksi=0:tol/kt:tol;
    for j=1:length(ksi)
    tem(j) = polyval(koe,ksi(j));
    end
Tz=poluchMassivTemp(na,ko,ksi,tem);
for j=1:length(Tz)
    dmt(q,j)=Tz(j);
end
q=q+1;
end
tm=dmt;
end

function [ tm ] = RasTem1000_x(Nt)
na=NachToch();
ko=KonechToch();
tol=30e-3;
%Ntn=[3e1,45,6e1,1e2,12e1,15e1,2e2];
q=1;
nac=1;
for k = nac:length(Nt)
    kt = Nt(k);
    ksi=0:tol/kt:tol;
xi=poluchMassivKoor(na,ko,ksi);
for j=1:length(xi)
    dmxr(q,j)=xi(j);
end
q=q+1;
end
tm=dmxr;
end

function [ ar ] = poluchMassivPTP800(no,nac,kon,x,kt)
qr=0;
switch (kt)
        case 30
        qr = tmp3330_800;
        case 60
        qr = tmp3360_800;
        case 120
        qr = tmp33120_800;
        case 240
        qr = tmp33240_800;
        case 480
        qr = tmp33480_800;
end
nk=PoiskNacKon(x,nac,kon);
qrr=VydelPodmas(qr,nk(1),nk(2));
ar=qrr;
end

function [ ar ] = poluchMassivPTP1000(no,nac,kon,x,kt)
qr=0;
switch (kt)
        case 30
        qr = tmp3330_1000();
        case 45
        qr = tmp3345_1000();
        case 60
        qr = tmp3360_1000();
        case 100
        qr = tmp33100_1000();
        case 120
        qr = tmp33120_1000();
        case 150
        qr = tmp33150_1000();
        case 200
        qr = tmp33200_1000();
        case 240 
        qr = tmp33240_1000();
end
nk1000=PoiskNacKon(x,nac,kon);
qrr=VydelPodmas(qr,nk1000(1),nk1000(2));
ar=qrr;
end

function pg = postrGrafik(xi,y1,y2,y3,y30,ksi)
%te0=273.15;
te0=0;
%disp('tem eff 600');
%nn=vyvodZnach(xi);
%ktp4=MasLamEff(xi);
%disp('lambda eff 600');
%nn=vyvodZnach(ktp4);
ti2=opredPolTem(2,ksi);
%disp('tem eff 800');
%nn=vyvodZnach(ti2);
ktp5=MasLamEff(ti2);
%disp('lambda eff 800');
%nn=vyvodZnach(ktp5);
ti3=opredPolTem(3,ksi);
%disp('tem eff 1000');
%nn=vyvodZnach(ti3);
ktp6=MasLamEff(ti3);
%disp('lambda eff 1000');
%nn=vyvodZnach(ktp6);
%kt1=opredKoefTeplIzPlotTeplPotN(1,xi,y1)
%kt2=opredKoefTeplIzPlotTeplPotN(2,xi,y2)
%kt3=opredKoefTeplIzPlotTeplPotN(3,xi,y3)
%kt1=opredKoefTeplIzPlotTeplPot(xi,y1,ksi);
%kt2=opredKoefTeplIzPlotTeplPot(xi,y2,ksi);
%kt3=opredKoefTeplIzPlotTeplPot(xi,y3,ksi);
ktp1=opredKoefTeplIzPlotTeplPot1000(xi,y1,ksi);
%disp('lambda luch 600');
%nn=vyvodZnach(ktp1);
ktp2=opredKoefTeplIzPlotTeplPot1000(ti2,y2,ksi);
%disp('lambda luch 800');
%nn=vyvodZnach(ktp2);
ktp3=opredKoefTeplIzPlotTeplPot1000(ti3,y3,ksi);
%disp('lambda luch 1000');
%nn=vyvodZnach(ktp3);
ktp30=opredKoefTeplIzPlotTeplPot1000(ti3,y30,ksi);
%disp('lambda prib luch 1000');
%nn=vyvodZnach(ktp30);
%xk=opredMasSredTemp(xi);
tk=opredMasSredTemp1000(xi);
%xk=xi(2:length(xi))-te0
%xk=ksi(2:length(xi));
subplot(2,1,2);
tektpe=skleimas(tk,ti2,ti3);
%disp('tektpe');
%nn=vyvodZnach(tektpe);
ktpe=skleimas(ktp4,ktp5,ktp6);
%disp('ktpe');
%nn=vyvodZnach(ktpe);
nn=2; nx=length(ti2); ny=length(ktp2); nz=length(tektpe);
pl=plot(tk(nn:ny),ktp1(nn:ny),'-b',ti2(nn:nx-1),ktp2(nn:ny),'-k',ti3(nn:nx-1),ktp3(nn:ny),'-m',ti3(nn:nx-1),ktp30(nn:ny),'-r',tektpe(nn:nz-nn),ktpe(nn:nz-nn),'-g');
%pl=plot(ksi(nn:n),ktp1,'-b',ksi(nn:n),ktp2,'-k',ksi(nn:n),ktp3,'-m',ksi(nn:n),ktp4(nn:n),'-r',ksi(nn:n),ktp5(nn:n),'-c',ksi(nn:n),ktp6(nn:n),'-g');
set(pl,'LineWidth',2); 
hold on; grid on;
%xlabel({'Координата, м'}); 
xlabel({'Температура, К'}); 
%ylabel({'Температура, К'}); 
%xlabel({'Температура, °C'}); 
%ylabel({'Плотность теплового потока излучательной энергии'}); 
ylabel({'Коэффициент лучистой теплопроводности, Вт/(м*К)'}); 
title({'График зависимости ЛКТП от температуры'});
%title({'График зависимости ЛКТП от координаты'});
legend('626 К','773 К приб.','910 К','910 К приб.','ЭКТП');
%legend('ЛКТП (626 К)','ЭКТП');
%legend('626 К','773 К','910 К');
legend('location','best')
pg=0;
end

function pg = postrGrafik1000(tei,y,ksi)
%te0=273.15;
te0=0;
ktp2=MasLamEff(tei);
disp('tem eff 1000');
nn=vyvodZnach(tei);
disp('lambda eff 1000');
nn=vyvodZnach(ktp2);
ktp=opredKoefTeplIzPlotTeplPot1000(tei,y,ksi);
tek=opredMasSredTemp1000(tei);
disp('lambda luch 1000');
nn=vyvodZnach(ktp);
disp('tem luch 1000');
nn=vyvodZnach(tek);
nx=length(tek);
ny=length(ksi);
nn=2;
subplot(2,1,1);
pl=plot(tek(nn:nx),ktp(nn:nx),'-b',tek,ktp2(nn:ny),'-k');
%pl=plot(ksi(nn:n),ktp,'-b',ksi,ktp2,'-k');
set(pl,'LineWidth',2); 
hold on; grid on; 
xlabel({'Температура, К'}); 
%xlabel({'Координата, м'}); 
ylabel({'Коэффициент лучистой теплопроводности, Вт/(м*К)'}); 
title({'График зависимости ЛКТП от температуры'});
%title({'График зависимости ЛКТП от координаты'});
legend('ЛКТП (910 К)','ЭКТП');
hold off;
legend('location','best')
pg=0;
end

function vy = vyvodZnach(arr)
n=length(arr); ar=0; q=1;
for k=1:n
    if (rem(k,100)==0)
        ar(q)=arr(k);
        q=q+1;
    end
end
    ar=ar'
vy=0;
end

function [ le ] = MasLamEff(tem)
Te=arrTem_VVF2();
te0 = 273.15;
Te = Te + te0;
laef=arrKTP_VVF2();
p=opredPorPribl();
kola=polyfit(Te,laef,p);
q=1;
for k=1:length(tem)
    la(q)=polyval(kola,tem(k));
    q=q+1;
end
le=la;
end

function [ pr ] = PriblizhkBolChiToc(dmt,dmq,dlma,dmksi,Nt)
n=dlma(1);
for k=1:n
    ksi1(k)=dmksi(1,k);
    te1(k)=dmt(1,k); 
    ptp1(k)=dmq(1,k); 
end
n=dlma(2);
for k=1:n
    ksi2(k)=dmksi(2,k);
    te2(k)=dmt(2,k); 
    ptp2(k)=dmq(2,k); 
end
n=dlma(3);
for k=1:n
    ksi3(k)=dmksi(3,k);
    te3(k)=dmt(3,k); 
    ptp3(k)=dmq(3,k); 
end
n=dlma(4);
for k=1:n
    ksi4(k)=dmksi(4,k);
    te4(k)=dmt(4,k); 
    ptp4(k)=dmq(4,k); 
end
po=opredPorPribl();
ko1=polyfit(ksi1,ptp1,po);
ko2=polyfit(ksi2,ptp2,po);
ko3=polyfit(ksi3,ptp3,po);
ko4=polyfit(ksi4,ptp4,po);
h=opredShagKoo();
ksik=PoiskMin(ksi4):h:PoiskMax(ksi4);%te - температура, ksi - координата
%disp('nach koo = ');
%n=vyvodZnach(ksik);
koksi=polyfit(ksi4,te4,po);
n=length(ksik);
for k=1:n
tn(k)=polyval(koksi,ksik(k));
ptpn1(k)=polyval(ko1,ksik(k));
ptpn2(k)=polyval(ko2,ksik(k));
ot1(k)=ptpn2(k)/ptpn1(k); %60/30
ptpn3(k)=polyval(ko3,ksik(k));
ot2(k)=ptpn3(k)/ptpn1(k); %120/30
otn1(k)=(ptpn3(k)-ptpn2(k))/(ptpn2(k)-ptpn1(k));
ptpn4(k)=polyval(ko4,ksik(k));
ot3(k)=ptpn4(k)/ptpn1(k); %240/30
otn2(k)=(ptpn4(k)-ptpn3(k))/(ptpn3(k)-ptpn2(k));
end
prib=PoiskYnBeskon600(otn1,otn2,ot1,ptpn1,ptpn3,ptpn4,6e2,Nt);
prib1=PoiskYnBeskon800();
prib20=RasTem1000_new();
%disp('prib 1000');
pri=0;
%pri=vyvodZnach(prib20);
prib2=PoiskYnBeskon1000();
%for k=1:length(prib)
    %ot4(k)=prib(k)/ptpn1(k);
%end
%pri=postrGrafik(tn,prib,prib1,prib2,prib20,ksik);
%n=vyvodZnach(ksin);
pr = [ pri ];
end

function [ pr ] = PriblizhkBolChiToc1000(dmt,dmq,dlma,dmksi,Nt)
n=dlma(1);
for k=1:n
    ksi1(k)=dmksi(1,k); 
    te1(k)=dmt(1,k); 
    ptp1(k)=dmq(1,k); 
end
n=dlma(2);
for k=1:n
    ksi2(k)=dmksi(2,k); 
    te2(k)=dmt(2,k); 
    ptp2(k)=dmq(2,k); 
end
n=dlma(3);
for k=1:n
    ksi3(k)=dmksi(3,k); 
    te3(k)=dmt(3,k); 
    ptp3(k)=dmq(3,k); 
end
n=dlma(4);
for k=1:n
    ksi4(k)=dmksi(4,k);
    te4(k)=dmt(4,k); 
    ptp4(k)=dmq(4,k); 
end
po=opredPorPribl();
ko1=polyfit(ksi1,ptp1,po);
ko2=polyfit(ksi2,ptp2,po);
ko3=polyfit(ksi3,ptp3,po);
ko4=polyfit(ksi4,ptp4,po);
h=opredShagKoo();
ksin=PoiskMin(ksi4):h:PoiskMax(ksi4);%xn - температура
koksi=polyfit(ksi4,te4,po);
n=length(ksin);
for k=1:n
ten(k)=polyval(koksi,ksin(k));
ptpn1(k)=polyval(ko1,ksin(k));
ptpn2(k)=polyval(ko2,ksin(k));
ot1(k)=ptpn2(k)/ptpn1(k);
ptpn3(k)=polyval(ko3,ksin(k));
ot2(k)=ptpn3(k)/ptpn1(k);
otn1(k)=(ptpn3(k)-ptpn2(k))/(ptpn2(k)-ptpn1(k));
ptpn4(k)=polyval(ko4,ksin(k));
ot3(k)=ptpn4(k)/ptpn1(k);
otn2(k)=(ptpn4(k)-ptpn3(k))/(ptpn3(k)-ptpn2(k));
end
prib=PoiskYnBeskon600(otn1,otn2,ot1,ptpn1,ptpn3,ptpn4,1e3);
for k=1:length(prib)
    ot4(k)=ProvAd(prib(k)/ptpn1(k));
end
pri=postrGrafikOtn(ot1,ot2,ot3,ot4,1e3,Nt);
%pri=postrGrafik1000(ten,prib,ksin);
%n=vyvodZnach(ksin);
pr = prib;
end

function [ na ] = PoiskYnBeskon600(o1,o2,ot1,y1,y3,y4,te,Nt)
Nit=1e2; y30=y3; y40=y4;
f=1;
q=1;
ep=1e-3;
for j=1:Nit
    if (f>0)
for k=1:length(o1)
    o3(k)=ProvAd(o2(k)^2/o1(k));
    y5(k)=o3(k)*(y4(k)-y3(k))+y4(k);
    y5(k)=ProvAd(y5(k));
end
maxi=PoiskMax((y4-y5));
q=1;
if (maxi<ep)
    q=j;
    f=0;
end
o1=o2;
o2=o3;
y3=y4;
y4=y5;
    end
end
q=q';
for k=1:length(y30)
    ot3(k)=y30(k)/y1(k);
    ot3(k)=ProvAd(ot3(k));
end
for k=1:length(y40)
    ot4(k)=y40(k)/y1(k);
    ot4(k)=ProvAd(ot4(k));
end
for k=1:length(y5)
    ot5(k)=y5(k)/y1(k);
    ot5(k)=ProvAd(ot5(k));
end
if (te==6e2)
pri=postrGrafikOtn(ot1,ot3,ot4,ot5,te,Nt);
end
na = y5;
end

function [ na ] = PoiskYnBeskon800()
dmq1=RasTem800(0);
dmt1=RasTem800_t();
dmx1=RasTem800_x();
dm=RasTem800(1);
Nt=RasTem800(2);
disp(length(dm));
po=opredPorPribl();
n=dm(1); %30
for k=1:n
    ksi1(k)=dmx1(1,k); 
    te1(k)=dmt1(1,k); 
    ptp1(k)=dmq1(1,k); 
end
n=dm(2); %60
for k=1:n
    ksi2(k)=dmx1(2,k); 
    te2(k)=dmt1(2,k); 
    ptp2(k)=dmq1(2,k); 
end
n=dm(3); %120
for k=1:n
    ksi3(k)=dmx1(3,k); 
    te3(k)=dmt1(3,k); 
    ptp3(k)=dmq1(3,k); 
end
n=dm(4); %240
for k=1:n
    ksi4(k)=dmx1(4,k); 
    te4(k)=dmt1(4,k); 
    ptp4(k)=dmq1(4,k); 
end
if (length(dm)>4)
n=dm(5); %480
for k=1:n
    ksi5(k)=dmx1(5,k); 
    te5(k)=dmt1(5,k); 
    ptp5(k)=dmq1(5,k); 
end
ko5=polyfit(ksi5,ptp5,po);
end
ko1=polyfit(ksi1,ptp1,po);
ko2=polyfit(ksi2,ptp2,po);
ko3=polyfit(ksi3,ptp3,po);
ko4=polyfit(ksi4,ptp4,po);
h=opredShagKoo();
ksik=PoiskMin(ksi4):h:PoiskMax(ksi4);
koksi=polyfit(ksi4,te4,po);
n=length(ksik); ot1=0; ot2=0; ot3=0; ot4=0;
for k=1:n
tn(k)=polyval(koksi,ksik(k));
ptpn1(k)=polyval(ko1,ksik(k));
ptpn2(k)=polyval(ko2,ksik(k));
ot1(k)=ProvAd(ptpn2(k)/ptpn1(k)); %60/30
ptpn3(k)=polyval(ko3,ksik(k));
ot2(k)=ProvAd(ptpn3(k)/ptpn1(k)); %120/30
o1(k)=ProvAd((ptpn3(k)-ptpn2(k))/(ptpn2(k)-ptpn1(k)));
ptpn4(k)=polyval(ko4,ksik(k));
ot3(k)=ProvAd(ptpn4(k)/ptpn1(k)); %240/30
o2(k)=ProvAd((ptpn4(k)-ptpn3(k))/(ptpn3(k)-ptpn2(k)));
end
if (length(dm)>4)
    for k=1:n
        ptpn5(k)=polyval(ko5,ksik(k));
        ot4(k)=ProvAd(ptpn5(k)/ptpn1(k)); %480/30
        o3(k)=ProvAd((ptpn5(k)-ptpn4(k))/(ptpn4(k)-ptpn3(k)));
    end
end
%o1=o2;
%o2=o3;
%maxi=PoiskMax((ptpn4-ptpn5));
%ptpn3=ptpn4;
%ptpn4=ptpn5;
Nit=1e2;
f=1;
q=1;
ep=1e-3;
for j=1:Nit
    if (f>0)
        ptpn5=0; o3=0;
for k=1:length(o1)
    o3(k)=(o2(k)^2)/o1(k);
    o3(k)=ProvAd(o3(k));
    ptpn5(k)=o3(k)*(ptpn4(k)-ptpn3(k))+ptpn4(k);
    ptpn5(k)=ProvAd(ptpn5(k));
end
maxi=PoiskMax((ptpn4-ptpn5));
if (maxi<ep)
    q=j;
    f=0;
end
o1=o2;
o2=o3;
ptpn3=ptpn4;
ptpn4=ptpn5;
    end
end
for k=1:length(ptpn5)
    ot5(k)=ptpn5(k)/ptpn1(k);
    ot5(k)=ProvAd(ot5(k));
end
%pri=postrGrafikOtn800(ot1,ot2,ot3,ot4,ot5,8e2);
pri=postrGrafikOtn(ot1,ot2,ot3,ot5,8e2,Nt);
q=q';
na = ptpn5;
end

function mk = ProvAd(m)
e=1e-30; ms=1e10;
if (isnan(m))   
    m=0; 
end 
if (isinf(m))  
    m=0;  
end
if (abs(m)<e)
    m=0;
end
if (abs(m)>ms)
    m=0;
end
mk=real(abs(m));
end

function [ nat ] = PoiskYnBeskon1000()
Ntn=[3e1,45,6e1,1e2,12e1,15e1,2e2,24e1]; mi=1e2;
nNtn=length(Ntn);
for j1=2:nNtn-3
    for j2=j1+1:nNtn-2
        for j3=j2+1:nNtn-1
            for j4=j3+1:nNtn
Nt=[Ntn(j1),Ntn(j2),Ntn(j3),Ntn(j4)];
dmq2=RasTem1000(0,Nt);
dmt2=RasTem1000_t(Nt);
dmx2=RasTem1000_x(Nt);
dm=RasTem1000(1,Nt);
n=dm(1);
for k=1:n
    ksi1(k)=dmx2(1,k); 
    te1(k)=dmt2(1,k); 
    ptp1(k)=dmq2(1,k); 
end
n=dm(2);
for k=1:n
    ksi2(k)=dmx2(2,k); 
    te2(k)=dmt2(2,k); 
    ptp2(k)=dmq2(2,k); 
end
n=dm(3);
for k=1:n
    ksi3(k)=dmx2(3,k); 
    te3(k)=dmt2(3,k); 
    ptp3(k)=dmq2(3,k); 
end
n=dm(4);
for k=1:n
    ksi4(k)=dmx2(4,k); 
    te4(k)=dmt2(4,k);
    ptp4(k)=dmq2(4,k); 
end
po=opredPorPribl();
ko1=polyfit(ksi1,ptp1,po);
ko2=polyfit(ksi2,ptp2,po);
ko3=polyfit(ksi3,ptp3,po);
ko4=polyfit(ksi4,ptp4,po);
h=opredShagKoo();
ksik=PoiskMin(ksi4):h:PoiskMax(ksi4);%te - температура, ksi - координата%disp('nach koo = ');%n=vyvodZnach(ksik);
koksi=polyfit(ksi4,te4,po);
n=length(ksik);
for k=1:n
tn(k)=polyval(koksi,ksik(k));
ptpn1(k)=polyval(ko1,ksik(k));
ptpn2(k)=polyval(ko2,ksik(k));
ot1(k)=ptpn2(k)/ptpn1(k); %60/30
ptpn3(k)=polyval(ko3,ksik(k));
ot2(k)=ptpn3(k)/ptpn1(k); %120/30
o1(k)=(ptpn3(k)-ptpn2(k))/(ptpn2(k)-ptpn1(k));
ptpn4(k)=polyval(ko4,ksik(k));
ot3(k)=ptpn4(k)/ptpn1(k); %240/30
o2(k)=(ptpn4(k)-ptpn3(k))/(ptpn3(k)-ptpn2(k));
end
maxi=PoiskMax((ptpn4-ptpn3));
Nit=1e2;
f=1;
q=1;
ep=1e-3;
for j=1:Nit
    if (f>0)
        ptpn5=0; o3=0;
for k=1:length(o1)
    o3(k)=ProvAd(o2(k)^2/o1(k));
    ptpn5(k)=o3(k)*(ptpn4(k)-ptpn3(k))+ptpn4(k);
    ptpn5(k)=ProvAd(ptpn5(k));
end
maxi=PoiskMax((ptpn4-ptpn5));
if (maxi<ep)
    q=j;
    f=0;
end
o1=o2;
o2=o3;
ptpn3=ptpn4;
ptpn4=ptpn5;
    end
end
for k=1:length(ptpn5)
    ot4(k)=ProvAd(ptpn5(k)/ptpn1(k));
end
q=q';
q=postrGrafikOtn(ot1,ot2,ot3,ot4,1e3,Nt);
if (q<mi)
    mi=q;
    Ntk=Nt;
end
            end
        end
    end
end
nat = ptpn5;
end

function pm = PoiskMax(ar)
for k=1:length(ar)
ar(k)=abs(ar(k));
end
m = ar(1);
for k=1:length(ar)
    if (ar(k)>m)
        m=ar(k);
    end
end
pm = m;
end

function pm = PoiskMin(ar)
m = ar(1);
for k=1:length(ar)
    if (ar(k)<m)
        m=ar(k);
    end
end
pm = m;
end

function [ ar ] = poluchMassivPTP(no,nac,kon,x)
qr=0;
switch (no)
    case 1
        qr = tmp3330();
        case 2
        qr = tmp3360();
        case 3
        qr = tmp33120();
        case 4
        qr = tmp33240();
end
nk=PoiskNacKon(x,nac,kon);
qrr=VydelPodmas(qr,nk(1),nk(2));
ar=qrr;
end

function [ ar ] = poluchMassivTemp(nac,kon,x,tem)
nk=PoiskNacKon(x,nac,kon);
te=VydelPodmas(tem,nk(1),nk(2));
ar=te;
end

function [ ar ] = poluchMassivKoor(nac,kon,x)
nk=PoiskNacKon(x,nac,kon);
xi=VydelPodmas(x,nk(1),nk(2));
ar=xi;
end

function [ pm ] = VydelPodmas(x,nn,nk)
q=1;
for k=nn:nk
    xi(q)=x(k);
    q=q+1;
end
pm=xi;
end

function [ nako ] = PoiskNacKon(x,nac,kon)
fn=1; fk=1; nn=1; nk=1;
for k=1:length(x)
    if (x(k)>nac)
        if (fn>0)
            nn=k;
            fn=0;
        end
    end
    if (x(k)>kon)
        if (fk>0)
            nk=k-1;
            fk=0;
        end
    end
end
nako=[nn nk];
end

function [ kt ] = opredKoefTeplIzPlotTeplPot(te,qv,x)
na=NachToch();
ko=KonechToch();
tm=PoiskMin(te);
q=1; 
for k=1:length(te)
    if (te(k)==tm)
    else
    ktt(q)=qv(k)*abs(x(k)-na)/abs(tm-te(k));
    q=q+1;
    end
end
kt=ktt;
end

function [ kt ] = opredKoefTeplIzPlotTeplPot1000(te,qv,x)
q=1; xn=x(1); 
tn=PoiskMin(te);
for k=2:length(te)
    if (te(k)==tn)
    else
        dx=abs(x(k)-xn);
    ktt(q)=qv(k)*dx/abs(tn-te(k));
    q=q+1;
    end
    tn=te(k);
    xn=x(k);
end
kt=ktt;
end

function [ st ] = opredMasSredTemp(te)
Tn=PoiskMax(te);
q=1;
for k=1:length(te)
    if (te(k)==Tn)
    else
    stt(q)=(te(k)+Tn)/2;
    q=q+1;
    end
end
st=stt;
end

function [ st ] = opredMasSredTemp1000(te)
Tn=PoiskMax(te);
q=1;
for k=1:length(te)
    if (te(k)==Tn)
    else
    stt(q)=te(k);
    q=q+1;
    end
end
st=stt;
end

function [ kt ] = opredKoefTeplIzPlotTeplPotN(no,te,qv)
tol = 30e-3;
switch (no)
    case 1
Tk = 392.9;
T0 = 858.15;
grT=abs(T0-Tk)/tol;
    case 2
Tk = 473.15;
T0 = 1073.15;
grT=abs(T0-Tk)/tol;
    case 3
Tk = 546.15;
T0 = 1273.15;
grT=abs(T0-Tk)/tol;
end
q=1;
for k=1:length(te)
    ktt(k)=qv(k)/grT;
end
%ktt=ktt'
kt=ktt;
end

function [ tem ] = opredPolTem(no,x)
tol = 30e-3;
switch (no)
    case 1
Tk = 392.9;
T0 = 858.15;
Tsre=650.15;
    case 2
Tk = 473.15;
T0 = 1073.15;
Tsre=821.15;
    case 3
Tk = 546.15;
T0 = 1273.15;
Tsre=976.15;
end
tol=30e-3;
koe=polyfit([0 tol/2 tol],[T0 Tsre Tk],2);
for k=1:length(x)
    te(k)=polyval(koe,x(k));
end
tem=te;
end

function [ yre ] = skleimas(x1,x2,x3)
max1=PoiskMax(x1);
max2=PoiskMax(x2);
xr=podfuncskleivmas(x1,x2,max1);
xr=podfuncskleivmas(xr,x3,max2);
xr=-sort(-xr);
yre=xr;
end

function [ ar ] = podfuncskleivmas(xn,xd,ma)
xp=xn; f=1; le=length(xd); q=le;
for k=1:le
    if (f>0)
        if (xd(k)>ma)
            f=0;
            q=k;
        end
    end
end
p=length(xp)+1;
        for k=q:le
            xp(p)=xd(k);
            p=p+1;
        end
ar=xp;
end

function nato = NachToch()
nato=10e-3;
end

function koto = KonechToch()
koto=27e-3;
end

function nato = NachToch1000()
nato=10e-3;
end

function koto = KonechToch1000()
koto=27e-3;
end

function s = opredShagKoo()
s=1e-3*1e-2;
end

function op = opredPorPribl()
op=2;
end

function t = SkhodPosled800()
ct=[30,60,120,240,480];
h=30e-3; co=1; hh=h/ct(co);
x30=0:hh:h; co=co+1; hh=h/ct(co);
y30=0; y30=tmp3330_800();
q=1;
for k=1:length(y30)
    if (y30(k)>0)
       xa30(q)=x30(k);
       ya30(q)=y30(k);
       q=q+1;
    end
end
q30=q-1;
k30=koefPrib(ya30,xa30)
x60=0:hh:h; co=co+1; hh=h/ct(co);
y60=0; y60=tmp3360_800();
q=1;
for k=1:length(y60)
    if (y60(k)>0)
       xa60(q)=x60(k);
       ya60(q)=y60(k);
       q=q+1;
    end
end
q60=q-1;
k60=koefPrib(ya60,xa60)
x120=0:hh:h; co=co+1; hh=h/ct(co);
y120=0; y120=tmp33120_800();
q=1;
for k=1:length(y120)
    if (y120(k)>0)
       xa120(q)=x120(k);
       ya120(q)=y120(k);
       q=q+1;
    end
end
q120=q-1;
k120=koefPrib(ya120,xa120)
x240=0:hh:h; co=co+1; hh=h/ct(co);
y240=0; y240=tmp33240_800();
q=1;
for k=1:length(y240)
    if (y240(k)>0)
       xa240(q)=x240(k);
       ya240(q)=y240(k);
       q=q+1;
    end
end
q240=q-1;
k240=koefPrib(ya240,xa240)
x480=0:hh:h; co=co+1;
y480=0; y480=tmp33480_800();
q=1;
for k=1:length(y480)
    if (y480(k)>0)
       xa480(q)=x480(k);
       ya480(q)=y480(k);
       q=q+1;
    end
end
q480=q-1;
k480=koefPrib(ya480,xa480)
n=1e4; hh=h/n;
x=0:hh:h;
for k=1:length(x)
    s0=0; s1=0; s2=0; s3=0; s4=0;
    for l=1:3
    s0=s0+k30(l)*(x(k)^(l-1));
    s1=s1+k60(l)*(x(k)^(l-1));
    s2=s2+k120(l)*(x(k)^(l-1));
    s3=s3+k240(l)*(x(k)^(l-1));
    s4=s4+k480(l)*(x(k)^(l-1));
    end
    q30(k)=s0; q60(k)=s1; q120(k)=s2; q240(k)=s3; q480(k)=s4;
end
tem=opredPolTem(2,x); 
for k=2:length(x)
    gt=(tem(k)-tem(k-1))/(x(k)-x(k-1));
    lam30(k)=-q30(k)/gt;
    lam60(k)=-q60(k)/gt;
    lam120(k)=-q120(k)/gt;
    lam240(k)=-q240(k)/gt;
    lam480(k)=-q480(k)/gt;
end
%t=postrGraf5(x30,y30,x60,y60,x120,y120,x240,y240,x480,y480);
t=postrGraf5(x,lam30,x,lam60,x,lam120,x,lam240,x,lam480);
end

function t = SkhodPosled1000()
ct=[30,60,120,240];
h=30e-3; co=1; hh=h/ct(co);
x30=0:hh:h; co=co+1; 
hh=h/ct(co);
y30=0; y30=tmp3330_1000();
q=1;
for k=1:length(y30)
    if (y30(k)>0)
       xa30(q)=x30(k);
       ya30(q)=y30(k);
       q=q+1;
    end
end
q30=q-1;
k30=koefPrib(ya30,xa30)
x60=0:hh:h; co=co+1; 
hh=h/ct(co);
y60=0; y60=tmp3360_1000();
q=1;
for k=1:length(y60)
    if (y60(k)>0)
       xa60(q)=x60(k);
       ya60(q)=y60(k);
       q=q+1;
    end
end
q60=q-1;
k60=koefPrib(ya60,xa60)
x120=0:hh:h; co=co+1; 
hh=h/ct(co);
y120=0; y120=tmp33120_1000();
q=1;
for k=1:length(y120)
    if (y120(k)>0)
       xa120(q)=x120(k);
       ya120(q)=y120(k);
       q=q+1;
    end
end
q120=q-1;
k120=koefPrib(ya120,xa120)
x240=0:hh:h;
y240=0; y240=tmp33240_1000();
for k=1:length(y240)
    if (y240(k)>0)
       xa240(q)=x240(k);
       ya240(q)=y240(k);
       q=q+1;
    end
end
q240=q-1;
k240=koefPrib(ya240,xa240)
n=1e4; hh=h/n;
x=0:hh:h;
for k=1:length(x)
    s0=0; s1=0; s2=0; s3=0;
    for l=1:3
    s0=s0+k30(l)*(x(k)^(l-1));
    s1=s1+k60(l)*(x(k)^(l-1));
    s2=s2+k120(l)*(x(k)^(l-1));
    s3=s3+k240(l)*(x(k)^(l-1));
    end
    q30(k)=s0; q60(k)=s1; q120(k)=s2; q240(k)=s3;
end
tem=opredPolTem(2,x); 
for k=2:length(x)
    gt=(tem(k)-tem(k-1))/(x(k)-x(k-1));
    lam30(k)=-q30(k)/gt;
    lam60(k)=-q60(k)/gt;
    lam120(k)=-q120(k)/gt;
    lam240(k)=-q240(k)/gt;
end
t=postrGraf4(x,lam30,lam60,lam120,lam240);
end

function [ koef ] = koefPrib(ktp,te)
	yx2=0; yx=0; x4=0; x3=0; x2=0; x=0; y=0; de=0; de1=0; de2=0; de3=0; le=length(te);
for k=1:le
        yx2=yx2+ktp(k)*(te(k)^2); 
        yx=yx+ktp(k)*te(k); 
        y=y+ktp(k); 
        x4=x4+(te(k)^4); 
        x3=x3+(te(k)^3); 
        x2=x2+(te(k)^2); 
        x=x+te(k);
end
	b(1)=yx2; b(2)=yx; b(3)=y; 
    A(1,1)=x4; A(1,2)=x3; A(1,3)=x2; 
    A(2,1)=A(1,2); A(2,2)=A(1,3); A(2,3)=x; 
    A(3,1)=A(1,3); A(3,2)=A(2,3); A(3,3)=le;
	de=A(1,1)*(A(3,3)*A(2,2)-A(3,2)*A(2,3))-A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2));
	de1=b(1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))-A(1,2)*(b(2)*A(3,3)-b(3)*A(2,3))+A(1,3)*(b(2)*A(3,2)-b(3)*A(2,2));
	de2=A(1,1)*(b(2)*A(3,3)-b(3)*A(2,3))-b(1)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))+A(1,3)*(A(2,1)*b(3)-A(3,1)*b(2));
	de3=A(1,1)*(A(2,2)*b(3)-A(3,2)*b(2))-A(1,2)*(A(2,1)*b(3)-A(3,1)*b(2))+b(1)*(A(2,1)*A(3,2)-A(3,1)*A(2,2));
	ko(3)=de1/de; ko(2)=de2/de; ko(1)=de3/de; ko=ko';
    %ko=inv(A)*b'
koef=ko;
end

function t = postrGraf5(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5)
pl=plot(x1,y1,'-b',x2,y2,'-k',x3,y3,'-g',x4,y4,'-r',x5,y5,'-m');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Координата, м'}); 
ylabel({'Отношение'}); 
title({'График зависимости отношения от x'});
legend('30/30','60/30','120/30','240/30','480/30','location','best');
t=0;
end

function t = postrGraf4(x,y1,y2,y3,y4)
pl=plot(x,y1,'-b',x,y2,'-k',x,y3,'-g',x,y4,'-r');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Координата, м'}); 
ylabel({'Отношение'}); 
title({'График зависимости отношения от x'});
legend('30/30','60/30','120/30','240/30','480/30','location','best');
t=0;
end

function t = postrGraf3(x,y1,y2,y3)
pl=plot(x,y1,'-b',x,y2,'-k',x,y3,'-g');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Координата, м'}); 
ylabel({'Отношение'}); 
title({'График зависимости отношения от x'});
legend('60/30','120/30','240/30','location','best');
t=0;
end

function pg = postrGrafikOtn(otn60,otn120,otn240,otnb,v,Nt)
y=0; y=1:1:length(otn60); z2=trapz(y,otn60)/length(otn60); 
y=0; y=1:1:length(otn120); z3=trapz(y,otn120)/length(otn120); 
y=0; y=1:1:length(otn240); z4=trapz(y,otn240)/length(otn240);
y=0; y=1:1:length(otnb); z5=trapz(y,otnb)/length(otnb); 
v=v'
Nt=Nt
z=[z2,z3,z4,z5]
if (v==2e3)
x=[6e1,12e1,2e2,1e3];
%kxz=polyfit(x,z,2);
%p=30; x=0; z=0; nt=2e1;
%for k=1:nt
    %x(k)=k*p;
    %z(k)=polyval(kxz,x(k));
%end
pl=plot(x,z,'-b');
set(pl,'LineWidth',2); 
hold on; grid on;
xlabel({'Число точек'}); 
ylabel({'Среднее отношение'}); 
title({'График зависимости СО от ЧТ'});
end
pg=z5;
end

function pg = postrGrafikOtn800(otn60,otn120,otn240,otn480,otnb,v)
y=0; y=1:1:length(otn60); z2=trapz(y,otn60)/length(otn60); 
y=0; y=1:1:length(otn120); z3=trapz(y,otn120)/length(otn120); 
y=0; y=1:1:length(otn240); z4=trapz(y,otn240)/length(otn240);
y=0; y=1:1:length(otn480); z5=trapz(y,otn480)/length(otn480); 
y=0; y=1:1:length(otnb); z6=trapz(y,otnb)/length(otnb); 
v=v'
z=[z2,z3,z4,z5,z6]
%x=[60,120,240,480,960];
%pl=plot(x,z,'-b');
%set(pl,'LineWidth',2); 
%hold on; grid on;
%xlabel({'Число точек'}); 
%ylabel({'Среднее отношение'}); 
%title({'График зависимости СО от ЧТ'});
pg=0;
end
%номер 1
function [qr] = tmp3330()
qr = [1162.19783592224
           14064.068862915
          6346.05541419983
          4286.39482116699
          3472.70532226563
           3028.0277338028
          2659.69601249695
          2405.89868545532
          2192.42489528656
          2001.37012958527
          1829.24831676483
          1708.98024463654
          1586.55086994171
          1321.23091983795
          1162.19783687592
          1023.62346744537
          896.049344062805
          772.326431274414
          585.763700008392
          589.431444168091
          533.973300933838
          457.497436523438
           384.29551410675
          316.508916854858
          244.026017665863
          207.390486717224
          152.802119731903
          72.2684442996979
         -119.025223016739
         -1004.04144668579
           1162.1978366375];
end
%номер 1
function [qr] = tmp3360()
qr = [1552.0892868042
          9185.35584259033
          6035.39772033691
          5175.50819015503
          4785.06861495972
          4537.32759857178
          4348.20451545715
          4191.49641990662
          4067.01355934143
          4056.79320526123
          3817.47589206696
          3636.36074352264
          3490.86863517761
          3355.54878330231
          3225.53586196899
          3099.31887149811
          2976.33177757263
          2856.39799213409
          2739.59995079041
          2626.34840106964
          2517.77361869812
          2417.29804325104
          2338.86431980133
          2385.99097919464
           2220.2259683609
           1982.0785074234
          1852.70462989807
          1744.91805076599
          1645.86817836761
           1552.0892868042
          1462.24687290192
          1375.62002277374
          1291.57381534576
          1209.07296657562
          1125.28607082367
          1027.98107242584
          813.248772621155
           846.45681142807
          887.868816375732
          847.372228622437
          794.554096221924
          739.572588920593
          685.139139175415
          632.228363990784
          581.222216129303
          532.189800262451
          484.784804344177
          436.917844772339
          367.961401939392
          357.759227752686
          332.826685905457
          298.384731054306
           262.56122303009
            225.1837849617
          183.636612653732
          131.819003820419
          55.6087141036987
          -82.102441072464
         -408.404317378998
         -1685.50518465042
          1552.08928656578];
end
%номер 1
function [qr] = tmp33120()
qr = [1814.74082565308
         -7636.51644706726
         -2253.04853248596
          456.741205215454
          1988.66799545288
          2917.73682403564
          3507.84474563599
          3893.47611427307
          4148.69856071472
          4316.98822784424
          4425.27988433838
          4491.07739830017
          4526.28287506104
          4539.42804145813
          4537.15787887573
          4525.59770011902
          4512.67169761658
          4515.96632575989
          4626.43507766724
          4492.14957618713
          4344.27829837799
          4254.58579063416
          4178.18642520905
          4106.29050540924
          4036.11470127106
          3966.64147758484
          3897.49417209625
          3828.55362796783
          3759.80892562866
          3691.29484272003
          3623.06427669525
          3555.17665195465
          3487.69236183167
          3420.67242717743
          3354.17903614044
          3288.27887916565
            3223.047955513
          3158.57985877991
          3094.99973773956
          3032.48939704895
          2971.33554363251
          2912.02960777283
          2855.49212932587
          2803.63252162933
          2760.95423316956
          2740.42283916473
           2803.0500793457
          2962.44169235229
            2584.498878479
          2444.38677978516
          2351.87015628815
          2276.74510860443
          2209.79420185089
          2147.40243339539
          2087.92483901978
          2030.52563285828
          1974.74057769775
          1920.29135513306
          1866.99914264679
          1814.74082660675
          1763.42448139191
          1712.97414207459
          1663.31838130951
          1614.37917900085
          1566.05781555176
          1518.21189403534
          1470.61194324493
          1422.84848880768
          1374.11159133911
          1322.59992599487
          1263.66916370392
          1182.19770622253
          986.968205451965
          860.729664325714
          1084.33601951599
          1110.97339344025
          1099.91797304153
          1076.01752567291
          1046.81836462021
          1015.25959396362
          982.672360420227
           949.73131942749
           916.80982875824
          884.128980636597
          851.826354503632
          819.990152835846
          788.677283763886
          757.922862529755
          727.744901657104
          698.144588947296
          669.101470470428
           640.55796957016
          612.376822948456
          584.216162204742
          555.102318286896
          521.452716827393
           455.58375120163
          451.684728145599
          464.113433361053
          449.725965023041
          430.617165088654
          410.186892032623
          389.439684152603
          368.727280139923
           348.15896987915
          327.724592447281
          307.330778837204
          286.802506923676
          265.863391876221
          244.094081163406
          220.859271287918
          195.179604291916
          165.499279737473
           129.23649263382
          81.8357443809509
          14.5360169410706
         -91.6901812553406
         -286.029559612274
          -740.31570982933
         -2052.89181065559
          1814.74082636833];
end
%номер 1
function [qr] = tmp33240()
qr = [1970.20248031616
         -31220.7749214172
         -19716.5008220673
          -12596.164560318
         -8046.89838790894
         -4958.72407150269
         -2767.79908180237
         -1161.32359695435
          47.5566635131836
          976.589519500732
          1703.03269195557
          2279.25423812866
          2741.73358535767
          3116.49446105957
          3422.49774360657
          3673.82295036316
          3881.10730934143
          4052.51424026489
          4194.39981842041
          4311.78012084961
          4408.66309738159
          4488.28875541687
          4553.30612945557
          4605.90591812134
          4647.92175865173
          4680.90980529785
          4706.21389007568
          4725.02264213562
          4738.42486381531
          4747.47177696228
          4753.26040458679
           4757.0667591095
          4760.59741973877
          4766.55429649353
           4780.2092628479
          4815.52501869202
            4940.104139328
          4949.33024215698
          4766.13354110718
          4686.36420822144
          4631.60130119324
          4586.81535243988
          4546.84271907806
          4509.39011383057
          4473.29465675354
          4437.92132091522
          4402.90920162201
          4368.04951953888
          4333.22207260132
          4298.35969638824
           4263.4276676178
          4228.41139030457
          4193.30872058868
            4158.125207901
          4122.87106704712
          4087.55907535553
          4052.20349121094
          4016.81917190552
          3981.42101669312
          3946.02372360229
          3910.64149475098
          3875.28815841675
          3839.97689628601
          3804.72034358978
          3769.53068637848
          3734.41972637177
          3699.39887428284
          3664.47936820984
          3629.67243766785
          3594.98933982849
          3560.44167613983
          3526.04159355164
          3491.80209445953
          3457.73739433289
          3423.86344909668
          3390.19858837128
          3356.76440143585
          3323.58691310883
          3290.69825935364
          3258.13906097412
          3225.96183681488
          3194.23612689972
          3163.05612659454
          3132.55262947083
           3102.9120054245
          3074.40742969513
          3047.45194244385
          3022.69323348999
          3001.19326496124
          2984.79844284058
           2977.0055141449
          2985.42952346802
          3031.61229324341
          3225.95837306976
          3335.63265037537
          2952.88038825989
          2813.04871463776
          2727.66100311279
          2664.37773227692
          2612.51476097107
          2567.33965396881
          2526.40474224091
          2488.32068347931
          2452.24782848358
          2417.65445804596
          2384.19088363647
           2351.6197013855
          2319.77520275116
          2288.53873634338
          2257.82342243195
          2227.56434440613
          2197.71199989319
          2168.22799015045
          2139.08204555511
          2110.24983596802
          2081.71152877808
          2053.45066928864
            2025.453332901
          1997.70747375488
          1970.20248031616
          1942.92867946625
          1915.87703609467
          1889.03879642487
          1862.40519428253
          1835.96702480316
          1809.71438884735
          1783.63613319397
          1757.71933364868
          1731.94851875305
          1706.30462741852
           1680.7635974884
          1655.29427337646
          1629.85530090332
          1604.39050388336
          1578.82164669037
          1553.03687095642
          1526.87164497375
          1500.07653331757
           1472.2600440979
          1442.78264045715
           1410.5472240448
          1373.54495429993
          1327.72224903107
          1263.42399978638
          1148.89408874512
          768.239392757416
          943.190466403961
          1164.33655261993
          1223.22218418121
          1243.29546785355
          1247.07407188416
          1242.38956785202
           1232.8266453743
          1220.27385044098
          1205.81351709366
          1190.10302305222
          1173.55989933014
          1156.45930051804
          1138.98845148087
          1121.27838420868
          1103.42310810089
          1085.49159526825
          1067.53546047211
          1049.59401702881
          1031.69765043259
          1013.87015628815
          996.130392551422
          978.493411064148
          960.971312999725
          943.573853969574
          926.308910369873
          909.182807445526
          892.200545787811
          875.365994930267
          858.681985378265
          842.150378704071
          825.772023200989
          809.546695232391
          793.472922801971
          777.547636508942
          761.765697479248
          746.118982791901
          730.595011234283
          715.174569606781
          699.827630996704
          684.506132125854
          669.130552768707
          653.563487052917
           637.55349111557
          620.599675655365
          601.548810005188
          576.853490352631
          524.652067661285
          454.400665283203
          534.574836730957
          542.930734634399
          539.416782855988
          531.799297332764
           522.30883026123
          511.861357688904
          500.906450748444
          489.687959194183
          478.346875190735
          466.968427181244
          455.605731487274
          444.292451858521
          433.049898862839
          421.891123533249
          410.823327541351
          399.849308490753
          388.968273162842
          378.176292181015
          367.466437101364
          356.828742265701
           346.24994134903
          335.713087320328
          325.196878671646
          314.674852609634
          304.114233970642
          293.474433660507
          282.705011367798
          271.742968559265
            260.5090675354
          248.902766704559
          236.795083761215
          224.018395185471
          210.351449489594
          195.496812105179
          179.045773983002
          160.421983718872
          138.787223815918
          112.876635313034
          80.6964244842529
           38.931613445282
         -18.3179869651794
         -102.681166410446
         -243.271758317947
         -555.190999746323
         -832.873386383057
         -1765.83007764816
          1970.20247983933];
end
%
function [qr] = tmp3360_800()
qr = [4940.62083053589
          116.249317169189
          9026.36068344116
          11111.8548088074
          11731.6630821228
          11813.7915668488
          11615.5649452209
          11216.7193107605
          10573.8286209106
          9270.23252868652
          4648.02250289917
          8488.08448410034
          9029.62323760986
          9031.69982337952
          8877.30480384827
          8661.40188217163
          8417.02228927612
          8158.99861526489
          7895.95744895935
          7635.72604942322
          7394.14087677002
          7255.56439590454
          7144.02194976807
          6569.84770011902
          6239.71360206604
          5953.37840652466
          5683.98192214966
          5425.49061584473
          5176.96906280518
          4940.62083053589
          4725.38876724243
          4570.84544944763
          4793.49223136902
          4079.42044734955
          3714.16229820251
          3462.15534305573
          3242.30223369598
          3039.55582046509
           2850.6335401535
          2678.58254337311
          2544.24087905884
          2624.98449611664
          2313.28330421448
           1920.6060667038
          1726.48462581635
          1570.68892478943
          1430.55772399902
          1298.90544033051
           1169.2837972641
          1018.66203212738
          668.085024356842
          905.660008907318
          857.969046115875
          771.452642917633
          671.360111236572
          555.416356086731
          404.382448196411
          159.550666809082
         -387.515203475952
         -2146.03608226776
          4940.62083005905];
end
%
function [qr] = tmp3330_800()
qr = [3712.99695968628
          25242.5880432129
          14908.4313716888
          11200.8934841156
          8962.72780990601
          5140.19753265381
          6836.37640953064
          6610.66784667969
          6212.23237037659
          5793.05110740662
          5412.17974662781
           5131.3022518158
          4518.83023643494
          4094.99907493591
          3712.99695777893
          3375.99701118469
          3274.16922187805
          2630.11152553558
           2270.5757484436
          1982.89824104309
          1788.73041248322
          1570.86688327789
          1164.01247596741
          947.008038520813
          759.911653518677
          501.979877471924
          531.469252586365
          369.853533744812
          37.5314788818359
         -1349.32714080811
          3712.99695849419];
end
%
function [qr] = tmp3360_1000()
qr = [6394.99278259277
         -33046.8503341675
          5989.45297622681
          15751.3901405334
          19502.8326301575
          21424.4204864502
          22677.0461959839
          23933.4703330994
          26974.2777175903
          31352.8488998413
          25314.7089996338
          23685.6282577515
          22768.7478561401
          22054.3345718384
          21445.8733406067
          20956.3416862488
          20687.9331359863
          20974.9811096191
          23203.5645637512
          37327.3307723999
          20201.9854545593
          15568.5815544128
          13557.8446846008
          12206.0589637756
          11087.7567672729
          10015.4947357178
          8753.07412719727
          6449.02854537964
         -4485.12927246094
          6394.99278068542
           8188.0920791626
          8464.16403579712
           8331.4283618927
          8052.94213676453
          7726.77498817444
          7433.71738052368
          7551.66881942749
          6582.39527320862
          6058.27027511597
          5650.21768951416
             5283.85455513
          4952.22469520569
          4684.63697814941
          4756.69906234741
          4233.12545204163
          3552.00262355804
          3203.72363471985
          2923.92358112335
          2681.48833656311
          2491.63113307953
          2580.04298400879
           2052.1155500412
          1658.19568729401
          1434.84227657318
          1238.36446857452
          1029.41250228882
          700.051002502441
          305.497896671295
         -92.8370852470398
         -2361.44493865967
          6394.99278068542];
end
%
function [qr] = tmp3330_1000()
qr = [-167.5641746521
          34853.6558799744
          26096.9808158875
          22215.9173736572
          23262.8381614685
          21320.8252067566
           17935.626750946
          16524.8756942749
          15958.9887466431
          19898.0620384216
          17154.7048492432
          10179.4129066467
          7996.81916236877
          6103.42633056641
         -167.564176559448
          5559.40611457825
          5772.69713020325
          5383.35733032227
          5105.89642715454
          4206.68611526489
          3645.88787460327
           3263.3914642334
          2860.85869979858
          2156.79219150543
          1780.01635932922
          1602.83222198486
           1067.8166885376
          717.150436401367
          171.456097602844
         -1520.63088703156
         -167.564176559448];
end
%
function [qr] = tmp33120_1000()
qr = [9420.73093032837
         -155183.781837463
         -66936.7010231018
         -28365.4500083923
         -8831.11661911011
          2247.69164276123
           9107.5241317749
          13659.7352638245
          16850.1640090942
          19185.7496452332
           20960.018081665
          22358.3786621094
          23513.0329551697
          24538.9919052124
          25572.7650222778
          26850.1250267029
          28964.4346923828
          34487.2041511536
           35830.237537384
          30027.5343742371
          28107.8708000183
          27145.8543014526
           26534.771686554
          26074.5751037598
          25686.1667633057
          25335.5980949402
          25008.3617172241
          24699.4909362793
          24409.6317520142
          24143.8322639465
          23912.0516967773
          23731.5120582581
          23632.1996536255
          23669.1872367859
          23951.8962516785
          24721.1854858398
          26587.0143814087
          31486.2938613892
          49678.1666564941
          30591.1053199768
          22879.4136886597
          19696.3968544006
          17907.5336036682
          16698.3344535828
          15773.6047210693
          15005.5684051514
           14330.621761322
          13712.3770866394
          13126.1796951294
           12551.068359375
          11963.4977550507
          11328.6775283813
          10582.1460762024
          9579.41772842407
          7926.60000038147
           4217.9802570343
         -10579.8073883057
          3909.01582717896
          7953.18844604492
          9420.73092842102
           10015.272359848
          10234.8827724457
          10268.6153793335
          10200.5404930115
          10072.8398113251
          9908.96371459961
          9723.33054161072
          9526.04271697998
          9325.95534706116
          9134.36800384521
          8974.02008628845
           8915.6855430603
          9379.95199966431
          8472.79160499573
           7949.3731918335
          7619.75630950928
          7347.82665634155
          7101.57856559753
          6870.05511283875
          6648.92259216309
          6436.54277801514
           6232.8129940033
          6039.14858627319
          5859.61838722229
          5704.95534896851
          5607.98535346985
          5721.09216499329
          6175.49771308899
          5032.60343837738
          4622.24172973633
          4358.81524372101
          4147.92174816132
          3962.19000434875
            3791.486577034
          3631.59986209869
          3480.98474121094
          3340.00006866455
          3211.74836921692
          3106.34817123413
          3060.46950912476
          3301.99976825714
           3151.6142282486
          2515.00644779205
           2270.6501750946
          2107.11848831177
          1972.93270969391
          1852.66237640381
          1739.97707462311
          1631.29852962494
          1523.53680896759
          1412.33931732178
          1288.19648075104
          1119.70757007599
          696.982045650482
          725.543639183044
          864.819704532623
          717.273613929749
          392.544324398041
         -309.508699417114
         -2496.40636491776
          9420.73092937469];
end
%
function [qr] = tmp33200_1000()
qr = [11216.0245895386
         -146484.314178467
         -64455.2052383423
         -31363.6900253296
         -15381.5334472656
         -6564.07315063477
          -1176.6195602417
          2383.59459686279
          4879.35047149658
           6706.9953918457
          8089.07841491699
          9159.38584136963
          10003.2061843872
          10677.6609039307
          11222.5688476563
          11666.6048049927
          12030.9618530273
          12331.6346206665
          12580.9004745483
          12788.3155288696
          12961.4045181274
          13106.1500244141
           13227.348815918
          13328.8758468628
          13413.8794021606
          13484.9323348999
          13544.1469268799
          13593.2664031982
          13633.7339019775
          13666.7497940063
          13693.3138580322
          13714.2626113892
          13730.2975616455
          13742.0084991455
          13749.8924560547
          13754.3691864014
          13755.7934417725
           13754.466217041
          13750.6432266235
          13744.5416183472
          13736.3475112915
          13726.2190246582
           13714.291267395
           13700.679901123
          13685.4842605591
          13668.7895278931
          13650.6689147949
          13631.1854400635
          13610.3936157227
          13588.3403625488
          13565.0662460327
          13540.6067352295
          13514.9924316406
          13488.2499847412
          13460.4027328491
          13431.4710540771
          13401.4728927612
          13370.4240264893
          13338.3384552002
          13305.2286148071
          13271.1055870056
          13235.9793586731
          13199.8589248657
          13162.7524719238
          13124.6674880981
          13085.6108665466
          13045.5890312195
          13004.6079406738
          12962.6732597351
          12919.7903556824
          12875.9643440247
          12831.2001609802
           12785.502620697
          12738.8763923645
           12691.326084137
          12642.8562316895
          12593.4713478088
          12543.1759185791
          12491.9744567871
          12439.8714599609
          12386.8714981079
          12332.9791603088
          12278.1990699768
          12222.5359611511
          12165.9946250916
          12108.5798988342
          12050.2967567444
          11991.1502189636
          11931.1454315186
           11870.287612915
          11808.5821304321
           11746.034412384
            11682.65001297
           11618.434627533
          11553.3940391541
          11487.5341720581
          11420.8610687256
          11353.3808822632
          11285.0999221802
          11216.0245895386
          11146.1614723206
          11075.5172309875
          11004.0987167358
          10931.9129142761
          10858.9669113159
            10785.26795578
          10710.8234291077
          10635.6409263611
          10559.7280960083
           10483.092792511
           10405.743019104
          10327.6869239807
           10248.932800293
          10169.4891052246
          10089.3644943237
          10008.5677375793
          9927.10776138306
          9844.99373245239
          9762.23492431641
          9678.84077453613
          9594.82095909119
           9510.1852684021
          9424.94370269775
          9339.10641098022
          9252.68379592896
          9165.68639755249
          9078.12491607666
          8990.01033973694
          8901.35376358032
          8812.16653633118
          8722.46019744873
          8632.24646377563
          8541.53729057312
          8450.34485435486
          8358.68148040771
          8266.55982971191
          8173.99266052246
          8080.99304580688
          7987.57424354553
          7893.74973487854
          7799.53328514099
          7704.93884277344
          7609.98059654236
          7514.67300033569
          7419.03072166443
          7323.06866836548
          7226.80202484131
          7130.24614906311
           7033.4166431427
           6936.3293838501
          6839.00044059753
          6741.44609260559
          6643.68284225464
           6545.7274017334
          6447.59657859802
          6349.30738449097
           6250.8769493103
          6152.32249832153
          6053.66131210327
          5954.91069602966
            5856.087890625
             5757.21002388
          5658.29402923584
          5559.35658836365
          5460.41390228271
          5361.48163223267
          5262.57469844818
          5163.70705223084
          5064.89145946503
           4966.1391658783
          4867.45950984955
          4768.85945796967
          4670.34303569794
          4571.91058635712
          4473.55801200867
           4375.2755613327
          4277.04656791687
          4178.84566497803
          4080.63671588898
          3982.37002277374
          3883.97884941101
          3785.37466716766
          3686.44110488892
          3587.02578830719
          3486.92912960052
          3385.88877391815
          3283.55739593506
          3179.47006988525
           3072.9947052002
          2963.25356388092
          2848.99364280701
          2728.35952186584
          2598.47651195526
          2454.63258743286
          2288.56669902802
           2084.6130900383
          1810.22138929367
          1389.95161962509
          622.526295185089
         -1185.07016706467
          11216.0245909691];
end
%
function [qr] = tmp33120_800()
qr=[5906.48500061035
         -56234.5102043152
         -22517.0593032837
         -7502.68809127808
          283.674222946167
          4762.09992218018
          7527.73380851746
          9317.89939498901
          10507.8760509491
          11304.3392295837
          11828.7574367523
          12156.1579837799
          12334.1821994781
           12392.576417923
          12347.2657756805
          12200.0282917023
          11931.6858730316
          11479.4782924652
          10660.1860561371
          8801.14775848389
           3487.4055557251
          8647.83493232727
          10012.5517253876
          10524.2363948822
          10721.4740200043
          10774.9334354401
          10753.0019493103
          10687.3333702087
           10594.247127533
          10482.9356765747
           10358.982339859
          10226.0069656372
          10086.4958076477
          9942.25247383118
          9794.66342353821
          9644.87369728088
          9493.93020439148
          9342.93741416931
          9193.30318641663
          9047.23830223084
          8908.98868751526
          8788.50604057312
          8715.94795799255
          8851.16163063049
          8677.65129852295
          8202.16517829895
          7948.22588348389
          7747.27297973633
          7567.13365936279
           7397.3409576416
          7233.66926574707
          7074.18988418579
          6917.95935249329
          6764.51524353027
          6613.66757392883
           6465.4163608551
          6319.94137763977
          6177.65207672119
          6039.33289909363
           5906.4850025177
          5782.16233062744
          5673.23263931274
          5597.57726287842
          5614.77782821655
          6090.17800140381
           5429.9732131958
          4991.43873119354
          4756.27224731445
          4577.90929412842
          4423.75815200806
          4282.41945743561
           4149.0540971756
          4021.34852600098
          3898.12018299103
           3778.7738904953
          3663.08785152435
          3551.16466140747
          3443.52314472198
          3341.42401885986
          3247.78385257721
          3169.95864391327
           3130.3220949173
          3238.52266597748
          3506.18637752533
          2810.09126472473
           2568.7318944931
          2417.01377010345
           2297.3727645874
          2192.76275062561
          2096.65807914734
          2006.09561157227
          1919.52185821533
          1835.99187469482
          1754.80836868286
          1675.29482841492
          1596.53894233704
          1516.87781810761
          1432.45088386536
          1332.00006532669
          1168.75635480881
          692.405408382416
          1139.32972192764
          1192.52253818512
          1174.62128448486
          1133.93118047714
          1084.23028612137
          1030.49906158447
           974.72641658783
          917.626659393311
          859.199300289154
          798.866919517517
          735.380692005157
          666.471543312073
          588.000980854034
          491.644455432892
          354.306290149689
          127.660443782806
         -66.9841585159302
         -637.346429824829
         -2426.43390893936
          5906.48500227928];
end
%
function [qr] = tmp33240_800()
qr=[6170.37331771851
         -116459.676094055
         -75176.2152671814
         -48982.1040916443
         -32271.2683601379
          -21061.428150177
         -13217.4876670837
          -7539.6887512207
         -3314.40962982178
          -96.793025970459
           2401.2052116394
          4372.05837249756
          5947.91414260864
          7221.70766830444
          8260.21736907959
          9112.38132667542
          9814.71392250061
          10394.9097881317
           10874.289522171
          11269.4843425751
          11593.6145553589
          11857.1265659332
          12068.3910675049
           12234.131570816
           12359.729598999
          12449.4351005554
          12506.5000896454
          12533.2401313782
          12531.0206604004
          12500.1506690979
          12439.6463394165
          12346.7907772064
          12216.3495483398
          12039.1634044647
          11799.5329399109
          11470.0574455261
          11000.4858779907
          10290.0667114258
          9101.22247314453
          6630.08684158325
          1542.17343139648
          7058.34036064148
          8947.79386901855
          9891.13684463501
          10430.8250751495
          10759.2680206299
           10964.187379837
          11091.4834880829
          11167.3659305573
            11207.74858284
          11222.7576217651
          11219.0852622986
          11201.3000068665
          11172.6103382111
          11135.3280334473
          11091.1580638886
          11041.3850841522
          10986.9963741302
          10928.7645606995
          10867.3046016693
          10803.1136951447
          10736.5997257233
          10668.1018695831
          10597.9060497284
          10526.2563648224
          10453.3640804291
          10379.4145946503
           10304.573135376
          10228.9894828796
          10152.8027610779
          10076.1450138092
          9999.14554786682
           9921.9359588623
          9844.65574836731
          9767.46020317078
          9690.53105926514
          9614.09212303162
          9538.43360900879
          9463.95032501221
          9391.20531463623
          9321.03959846497
          9254.77200889587
          9194.58765220642
           9144.3638458252
          9111.66882133484
          9113.65755844116
          9201.34301948547
           9652.9100189209
          9021.86739730835
          8755.89543151855
          8582.95470428467
          8447.71223831177
          8331.74147224426
          8226.90579223633
          8129.03923797607
          8035.80547332764
          7945.81003952026
           7858.1811542511
          7772.35375785828
          7687.95070457459
          7604.71421432495
          7522.46493148804
          7441.07660865784
          7360.45986747742
          7280.55187988281
          7201.30966758728
          7122.70561027527
           7044.7246761322
          6967.36276435852
          6890.62608337402
          6814.53127861023
          6739.10648727417
           6664.3932056427
          6590.44968223572
          6517.35591506958
          6445.22141075134
          6374.19689178467
          6304.49256324768
          6236.40690803528
          6170.37331771851
          6107.03818702698
          6047.39669799805
          5993.04018592834
          5946.63592147827
          5912.94114112854
          5901.23389816284
          5932.42294502258
          6068.20018005371
          6646.11264228821
          6292.46021461487
          5705.52316856384
          5446.64347934723
          5277.43813991547
          5148.05821418762
          5040.12673950195
           4945.0497303009
          4858.26289463043
          4777.13983535767
          4700.08052539825
          4626.06620883942
          4554.42443561554
          4484.69792556763
          4416.56801891327
          4349.80839824677
          4284.25636577606
          4219.79456710815
          4156.33930397034
          4093.83299827576
          4032.23976802826
          3971.54295825958
           3911.7447977066
          3852.86775302887
          3794.95813083649
          3738.09276485443
          3682.39014053345
          3628.02879428864
          3575.27772426605
          3524.54814147949
          3476.48398113251
          3432.12757968903
          3393.24034881592
          3362.97324752808
           3347.4417963028
          3360.19248867035
          3439.88349342346
           3789.4089345932
          3973.54220962524
          3262.25332736969
          3013.43877124786
          2866.82316207886
          2761.53107357025
          2677.32922935486
          2605.30400753021
           2540.9155664444
          2481.63559436798
          2425.96311473846
          2372.96064281464
          2322.01679229736
          2272.71716403961
          2224.77073287964
          2177.96598434448
          2132.14372062683
          2087.17958068848
          2042.97212982178
          1999.43423557281
          1956.48628997803
          1914.04997825623
           1872.0417137146
          1830.36455059052
          1788.89683532715
          1747.47461891174
          1705.86290645599
          1663.70536136627
          1620.43243122101
           1575.0832529068
          1525.93644714355
          1469.66894435883
          1399.13008642197
          1295.63987588882
          1081.89634561539
          467.026814937592
          1120.66627264023
          1265.27351093292
          1313.02754163742
          1324.75324058533
          1319.19061040878
          1304.25861787796
          1283.87868833542
          1260.18419694901
          1234.41350317001
           1207.3216381073
          1179.38711500168
          1150.92216825485
          1122.13438367844
          1093.16236019135
          1064.09695148468
          1034.99388647079
          1005.88131046295
          976.763751029968
          947.623646259308
          918.420662879944
          889.088834762573
          859.530854701996
          829.608528614044
          799.126836299896
          767.806762218475
          735.236644744873
          700.777718544006
          663.357471942902
          620.915256977081
          568.221965789795
          477.541446685791
          352.776769638062
          414.952528953552
          366.178522109985
            274.7876496315
          129.133636951447
          -117.53421831131
         -596.400901317596
         -1775.43962860107
          6170.37331843376];
end
%
function [qr] = tmp33480_800()
qr=[6365.85315322876
         -156126.400913239
         -132917.058761597
         -106347.449363709
         -85177.3639755249
         -68725.0702972412
         -55812.2950248718
         -45525.1691474915
         -37210.1021194458
           -30401.09608078
         -24761.0919342041
         -20041.8831748962
         -16057.4604492188
         -12666.1934547424
         -9758.70100021362
         -7249.44915008545
         -5070.80916976929
         -3168.77101516724
         -1499.81531524658
         -28.5915069580078
          1273.83344650269
          2431.30432891846
          3463.59963989258
          4387.22564315796
          5216.03828048706
           5961.7306022644
          6634.22212600708
          7241.96993637085
          7792.21823883057
          8291.20061683655
          8744.30418395996
           9156.2039642334
          9530.97281074524
          9872.17265510559
          10182.9295711517
          10465.9966602325
          10723.8060092926
          10958.5125656128
          11172.0308704376
          11366.0660915375
          11542.1402359009
           11701.614528656
           11845.708065033
          11975.5138931274
          12092.0126953125
          12196.0838088989
           12288.515417099
          12370.0123462677
          12441.2027759552
          12502.6436862946
          12554.8250026703
          12598.1726341248
          12633.0503158569
          12659.7602348328
          12678.5423793793
          12689.5725479126
          12692.9584980011
          12688.7342395782
          12676.8518600464
          12657.1704044342
          12629.4406204224
          12593.2847499847
          12548.1696434021
          12493.3706245422
          12427.9231128693
          12350.5567684174
          12259.6052265167
          12152.8803443909
           12027.493473053
          11879.5968112946
          11703.9978866577
          11493.5675888062
          11238.2942142487
          10923.6961898804
          10527.9878673553
          10016.5880718231
          9330.23612976074
          8354.95781898499
          6826.48393821716
          3872.22381019592
          198.883491516113
          5105.48974227905
           7197.6478805542
          8403.66369438171
          9196.41744422913
          9755.48052787781
          10166.8167514801
          10477.9383277893
          10717.7193241119
          10904.9023799896
          11052.2511043549
          11168.7788448334
          11261.0342159271
          11333.8856544495
          11391.0218601227
          11435.2826442719
          11468.8849582672
          11493.5805053711
          11510.7685623169
          11521.5779342651
          11526.9272003174
          11527.5702991486
          11524.1307983398
          11517.1282291412
          11506.9983711243
          11494.1092090607
          11478.7734184265
          11461.2582416534
          11441.7933578491
          11420.5773715973
          11397.7828273773
          11373.5603313446
          11348.0421199799
          11321.3447551727
          11293.5714302063
          11264.8139324188
           11235.154209137
          11204.6657714844
          11173.4147453308
          11141.4608726501
          11108.8583240509
          11075.6563644409
          11041.8999595642
          11007.6302967072
           10972.885187149
          10937.6995220184
          10902.1055717468
          10866.1332206726
          10829.8103523254
          10793.1629428864
          10756.2153873444
          10718.9906044006
          10681.5102443695
          10643.7948532104
          10605.8639736176
          10567.7362804413
          10529.4298763275
          10490.9621810913
          10452.3502292633
          10413.6106853485
           10374.760055542
           10335.814786911
          10296.7914810181
          10257.7070178986
          10218.5787467957
          10179.4247169495
           10140.263961792
          10101.1167545319
          10062.0050106049
          10022.9526882172
          9983.98636054993
          9945.13582801819
          9906.43500137329
           9867.9229259491
          9829.64502334595
          9791.65490341187
          9754.01652336121
          9716.80701637268
          9680.12055587769
          9644.07351493835
          9608.81142807007
          9574.51874732971
          9541.43245124817
          9509.86200904846
          9480.21868133545
          9453.06057929993
          9429.16440391541
          9409.64512825012
          9396.16814422607
           9391.3547039032
          9399.63754272461
          9429.32519721985
          9498.61080360413
          9658.73331451416
           10143.057800293
          9573.65035438538
          9314.22665977478
          9153.30680274963
          9034.77348709106
          8939.18001747131
          8857.68670845032
           8785.5747127533
          8720.04758644104
          8659.32454872131
          8602.21113777161
          8547.87443351746
          8495.71584510803
          8445.29479217529
          8396.28052330017
          8348.42071533203
          8301.52013587952
          8255.42595100403
          8210.01714897156
          8165.19709968567
          8120.88796234131
          8077.02652931213
          8033.56115341187
          7990.44935035706
          7947.65602111816
          7905.15195655823
          7862.91282081604
          7820.91823768616
           7779.1510848999
          7737.59701538086
          7696.24391555786
          7655.08162689209
          7614.10165023804
          7573.29684829712
          7532.66136169434
          7492.19038200378
          7451.88005828857
          7411.72739982605
          7371.73017311096
          7331.88692092896
          7292.19679832459
          7252.65964698792
          7213.27596092224
          7174.04683303833
          7134.97403335571
          7096.05998420715
          7057.30781555176
          7018.72143936157
          6980.30554580688
          6942.06579017639
          6904.00885772705
          6866.14254188538
          6828.47603034973
          6791.02004623413
           6753.7871055603
          6716.79185676575
          6680.05146789551
          6643.58612442017
          6607.41961097717
          6571.58012771606
          6536.10124206543
          6501.02307701111
          6466.39394569397
          6432.27231788635
          6398.72953224182
          6365.85315322876
          6333.75164031982
          6302.56044578552
          6272.45044326782
           6243.6393699646
          6216.40808868408
          6191.12378501892
          6168.27402877808
          6148.51831054688
            6132.768699646
          6122.32202148438
          6119.08880615234
          6126.01847076416
          6147.96545410156
          6193.68287277222
          6281.27379417419
          6457.56125450134
           6910.1462726593
          7288.86568832397
           6447.6789150238
          6107.58288383484
           5904.2577419281
          5760.52418804169
          5649.19534111023
          5557.84543704987
          5479.82123756409
          5411.17102909088
          5349.37343406677
          5292.73460483551
          5240.07199668884
          5190.53516483307
          5143.49786281586
          5098.48987865448
          5055.15241718292
           5013.2077331543
          4972.43820858002
          4932.67141246796
          4893.76937580109
          4855.62070846558
          4818.13470172882
          4781.23694229126
          4744.86596488953
          4708.97069835663
          4673.50843429565
          4638.44329166412
          4603.74501419067
          4569.38797187805
          4535.35043334961
          4501.61391639709
          4468.16271305084
           4434.9834985733
           4402.0650100708
          4369.39779663086
          4336.97401237488
           4304.7872428894
          4272.83243846893
          4241.10577583313
          4209.60462474823
          4178.32750797272
          4147.27411746979
          4116.44532966614
          4085.84326171875
          4055.47138023376
          4025.33457374573
          3995.43941783905
          3965.79433727264
          3936.40990638733
          3907.29923439026
            3878.478474617
          3849.96740531921
          3821.79025173187
          3793.97669792175
          3766.56322479248
          3739.59485626221
           3713.1274394989
          3687.23072338104
          3661.99241733551
          3637.52379608154
          3613.96737766266
          3591.50765705109
          3570.38645458221
          3550.92518138886
          3533.55812931061
          3518.88381958008
           3507.7472486496
          3501.37852954865
          3501.64178848267
          3511.51990509033
          3536.16443538666
          3585.51429367065
          3682.25043964386
          3894.39727878571
          4569.02601623535
          4222.19100952148
          3647.53976821899
           3405.7001991272
          3258.40218448639
          3153.49999046326
          3072.09933185577
          3005.37285327911
           2948.5138502121
          2898.63126659393
          2853.86418914795
          2812.95549201965
          2775.02496147156
           2739.4388961792
          2705.73076152802
          2673.55061435699
          2642.63171958923
          2612.76786899567
          2583.79759788513
          2555.59297943115
          2528.05154418945
          2501.09038448334
          2474.64178848267
          2448.64989566803
          2423.06826019287
          2397.85789394379
          2372.98581504822
          2348.42388725281
          2324.14787101746
           2300.1367483139
          2276.37207984924
          2252.83758544922
          2229.51870059967
          2206.40227985382
          2183.47629547119
           2160.7295923233
          2138.15169429779
          2115.73256874084
          2093.46246433258
          2071.33170604706
          2049.33049964905
          2027.44875335693
          2005.67578029633
          1984.00009536743
          1962.40908622742
          1940.88861560822
          1919.42255973816
          1897.99220657349
          1876.57548904419
          1855.14600658417
          1833.67171955109
          1812.11320018768
           1790.4212884903
          1768.53392791748
          1746.37170124054
           1723.8316488266
          1700.77836513519
          1677.03110122681
          1652.34450531006
          1626.37925815582
           1598.6554851532
           1568.4758682251
          1534.79120254517
          1495.94725751877
           1449.1587228775
          1389.26208305359
          1305.17123842239
          1166.70206546783
          844.581665039063
          166.759293556213
          970.171051502228
          1199.45637226105
          1302.98317193985
          1358.56501960754
          1390.17698431015
          1407.81709289551
          1416.51918411255
          1419.15170812607
          1417.49208259583
          1412.70785570145
           1405.5990152359
            1396.731341362
          1386.51504516602
          1375.25365209579
          1363.17560815811
          1350.45543432236
           1337.2282333374
          1323.59981632233
          1309.65396642685
          1295.45768117905
          1281.06501722336
          1266.51994514465
          1251.85849380493
          1237.11037111282
          1222.30023097992
          1207.44860315323
          1192.57266283035
          1177.68679189682
           1162.8030667305
          1147.93158245087
          1133.08076000214
          1118.25757694244
          1103.46772289276
           1088.7157497406
          1074.00516748428
          1059.33853244781
          1044.71748018265
          1030.14274549484
          1015.61418247223
          1001.13071918488
          986.690320491791
          972.289879322052
          957.925097465515
          943.590319633484
          929.278277397156
          914.979766368866
          900.683197975159
          886.373961448669
          872.033590316772
          857.638541221619
            843.1584815979
          828.553746700287
          813.771569252014
           798.74022436142
          783.359568595886
          767.484952449799
           750.89789056778
           733.24746465683
          713.918109416962
          691.680055618286
          663.539134979248
          619.492636680603
          498.588088035584
          450.190250873566
          578.302042484283
          601.225942134857
          602.268917560577
          594.571352958679
          581.744557380676
          564.918721199036
           544.19344329834
          518.961398124695
          487.845976352692
           448.36629486084
          396.253502845764
          324.221385478973
          219.903778553009
          61.6041145324707
         -222.153700351715
          6365.85315275192];
end

function [ t ] = tmp3345_1000()
t=[6384.987159729
          3334.41683959961
          19566.3413085938
          21501.2367477417
          22058.2024497986
          22713.9204368591
          25825.8963279724
          26968.8798332214
          22404.2151107788
          20961.2570228577
          20003.2338256836
          19278.8180351257
          18918.6730651855
          19755.7245521545
           28904.518447876
          20010.9447669983
           13538.800327301
           11369.873249054
          9901.12171554565
           8532.7559967041
          6553.14596557617
         -2397.09487533569
            6384.987159729
          7403.88428878784
           7336.1847820282
           7007.3120174408
          6633.24347877502
          6579.75969314575
          5627.21812057495
          5085.56938362122
          4642.66779518127
          4263.84671401978
           4119.9414730072
          3709.89892101288
          2989.90361309052
           2620.6066493988
          2329.29304695129
          2174.20138454437
          1962.94193553925
          1417.13361167908
          1152.98385620117
          899.786981582642
          489.738758087158
         -58.5798916816711
         -2031.76002550125
          6384.98715877533];
end

function [t] = tmp33150_1000()
t=[10731.9130325317
         -92260.8785171509
         -30345.7333068848
         -10060.2655639648
         -1293.06498718262
          3305.90059661865
          6068.66543579102
          7885.28453826904
          9151.94203186035
          10069.7653961182
          10752.4349899292
          11269.7065658569
          11667.0596084595
          11975.4793777466
          12216.7690429688
          12406.6478347778
          12556.6624603271
          12675.4247894287
          12769.4375915527
          12843.6599197388
          12901.9070892334
           12947.138458252
          12981.6644821167
          13007.2973937988
          13025.4676895142
          13037.3092422485
          13043.7247772217
          13045.4347381592
          13043.0166320801
          13036.9347991943
          13027.5644683838
          13015.2091522217
          13000.1141281128
          12982.4791030884
          12962.4668121338
          12940.2099838257
           12915.817565918
          12889.3787765503
          12860.9671630859
          12830.6437149048
          12798.4587783813
          12764.4543914795
          12728.6656723022
          12691.1222305298
          12651.8491821289
          12610.8680000305
          12568.1971931458
          12523.8529129028
          12477.8494262695
          12430.1994628906
          12380.9145889282
          12330.0054244995
          12277.4819030762
          12223.3534088135
          12167.6289558411
          12110.3172950745
          12051.4270591736
          11990.9667739868
          11928.9449768066
          11865.3702926636
          11800.2514762878
            11733.59740448
          11665.4172096252
           11595.720199585
          11524.5159645081
          11451.8144187927
          11377.6257362366
           11301.960395813
          11224.8292503357
          11146.2434692383
          11066.2146339417
          10984.7546424866
          10901.8758468628
          10817.5909385681
          10731.9130363464
          10644.8556938171
          10556.4328422546
          10466.6589202881
          10375.5487480164
           10283.117603302
          10189.3812522888
          10094.3559379578
          9998.05832290649
          9900.50564575195
          9801.71556472778
          9701.70626068115
          9600.49647140503
          9498.10539245605
          9394.55284500122
          9289.85911560059
          9184.04509162903
          9077.13223648071
          8969.14253425598
          8860.09864234924
          8750.02376747131
          8638.94176673889
          8526.87710380554
          8413.85489654541
          8299.90092277527
          8185.04158401489
          8069.30402374268
          7952.71597671509
          7835.30598831177
          7717.10325431824
           7598.1376914978
          7478.43997001648
            7358.041431427
          7236.97419929504
           7115.2710647583
           6992.9655380249
          6870.09189224243
          6746.68486022949
          6622.78001785278
          6498.41335296631
          6373.62146568298
          6248.44128417969
          6122.90991783142
          5997.06462860107
          5870.94254684448
          5744.58041000366
          5618.01420211792
          5491.27868843079
          5364.40690994263
          5237.42940235138
          5110.37321662903
          4983.26074123383
          4856.10816669464
          4728.92351436615
           4601.7039642334
          4474.43244457245
          4347.07307243347
          4219.56527996063
          4091.81563854218
          3963.68683052063
          3834.98265552521
          3705.42693901062
          3574.63320922852
          3442.06056976318
          3306.94762706757
          3168.20955467224
           3024.2717590332
          2872.78714466095
          2710.12586307526
          2530.38751411438
           2323.3215303421
          2069.52259922028
           1728.0211482048
          1199.58979129791
          196.002087593079
         -2396.32327985764
          10731.9130349159];
end

function [t] = tmp33100_1000()
t=[9774.25926971436
         -28045.7534790039
          658.018371582031
          6514.32192993164
          8593.59822845459
          9648.82284545898
           10287.755317688
          10706.7679824829
          10992.8926467896
           11193.063293457
          11335.3257446289
          11437.3732452393
          11510.7461700439
           11563.142578125
          11599.7892990112
          11624.2916641235
          11639.1606445313
          11646.1661224365
          11646.5650558472
          11641.2543106079
          11630.8802185059
          11615.9133148193
          11596.6927490234
          11573.4669952393
          11546.4176101685
          11515.6781845093
          11481.3477783203
           11443.499671936
          11402.1887435913
          11357.4564743042
           11309.334815979
          11257.8488235474
          11203.0188293457
          11144.8619270325
          11083.3931312561
           11018.626335144
          10950.5748329163
          10879.2519607544
          10804.6714096069
          10726.8475341797
          10645.7956466675
          10561.5321922302
          10474.0748825073
          10383.4428634644
          10289.6568336487
          10192.7391014099
          10092.7136764526
          9989.60641860962
          9883.44505691528
          9774.25926971436
          9662.08074569702
          9546.94329452515
          9428.88300704956
          9307.93801116943
          9184.14900588989
          9057.55893325806
          8928.21332550049
          8796.16022491455
          8661.45036697388
          8524.13731002808
           8384.2773475647
           8241.9298324585
          8097.15703964233
          7950.02451133728
          7800.60098838806
          7648.95856285095
          7495.17269325256
          7339.32237434387
          7181.49018669128
           7021.7622718811
          6860.22838973999
          6696.98192024231
          6532.11960220337
          6365.74139595032
          6197.95004272461
          6028.85055160522
          5858.54933357239
          5687.15301704407
          5514.76620483398
          5341.48929405212
          5167.41464614868
          4992.62115859985
          4817.16694164276
          4641.07797241211
          4464.33230495453
          4286.83711338043
          4108.39417648315
          3928.64802646637
          3747.00707149506
           3562.5187330246
          3373.66857624054
          3178.04226875305
          2971.72785282135
          2748.19042110443
          2495.96167469025
          2193.33870029449
          1794.32859325409
          1183.87468719482
         -3.48129463195801
         -3447.57680749893
          9774.25926971436];
end

function [ w ] = tmp33240_1000()
w=[12852.0129585266
         -281698.442531586
         -189511.249755859
         -127453.139732361
         -87093.0190620422
         -59847.9290962219
         -40761.3816452026
         -26955.3270454407
          -16690.091381073
         -8872.22232437134
         -2791.01805496216
          2028.84714126587
          5913.20347213745
          9090.34608078003
          11723.4305915833
          13931.2292556763
           15801.712928772
          17401.1449012756
          18780.2739334106
          19978.6308746338
          21027.5673751831
          21952.4545707703
          22774.3226890564
          23511.1363487244
          24178.8625221252
          24792.4609794617
          25366.9416160584
          25918.6819763184
          26467.3185005188
          27038.8173942566
          27670.9892234802
          28424.4748916626
          29407.5757026672
          30843.6616439819
           33319.212852478
          39397.3568878174
          34690.6937103271
          31907.4104385376
          30543.0940246582
          29724.8350830078
          29179.4195747375
          28786.6582794189
          28484.7501106262
          28238.9565696716
          28028.5936775208
          27840.8471412659
          27667.5533332825
          27503.4113464355
          27344.9407348633
          27189.8496208191
          27036.6410331726
          26884.3613853455
          26732.4381904602
          26580.5735855103
          26428.6753005981
          26276.8134765625
          26125.1956100464
          25974.1563491821
          25824.1590614319
          25675.8103790283
          25529.8891944885
          25387.3927612305
          25249.6079216003
          25118.2183074951
          24995.4660682678
           24884.398349762
          24789.2510604858
          24716.0609893799
          24673.6719017029
          24675.4517021179
          24742.3643493652
          24908.8086738586
          25234.6694908142
          25833.2697181702
          26948.9545173645
          29249.7931556702
           35845.395942688
          30397.1553688049
          26327.2106971741
           24357.467792511
          23136.2052574158
          22282.6651000977
          21636.3998069763
          21116.8075752258
          20679.1333465576
          20296.8154983521
          19953.2997932434
          19637.8629570007
          19343.3310241699
          19064.7680835724
          18798.6902942657
          18542.5802173615
          18294.5768127441
          18053.2731342316
          17817.5799827576
          17586.6323184967
          17359.7222118378
          17136.2484397888
          16915.6764774323
          16697.5024604797
          16481.2174415588
          16266.2645092011
          16051.9826068878
          15837.5228157043
          15621.7143783569
          15402.8362102509
           15178.202299118
          14943.3565330505
          14690.3655185699
           14403.729101181
          14048.5648422241
          13523.5734119415
          12300.9072399139
          13179.4621219635
          13402.0091457367
          13421.1008338928
          13359.5559024811
          13257.8821754456
          13133.9780521393
          12996.9864711761
          12852.0129585266
           12702.090265274
          12549.1059131622
          12394.2790184021
          12238.4212207794
           12082.086101532
           11925.657989502
          11769.4063339233
          11613.5192890167
          11458.1244754791
          11303.3014831543
          11149.0881958008
           10995.481508255
          10842.4328651428
          10689.8367900848
           10537.508682251
          10385.1449604034
          10232.2510948181
          10078.0076770782
           9921.0063419342
          9758.68676567078
          9585.98136520386
          9391.37822914124
          9141.07387924194
          8654.11965179443
          8761.63831520081
          8846.99244499207
           8803.0346736908
          8719.93364524841
           8619.5057888031
          8510.34871101379
          8396.62931060791
          8280.62519454956
          8163.69611740112
          8046.72121810913
          7930.31862258911
          7814.96675872803
          7701.07833862305
           7589.0515460968
          7479.31279563904
          7372.36010742188
          7268.81570625305
          7169.49871444702
          7075.53389549255
          6988.52400970459
          6910.83512687683
          6846.08975791931
          6800.06115531921
          6782.38871574402
          6810.12779808044
          6915.93641281128
          7170.48468589783
          7765.38715553284
          9594.62506866455
          11108.1646633148
          7329.38069438934
          6107.15048503876
          5468.04528713226
          5064.21235752106
           4778.3440694809
          4559.04200744629
          4380.43989181519
          4228.30378437042
          4094.26153945923
          3973.11202907562
          3861.46524620056
          3757.01057910919
          3658.10119915009
          3563.50445652008
          3472.24170017242
          3383.47424983978
          3296.40821456909
          3210.19498729706
          3123.79949569702
          3035.78952217102
          2943.95243549347
          2844.52548694611
          2730.49209690094
          2587.30306148529
          2379.74328804016
            1993.715883255
          670.335522651672
          1742.19987392426
          2256.43843364716
          2421.45724010468
          2479.84534645081
          2489.90452289581
          2474.47892284393
          2444.45426845551
          2405.61053276062
          2361.25047683716
          2313.37131977081
          2263.23856067657
          2211.68681240082
          2159.28674411774
          2106.44205474854
          2053.44942951202
          2000.53854656219
          1947.90250587463
          1895.72730350494
          1844.23021221161
          1793.72414779663
          1744.74819803238
          1698.37376213074
          1657.06845664978
          1628.07951593399
          1650.92876338959
          1697.64212942123
          1459.10256576538
          1356.08037662506
          1277.39448738098
           1205.5604891777
          1133.77845478058
          1057.23715400696
          970.303649425507
          863.728557109833
          718.829514980316
          491.240943431854
          55.7967891693115
         -1037.36587667465
          12852.0129590034];
end