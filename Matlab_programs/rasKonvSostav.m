function koTeplop = rasKonvSostav()
format long g;
r=rasPronits()
koTepl=rasKonSos207(r);
%koTepl=rasKonSos2073();
%koTepl=rasKonSos84();
%koTepl=rasKonSos842();
koTeplop=0;
end

function koef = koefoslab(wmg, wsi, wal, T)
tere=3e2:1e2:16e2; 
mgo=tmp49(0,tere);
alo=tmp49(2,tere);
sio=tmp49(1,tere);
rt=length(tere);
    f=1; j=rt; 
    for q=1:rt
if ((tere(q)>T) && (f>0))
        j=q; f=0; break;
end
    end
    if ((f == 0) && (j==1))
    j=2;
    end
if ((f == 0) && (j>1))
    koefm=mgo(j-1)+(mgo(j)-mgo(j-1))*(T-tere(j-1))/(tere(j)-tere(j-1)); 
    koefa=alo(j-1)+(alo(j)-alo(j-1))*(T-tere(j-1))/(tere(j)-tere(j-1));
    koefs=sio(j-1)+(sio(j)-sio(j-1))*(T-tere(j-1))/(tere(j)-tere(j-1)); 
end
wo=wmg+wsi+wal;
kuo=(koefm*wmg+koefa*wal+koefs*wsi)/wo;
koef=kuo;
end

function t = rasPronits()
vfv=0;
srpv=srRaPorVerm(0);
%srps=SreRazPorSha()*1e-6;
%lam1=tmp55();
r=srpv;
k1p=r^2/32;
if (vfv==0) 
    por = 0.6635; %пористость вермикулита фракции 2-0,7 мм
    r = (2+0.7)*1e-3/2;
elseif (vfv==1) 
    por=0.5575; %пористость вермикулита фракции 8-4 мм
    r = (8+4)*1e-3/2;
end
d=2*r;
%npp=Kramers_n_uk(); te0=273.15;
%tem = (2e2:1e2:12e2)+te0;
%n=length(tem); %wsio=368e-3; walo=132e-3; wmgo=217e-3;
%for k=1:n
%kosc(k)=koefoslab(wmg,wsi,wal,ts); eps(k)=epsisred(tem(k),npp); eps(k)=eps(k)*kosc(k);
%lavo(k)=opredTeploprovVozd(tem(k));
%end
k1p=r^2/32;
c=opredelc(1-por);
delt=d*((1/c-1)^(-1));
kp=k1p*(c^2);
%for k=1:n
%rpd(k)=DulnevSloFor(por,k1p);
%end
p=DulnevSloFor(por,k1p,kp);
t=p;
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

function t = opredelc(m1)
nko=1e5; lada=-1e1; ladb=1e5; ra=1e2; ep=1e-6; h=0;
while ((ra>ep) && (h<nko))
ladc=(lada+ladb)/2;
fa=2*(lada^3)-3*(lada^2)+m1;
fb=2*(ladb^3)-3*(ladb^2)+m1;
fc=2*(ladc^3)-3*(ladc^2)+m1;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); 
ra=abs(fa-fb);
h=h+1;
end
t=ladc;
end

function SRP = srRaPorVerm(vfv)
legr = LevGranPorom();
prgr = PravGranPorom();
sr = (legr + prgr) / 2; %средний размер пор
if (vfv==1)
rp = rasPorpoRazm84(); %только для фракции 8-4 мм
elseif (vfv==0) 
    rp = rasPorpoRazm207(); %только для фракции 2-0,7 мм
end
s = 0;
for k=1:length(rp)
    s = s + sr(k) * rp(k);
end
SRP = s;
end

function srp = SreRazPorSha()
rpr(1) = 0;
raspr=[21.8,21.75,21.625,21.25,20.75,20.25,19.75,19,18.25,17.5,12,9.875,8.75,8.25,7.25,6.25,5,4.25,3.25,0.875,0.325,0];
for k=2:length(raspr)
    rpr(k)=(raspr(k-1)-raspr(k))*1e2/raspr(1);
end
prgr=[0.05,0.1,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50]*1e-6; %в метрах
legr(1)=0;
for k=2:length(prgr)
    legr(k)=prgr(k-1);
end
sr = (legr + prgr) / 2; %средний размер пор
le = length(sr);
srpo = 0;
for k=1:le
    if (rpr(k) > 0)
    srpo = srpo + rpr(k) * sr(k) / 1e2; %объемная доля поры заданного размера в полном (во всем) объеме, в процентах
    end
end
srp = srpo;
end

function lam4 = DulnevSloFor(m2,k1p,kp)
Nk=sqrt(m2^2-10*m2+9); Nk=(m2+3+Nk)/2/m2; %y2=3.3e-3*(1-m2)^(-2/9); %pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2); %y2=y2*(pud+9.8*ro1*(1-m2)*hsl)^(1/3); %y1=1e-20; y2=y1;
eta=1e-4; y1=10e-4; y2=y1/sqrt(eta); %y1=(10+50)*1e-4/2;
y3=2*sqrt(Nk-1)/Nk;
y4=y3/((1-m2)^(1/3));
kp=kp*(y4^2)/(y4^2-y3^2);
D=sqrt(1-y3^2);
E=y4^2-y3^2;
nu=kp/k1p;
w=(1-nu*D)/(1-nu);
DF=abs(w-D)/abs(w-1);
DF=D-1-w*log(DF); DF=1/DF;
DF=(1-nu)*DF/2/nu;
DF=DF+D/(y3^2); DF=1/DF;
lam4=DF*k1p/(y4^2)+kp*E/(y4^2);
end

function vyb = rasKonSos207(ka)
%Фракция 2-0.7 мм
y0=30e-3;
l=1e-3;
x=0:l:y0;
te0=273.15;
por=0.6635;
g=9.81;
ro=25e1;
specarea=3.377e3;
apzeos=ro*specarea;
deq=4*por/apzeos
temvh=arrTemHigh207();
temvh=temvh+te0; 
temvc=arrTemCold207(); 
temvc=temvc+te0; 
temvs=(temvh+temvc)/2;
tepo=arrTepPot207();
qv=(1e4/13.85)*tepo; 
tepv=koeftepv(qv,y0,temvh,temvc,length(temvs));
vtem=5;
%a=lambda/cp/ro;
Thi=arrTem1VVF1();
Thi=Thi+te0;
Th=Thi(vtem);
Tco=arrTem2VVF1();
Tco=Tco+te0;
Tc=Tco(vtem);
Tsred=arrTem3VVF1();
Tsred=Tsred+te0;
Ts=Tsred(vtem);
Tsr=(Thi+Tco)/2;
Tsredn=(Th+Tc)/2;
ktp=arrKTP_VVF1();
lambda=ktp(vtem);
%udte=arrUdTepl_VVF1();
%koudte=polyfit(Tsr,udte,2);
%ute=polyval(koudte,Tsredn);
%T=nachPribTem(Th);
%T=pribltemkTc(T,Th,Tc,x);
Tm=NachPribN1(Th,Tc,Ts,x);
Th=Tm(1);
Nt=length(Tm);
ktpm=polyfit([Tsr(1) Tsr(4) Tsr(5)],[ktp(1) ktp(4) ktp(5)],2);
koTep=0; koTe=0; ktpr=0; Tsre=0; Gras=0; kri=0; srch=(2+0.7)*1e-3/2;
for k=1:Nt
    Tsre(k)=(Tm(k)+Th)/2;
    beta0=1/(Tsre(k));
    ktpr(k)=polyval(ktpm,Tsre(k));
PrVo=vychKoefVoz(1,Tsre(k));
roVoz=vychKoefVoz(2,Tsre(k));
muVoz=vychKoefVoz(3,Tsre(k));
lambVoz=vychKoefVoz(4,Tsre(k));
cp=vychKoefVoz(5,Tsre(k));
atem=lambVoz/cp/roVoz;
nuVo=muVoz/roVoz;
kri(k)=g*beta0*(Th-Tm(k))*srch*ka/atem/nuVo;
Gras(k)=beta0*(Th-Tm(k))*g*(deq^3)*(roVoz^2)/(muVoz^2);
psi1(k)=1+2*beta0*(Tm(k)-Th)*g*(roVoz^2)*(por^3)*cp*x(k)/18/(apzeos^2)/muVoz/ktpr(k); %psi1
psi2(k)=1-(16e-5)*PrVo*Gras(k)*((por*x(k))^2)*lambVoz/(deq^2)/ktpr(k); %psi2
end
krim=max(kri)
Tsre=Tsre';
ktpr=ktpr';
kot1=0; kot2=0; psi11=0; psi12=0;
for k=1:Nt
    kot1(k)=(psi1(k)-1)*ktpr(k);
    kot2(k)=(psi2(k)-1)*ktpr(k);
    psi11(k)=psi1(k)-1;
    psi12(k)=psi2(k)-1;
end
psi3=0; psi3=(psi11+psi12)/2;
kot3=0; kot3=(kot1+kot2)/2;
%beta0=vyvod(ktpr);
disp('Tm');
beta0=vyvod(Tm);
disp('Gras');
beta0=vyvod(Gras);
disp('xi');
beta0=vyvod(x);
%beta0=vyvod(kot3);
%beta0=vyvod(psi3);
disp('delta1');
beta0=vyvod(kot1);
disp('delta2');
beta0=vyvod(kot2);
disp('psi1');
beta0=vyvod(psi11);
disp('psi2');
beta0=vyvod(psi12);
Tsredn=(Tc+Th)/2;
beta0=1/(Tsredn);
PrVo=vychKoefVoz(1,Tsredn);
roVoz=vychKoefVoz(2,Tsredn);
muVoz=vychKoefVoz(3,Tsredn);
a=1e20;
b=1e-2;
e=1e-5;
p=1;
while (abs(a-b)>e)
    c=(a+b)/2;
Grasa=beta0*a*g*(deq^3)*(roVoz^2)/(muVoz^2);
Grasb=beta0*b*g*(deq^3)*(roVoz^2)/(muVoz^2);
Grasc=beta0*c*g*(deq^3)*(roVoz^2)/(muVoz^2);
Raa=Grasa*PrVo-680;
Rab=Grasb*PrVo-680;
Rac=Grasc*PrVo-680;
if (Raa*Rac<0)
    b=c;
end
if (Rab*Rac<0)
    a=c;
end
if (p>1e3)
    break;
end
p=p+1;
end
c=1*c;
%pl=plot(1e3*x,koTep,'-b',1e3*x,koTe,'-k');
%pl=plot(1e3*x,koTep,'-b');
%set(pl,'LineWidth',2); 
%hold on; 
%grid on; 
%xlabel({'Координата, мм'}); 
%ylabel({'Коэффициент теплопроводности без конвекции, Вт/(м*К)'}); 
%title({'График зависимости температуры от координаты'});
vyb=0;
end

function vyb = rasKonSos84()
%Фракция 8-4 мм
y0=30e-3;
l=1e-3;
x=0:l:y0;
te0=273.15;
por=0.5575;
g=9.81;
ro=2e2;
specarea=3.578e3;
apzeos=ro*specarea;
deq=4*por/apzeos
temvh=arrTemHigh(); 
temvh=temvh+te0; 
temvc=arrTemCold(); 
temvc=temvc+te0; 
temvs=(temvh+temvc)/2;
tepo=arrTepPot84(); 
qv=(1e4/13.85)*tepo; 
tepv=koeftepv(qv,y0,temvh,temvc,length(temvs));
vtem=5;
ktpm=danPoTemTepl(temvs,tepv);
kTh=danPoTemH(temvs,temvh);
tems1=(temvs(16)+temvs(18))/2;
%Th=polyval(kTh,tems1)
%Th=(temvh(15)+temvh(17))/2;
Th=(temvh(16)+temvh(18))/2;
kTc=danPoTemC(temvs,temvc);
%Tc=polyval(kTc,tems1)
%Tc=(temvc(15)+temvc(17))/2;
Tc=(temvc(16)+temvc(18))/2;
koTe=polyval(ktpm,tems1)*(Th-Tc)/y0;
Tsr=(Th+Tc)/2;
T=nachPribTem(Th);
Tm=pribltemkTc(T,Th,Tc,x);
Tm=Tm';
Nt=length(Tm);
%udte=arrUdTepl_VVF1();
%koudte=polyfit(Tsr,udte,2);
koTep=0; koTe=0; ktpr=0; Gras=0; Tsre=0;
for k=1:Nt
    Tsre(k)=(Tm(k)+Th)/2;
    beta0=1/(Tsre(k));
    ktpr(k)=polyval(ktpm,Tsre(k));
PrVo=vychKoefVoz(1,Tsre(k));
roVoz=vychKoefVoz(2,Tsre(k));
muVoz=vychKoefVoz(3,Tsre(k));
lambVoz=vychKoefVoz(4,Tsre(k));
cp=vychKoefVoz(5,Tsre(k));
Gras(k)=beta0*(Th-Tm(k))*g*(deq^3)*(roVoz^2)/(muVoz^2);
psi1(k)=1+2*beta0*(Tm(k)-Th)*g*(roVoz^2)*(por^3)*cp*x(k)/18/(apzeos^2)/muVoz/ktpr(k); %psi1
psi2(k)=1-(16e-5)*PrVo*Gras(k)*((por*x(k))^2)*lambVoz/(deq^2)/ktpr(k); %psi2
end
Tsre=Tsre';
ktpr=ktpr';
kot1=0; kot2=0; psi11=0; psi12=0;
for k=1:Nt
    kot1(k)=(psi1(k)-1)*ktpr(k);
    kot2(k)=(psi2(k)-1)*ktpr(k);
    psi11(k)=psi1(k)-1;
    psi12(k)=psi2(k)-1;
end
psi13=0; psi13=(psi12+psi11)/2;
kot3=0; kot3=(kot1+kot2)/2;
%beta0=vyvod(Tm);
%beta0=vyvod(Tsre);
%beta0=vyvod(ktpr);
disp('Tm');
beta0=vyvod(Tm);
disp('Gras');
beta0=vyvod(Gras);
disp('xi');
beta0=vyvod(x);
%beta0=vyvod(kot3);
%beta0=vyvod(psi3);
disp('delta1');
beta0=vyvod(kot1);
disp('delta2');
beta0=vyvod(kot2);
disp('psi1');
beta0=vyvod(psi11);
disp('psi2');
beta0=vyvod(psi12);
Tsre=(Tc+Th)/2;
beta0=1/(Tsre);
PrVo=vychKoefVoz(1,Tsre);
roVoz=vychKoefVoz(2,Tsre);
muVoz=vychKoefVoz(3,Tsre);
a=1e20;
b=1e-2;
e=1e-5;
p=1;
while (abs(a-b)>e)
    c=(a+b)/2;
Grasa=beta0*a*g*(deq^3)*(roVoz^2)/(muVoz^2);
Grasb=beta0*b*g*(deq^3)*(roVoz^2)/(muVoz^2);
Grasc=beta0*c*g*(deq^3)*(roVoz^2)/(muVoz^2);
Raa=Grasa*PrVo-680;
Rab=Grasb*PrVo-680;
Rac=Grasc*PrVo-680;
if (Raa*Rac<0)
    b=c;
end
if (Rab*Rac<0)
    a=c;
end
if (p>1e3)
    break;
end
p=p+1;
end
c=1*c;
%pl=plot(1e3*x,koTep,'-b',1e3*x,koTe,'-k');
%pl=plot(1e3*x,koTep,'-b');
%set(pl,'LineWidth',2); 
%hold on; 
%grid on; 
%xlabel({'Координата, мм'}); 
%ylabel({'Коэффициент теплопроводности без конвекции, Вт/(м*К)'}); 
%title({'График зависимости температуры от координаты'});
vyb=0;
end

function t = vyvod(arr)
Nt=length(arr);
ar=0;
q=1;
for k=1:Nt
    if (rem(k-1,5)==0)
        ar(q)=arr(k);
        q=q+1;
    end
end
ar=ar'
t=0;
end

function vyb = rasKonSos2073()
%Фракция 2-0.7 мм
y0=30e-3;
l=1e-3;
x=0:l:y0;
te0=273.15;
por=0.6635;
g=9.81;
ro=25e1;
specarea=3.377e3;
apzeos=ro*specarea;
deq=4*por/apzeos
temvh=arrTemHigh207();
temvh=temvh+te0; 
temvc=arrTemCold207(); 
temvc=temvc+te0; 
temvs=(temvh+temvc)/2;
tepo=arrTepPot207();
qv=(1e4/13.85)*tepo; 
tepv=koeftepv(qv,y0,temvh,temvc,length(temvs));
ktpm=danPoTemTepl2073(temvs,tepv);
kTc=danPoTemC2073(temvs,temvc);
kTh=danPoTemH2073(temvs,temvh);
tems1=(temvs(8)+temvs(10))/4;
%Th=polyval(kTh,tems1)
%Th=(temvh(15)+temvh(17))/2;
Th=(temvh(8)+temvh(10))/2;
%Tc=polyval(kTc,tems1)
%Tc=(temvc(15)+temvc(17))/2;
Tc=(temvc(8)+temvc(10))/2;
Tsr=(Th+Tc)/2;
T=nachPribTem(Th);
Tm=pribltemkTc(T,Th,Tc,x);
Th=Tm(1);
Nt=length(Tm);
koTep=0; koTe=0; ktpr=0; Tsre=0; Gras=0;
for k=1:Nt
    Tsre(k)=(Tm(k)+Th)/2;
    beta0=1/(Tsre(k));
    ktpr(k)=polyval(ktpm,Tsre(k));
PrVo=vychKoefVoz(1,Tsre(k));
roVoz=vychKoefVoz(2,Tsre(k));
muVoz=vychKoefVoz(3,Tsre(k));
lambVoz=vychKoefVoz(4,Tsre(k));
cp=vychKoefVoz(5,Tsre(k));
Gras(k)=beta0*(Th-Tm(k))*g*(deq^3)*(roVoz^2)/(muVoz^2);
psi1(k)=1+2*beta0*(Tm(k)-Th)*g*(roVoz^2)*(por^3)*cp*x(k)/18/(apzeos^2)/muVoz/ktpr(k); %psi1
psi2(k)=1-(16e-5)*PrVo*Gras(k)*((por*x(k))^2)*lambVoz/(deq^2)/ktpr(k); %psi2
end
Tsre=Tsre';
ktpr=ktpr';
kot1=0; kot2=0; psi11=0; psi12=0;
for k=1:Nt
    kot1(k)=(psi1(k)-1)*ktpr(k);
    kot2(k)=(psi2(k)-1)*ktpr(k);
    psi11(k)=psi1(k)-1;
    psi12(k)=psi2(k)-1;
end
psi3=0; psi3=(psi11+psi12)/2;
kot3=0; kot3=(kot1+kot2)/2;
disp('Tm');
beta0=vyvod(Tm);
disp('ktp');
beta0=vyvod(ktpr);
disp('Gras');
beta0=vyvod(Gras);
disp('xi');
beta0=vyvod(x);
%beta0=vyvod(kot3);
%beta0=vyvod(psi3);
disp('delta1');
beta0=vyvod(kot1);
disp('delta2');
beta0=vyvod(kot2);
disp('psi1');
beta0=vyvod(psi11);
disp('psi2');
beta0=vyvod(psi12);
Tsredn=(Tc+Th)/2;
beta0=1/(Tsredn);
PrVo=vychKoefVoz(1,Tsredn);
roVoz=vychKoefVoz(2,Tsredn);
muVoz=vychKoefVoz(3,Tsredn);
a=1e20;
b=1e-2;
e=1e-5;
p=1;
while (abs(a-b)>e)
    c=(a+b)/2;
Grasa=beta0*a*g*(deq^3)*(roVoz^2)/(muVoz^2);
Grasb=beta0*b*g*(deq^3)*(roVoz^2)/(muVoz^2);
Grasc=beta0*c*g*(deq^3)*(roVoz^2)/(muVoz^2);
Raa=Grasa*PrVo-680;
Rab=Grasb*PrVo-680;
Rac=Grasc*PrVo-680;
if (Raa*Rac<0)
    b=c;
end
if (Rab*Rac<0)
    a=c;
end
if (p>1e3)
    break;
end
p=p+1;
end
c=1*c;
%pl=plot(1e3*x,koTep,'-b',1e3*x,koTe,'-k');
%pl=plot(1e3*x,koTep,'-b');
%set(pl,'LineWidth',2); 
%hold on; 
%grid on; 
%xlabel({'Координата, мм'}); 
%ylabel({'Коэффициент теплопроводности без конвекции, Вт/(м*К)'}); 
%title({'График зависимости температуры от координаты'});
vyb=0;
end

function vyb = rasKonSos842()
%Фракция 8-4 мм
y0=30e-3;
l=1e-3;
x=0:l:y0;
te0=273.15;
por=0.5575;
g=9.81;
ro=2e2;
specarea=3.578e3;
apzeos=ro*specarea;
deq=4*por/apzeos
temvh=arrTemHigh(); 
temvh=temvh+te0; 
temvc=arrTemCold(); 
temvc=temvc+te0; 
temvs=(temvh+temvc)/2;
tepo=arrTepPot84(); 
qv=(1e4/13.85)*tepo; 
tepv=koeftepv(qv,y0,temvh,temvc,length(temvs));
vtem=5;
ktpm=danPoTemTepl2(temvs,tepv);
kTh=danPoTemH2(temvs,temvh);
tems1=(temvs(6)+temvs(12)+temvs(8))/3;
%Th=polyval(kTh,tems1)
%Th=(temvh(15)+temvh(17))/2;
Th=(temvh(6)+temvh(8)+temvh(12))/3;
kTc=danPoTemC2(temvs,temvc);
%Tc=polyval(kTc,tems1)
%Tc=(temvc(15)+temvc(17))/2;
Tc=(temvc(6)+temvc(8)+temvc(12))/3;
Tsr=(Th+Tc)/2;
T=nachPribTem(Th);
Tm=pribltemkTc(T,Th,Tc,x);
Tm=Tm';
Nt=length(Tm);
%udte=arrUdTepl_VVF1();
%koudte=polyfit(Tsr,udte,2);
koTep=0; koTe=0; ktpr=0; Gras=0; Tsre=0;
for k=1:Nt
    Tsre(k)=(Tm(k)+Th)/2;
    beta0=1/(Tsre(k));
    ktpr(k)=polyval(ktpm,Tsre(k));
PrVo=vychKoefVoz(1,Tsre(k));
roVoz=vychKoefVoz(2,Tsre(k));
muVoz=vychKoefVoz(3,Tsre(k));
lambVoz=vychKoefVoz(4,Tsre(k));
cp=vychKoefVoz(5,Tsre(k));
Gras(k)=beta0*(Th-Tm(k))*g*(deq^3)*(roVoz^2)/(muVoz^2);
psi1(k)=1+2*beta0*(Tm(k)-Th)*g*(roVoz^2)*(por^3)*cp*x(k)/18/(apzeos^2)/muVoz/ktpr(k); %psi1
psi2(k)=1-(16e-5)*PrVo*Gras(k)*((por*x(k))^2)*lambVoz/(deq^2)/ktpr(k); %psi2
end
Tsre=Tsre';
ktpr=ktpr';
kot1=0; kot2=0; psi11=0; psi12=0;
for k=1:Nt
    kot1(k)=(psi1(k)-1)*ktpr(k);
    kot2(k)=(psi2(k)-1)*ktpr(k);
    psi11(k)=psi1(k)-1;
    psi12(k)=psi2(k)-1;
end
psi13=0; psi13=(psi12+psi11)/2;
kot3=0; kot3=(kot1+kot2)/2;
%beta0=vyvod(Tm);
%beta0=vyvod(Tsre);
%beta0=vyvod(ktpr);
disp('Tm');
beta0=vyvod(Tm);
disp('Gras');
beta0=vyvod(Gras);
disp('xi');
beta0=vyvod(x);
%beta0=vyvod(kot3);
%beta0=vyvod(psi3);
disp('delta1');
beta0=vyvod(kot1);
disp('delta2');
beta0=vyvod(kot2);
disp('psi1');
beta0=vyvod(psi11);
disp('psi2');
beta0=vyvod(psi12);
Tsre=(Tc+Th)/2;
beta0=1/(Tsre);
PrVo=vychKoefVoz(1,Tsre);
roVoz=vychKoefVoz(2,Tsre);
muVoz=vychKoefVoz(3,Tsre);
a=1e20;
b=1e-2;
e=1e-5;
p=1;
while (abs(a-b)>e)
    c=(a+b)/2;
Grasa=beta0*a*g*(deq^3)*(roVoz^2)/(muVoz^2);
Grasb=beta0*b*g*(deq^3)*(roVoz^2)/(muVoz^2);
Grasc=beta0*c*g*(deq^3)*(roVoz^2)/(muVoz^2);
Raa=Grasa*PrVo-680;
Rab=Grasb*PrVo-680;
Rac=Grasc*PrVo-680;
if (Raa*Rac<0)
    b=c;
end
if (Rab*Rac<0)
    a=c;
end
if (p>1e3)
    break;
end
p=p+1;
end
c=1*c;
%pl=plot(1e3*x,koTep,'-b',1e3*x,koTe,'-k');
%pl=plot(1e3*x,koTep,'-b');
%set(pl,'LineWidth',2); 
%hold on; 
%grid on; 
%xlabel({'Координата, мм'}); 
%ylabel({'Коэффициент теплопроводности без конвекции, Вт/(м*К)'}); 
%title({'График зависимости температуры от координаты'});
vyb=0;
end