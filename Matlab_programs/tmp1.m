%определяет ПТП, КТП, температуры для шамота
function t = tmp1()
format long g;
te0=273.15;
t1=arrTem1VVF1(); t1=t1+te0;
t2=arrTem2VVF1(); t2=t2+te0;
t3=arrTem3VVF1(); t3=t3+te0;
ktp=arrKTP_VVF1();
tol=3e1*1e-3;
dt=(t2-t1);
dtdx=dt/tol;
tsr=(t1+t2)/2;
ts=[(tsr(1)+tsr(2)+tsr(3))/3,tsr(4),(tsr(5)+tsr(6))/2];
tgo=[(t1(1)+t1(2)+t1(3))/3,t1(4),(t1(5)+t1(6))/2];
kotgo=polyfit(ts,tgo,2);
tho=[(t2(1)+t2(2)+t2(3))/3,t2(4),(t2(5)+t2(6))/2];
kotho=polyfit(ts,tho,2);
tce=[(t3(1)+t3(2)+t3(3))/3,t3(4),(t3(5)+t3(6))/2];
kotce=polyfit(ts,tce,2);
q=0;
for k=1:length(dt)
q(k)=-ktp(k)*dtdx(k);
end
koq=polyfit(tsr',q,2);
n=length(tce);
koktp=polyfit(tsr,ktp,2);
Nt=1e3;
h=tol/Nt;
x=0:h:tol;
p=length(x);
for l=1:n
tgor=tgo(l);
thol=tho(l);
tcen=tce(l);
kote=polyfit([0,tol/2,tol],[tgor,tcen,thol],2);
for k=1:p-1
    te(l,k)=polyval(kote,x(k));
    telk1=polyval(kote,x(k+1));
    dtdxn=(telk1-te(l,k))/(x(k+1)-x(k));
    tet=(telk1+te(l,k))/2;
    ktp(l,k)=polyval(koktp,tet);
    ktpn=ktp(l,k);
    qo(l,k)=-ktpn*dtdxn;
end
qob=mean(qo(l,:))
ktps=mean(ktp(l,:))
tesre=mean(te(l,:))
end
ts=round(ts');
t=postgrafik(ktp,ts,1e3*x(1:p-1),te,qo);
end
function t = postgrafik(ktp,ts,x,te,qo)
%pl=plot(te(1,:),qo(1,:),'-b',te(2,:),qo(2,:),'-k',te(3,:),qo(3,:),'-m');
set(gca,'FontName','Times New Roman','Fontsize',14);
pl=plot(x,te(1,:),'-b',x,te(2,:),'-k',x,te(3,:),'-m');
set(pl,'LineWidth',3); 
hold on; grid on; 
xlabel({'Координата, мм'}); 
%xlabel({'Температура, К'}); 
%ylabel({'ПТП, Вт/м2'});
ylabel({'Температура, К'}); 
%title({'График зависимости температуры от координаты'});
%title({'График зависимости ПТП от температуры'});
sts1=sprintf('%d',ts(1));
sts2=sprintf('%d',ts(2));
sts3=sprintf('%d',ts(3));
legend(sts1,sts2,sts3,'location','best');
t=0;
end
function vy = vychSham()
format long g;
te0=273.15;
Thde=1e2;
tem=0; tem=1e2:Thde:12e2; 
koef=[-0.435e-9,0.685e-6,0.134e-3,0.725]; 
koef=[-0.397e-9,0.71e-6,0.011e-3,0.851];
koef=[-0.377e-9,0.918e-6,-0.338e-3,0.77];
n=length(tem);
h=65e-3;
Thna=3e2-te0;
qob=[2723,4235,5796];
te=[633,773,907]-te0;
ep=1e-1; nit=1e6; d=1e-2;
koeq=polyfit(te,qob,length(te)-1);
for k=1:n
del=1;
j=1;
ts=tem(k);
laef(k)=polyval(koef,ts);
while ((del>ep) && (j<nit))
Th(k)=Thna+(k-1)*Thde+(j-1)*d;
qo(k)=polyval(koeq,ts);
Tg(k)=Th(k)+qo(k)*h/laef(k);
del=2*tem(k)-Tg(k)-Th(k);
j=j+1;
end
end
%disp(length(Th));
%disp(length(Tg));
Th=Th'
qo=qo'
Tg=Tg'
tes=(Th+Tg)/2
laef=laef'
for k=1:n
    kote(k)=qo(k)*h/(Tg(k)-Th(k));
end
kote=kote';
tem=tem+te0; tem=tem';
for k=1:length(tem)
%asr(k)=sredRosSiegel(tem(k));
end
%asr=asr'
%t=itom();
vy=0;
end
function [ rs ] = sredRosSiegel(tem)
te0=273.15; npp=0; npp=Kramers_n(); 
dl=RasshDiapDlinVoln();
alfs=0; alfs=RasMasKoAbs(); 
%alfs=1e6*SredGraf();
%dl=dlivoln();
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4));
%Cons=(pi/2)*(C1*C2)/sig;
npk=nsreddvvsr(dl,npp,tem);
for k=1:length(npp)
c1m(k)=C1/(npp(k)^2);
c2m(k)=C2/npp(k);
Consk=(pi/2)*(c1m(k)*c2m(k))/sig/npk;
dlv(k)=dl(k)/npp(k);
%dlv(k)=dl(k);
%me=exp(C2/(dlv(k)*tem))-1;Ib(k)=2*pi*C1/(dlv(k)^5)/me;ote=(sig/Ib(k))^(1/4);
%eb=(npk^2)*sig*(tem^4);
%ote=(sig*(npk^2)/eb)^(1/4);
chi=exp(c2m(k)/dlv(k)/tem);
zna=(chi-1)^2;
Ibz(k)=Consk*chi/zna/(dlv(k)^6)/(tem^5); 
Ibc(k)=Ibz(k)/alfs(k);
end
chasc=trapz(dlv,Ibc);
chasz=trapz(dlv,Ibz);
dlsvprfo2=chasc/chasz;
dlsvprfo2=1/dlsvprfo2;
npk=npk';
rs=dlsvprfo2';
end

function rti = itom()
no=1; te0=273.15; %выбор номера ИТОМ
plotno=[440,620,860,1e3]; %кажущаяся плотность
po=[(80+82)/2,(75+78)/2,(65+68)/2,(62+65)/2]*1e-2; %пористость
sv=[60,50,25,20]*1e-2; %содержание вермикулита
rapo=[150,95,85,75,65,55,45,35,25,15,7.5,4,2,0.75,0.3,0.055,0.0055]; %в микронах
kapl=plotno(no); por=po(no); psv=sv(no); pss=1-psv;
te200=[0.09,0.12,0.183,0.23];
te380=[0.12,0.139,0.194,0.25];
te=[200,380]+te0;
kitom440=polyfit(te,[te200(1),te380(1)],1);
kitom620=polyfit(te,[te200(2),te380(2)],1);
kitom860=polyfit(te,[te200(3),te380(3)],1);
kitom1000=polyfit(te,[te200(4),te380(4)],1);
switch no
    case 1
        kti=kitom440;
    case 2
        kti=kitom620;
    case 3
        kti=kitom860;
    case 4
        kti=kitom1000;
end
srp = 10e-6;
srts = (1 - por) * srp / por;
%--------Моделирование процесса теплопереноса в ИТОМ-440
kote=ptpitom*h/(tgitom-thitom);
nom=1; %выбор Г и Х темп., ПТП для моделирования т/о в ИТОМ
no=1; te0=273.15; %выбор номера ИТОМ
h=30e-3; Thna=32e1-te0; Thde=1e-1;
temi=0; temi=1e2:1e2:5e2; n=length(temi); %моделирование процесса теплообмена в ИТОМ
qob=[2723,4235,5796]; te=[633,773,907]-te0;
ep=1e-6; nit=1e6; d=1e-2; koeq=polyfit(te,qob,2);
for k=1:n
del=1; j=1; ts=temi(k); laef(k)=polyval(kti,ts); qo(k)=polyval(koeq,ts); Thna=Thna+(k-1)*Thde;
while ((del>ep) && (j<nit))
Th(k)=Thna+(j-1)*d; Tg(k)=Th(k)+abs(qo(k)*h/laef(k)); del=(2*ts-Tg(k)-Th(k)); j=j+1;
end
tes(k)=(Th(k)+Tg(k))/2;
end
thitom=Th(nom); ptpitom=qo(nom); tgitom=Tg(nom); laef=laef';
dh=1e-4; x=0:dh:h; 
posh=21.8e-2; pov = 66.35e-2; %пористость шамота и вермикулита
hpv = srRaPorVerm(); htv = (1 - pov) * hpv / pov;
hpsh = SreRazPorSha(); htsh = (1 - posh) * hpsh / posh;
sv=[60,50,25,20]*1e-2; %содержание вермикулита
psv=sv(no); pssh=1-psv; csi=1e2; csv=round(psv*csi); cssh=round(pssh*csv); 
nsshl=round(cssh/2); nsshr=csi-csshl; tk=tgitom; titoms=0; titomss=0; x(1) = 0;
%for k=1:length(x)
for p=1:csi
    %dtdx=abs(thitom-tgitom)/h; %распределение температуры линейно
    %rastem(k)=t0-dtdx*k*dh;
    if ((p<=nsshl) || (p>=nsshr))
        ttsshk=opredToch(ttss,tk,tem);
        titoms(p) = tk - ptpitom * htsh / ttsshk;
        ktpvo = opredTeploprovVozd(tk);
        titomss(p) = tk - ptpitom * hpsh / ktpvo;
        hob = (hpsh + htsh);
        if (p > 1)
            x(p) = x(p-1) + hob;
        end
    else
        ttsvk=opredToch(ttsv,tk,tem);
        titoms(p) = tk - ptpitom * htv / ttsvk;
        ktpvo = opredTeploprovVozd(tk);
        titomss(p) = tk - ptpitom * hpv / ktpvo;
        hob = (htv + hpv);
        if (p > 1)
            x(p) = x(p-1) + hob;
        end
    end
end
%end
pl=plot(x,titoms,'-b',x+htsh,titomss,'-k');
set(pl,'LineWidth',3); hold on; grid on; 
xlabel({'Координата, м'}); 
ylabel({'Температура, К'}); 
title({'График зависимости ПТП от температуры'});
rti=0;
end