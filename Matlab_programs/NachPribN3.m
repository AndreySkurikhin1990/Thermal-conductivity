function NachPribN3()
format long g;
h=30e-3;
l=1e-3;
x=0:l:h;
%te0=273.15;
te0=0;
T0=1e3+te0;
heatt1=timeHeat1();
heatt2=timeHeat2();
heatt3=timeHeat3();
heatt4=timeHeat4();
temcs1=temCold1()+te0;
temcs2=temCold2()+te0;
temcs3=temCold3()+te0;
temcs4=temCold4()+te0;
temhs1=temHot1()+te0;
temhs2=temHot2()+te0;
temhs3=temHot3()+te0;
temhs4=temHot4()+te0;
ko1=vychKoefPribl(heatt1,temhs1,temcs1);
ko2=vychKoefPribl(heatt2,temhs2,temcs2);
ko3=vychKoefPrib(heatt3,temhs3,temcs3);
ko4=vychKoefPribl(heatt4,temhs4,temcs4);
tim(1)=(T0-ko1(2))/ko1(1);
tim(2)=(T0-ko2(2))/ko2(1);
tim(3)=(T0-ko3(2))/ko3(1);
tim(4)=(T0-ko4(2))/ko4(1);
tco(1)=ko1(3)*tim(1)+ko1(4);
tco(2)=ko2(3)*tim(2)+ko2(4);
tco(3)=ko3(3)*tim(3)+ko3(4);
tco(4)=ko4(3)*tim(4)+ko4(4);
%tco=tco'-te0
%tco=tco'+te0;
%tim=1*tim'
a=0;
for k=1:4
a(k)=opredTemperaturoprovVerm((tco(k)+T0)/2);
end
te1=RasTempPol(tim(1),a(1),ko1,x);
te2=RasTempPol(tim(2),a(2),ko2,x);
te3=RasTempPol(tim(3),a(3),ko3,x);
te4=RasTempPol(tim(4),a(4),ko4,x);
p=length(te1);
tsr=te2;
pl=plot(1e3*x,te2-te0,'-k');
set(pl,'LineWidth',2); 
hold on; grid on; 
xlabel({'Координата, мм'}); 
ylabel({'Температура, °С'}); 
title({'График зависимости температуры от координаты'});
%legend('te1', 'te2', 'te3', 'te4');
end

function [ te ] = RasTempPol(vre,hi,ko,x)
n=length(x);
for k=1:n
v(k)=0; 
end
a0=ko(1); b0=ko(2); v(1)=a0*vre+b0;
for k=1:n
    %hi=opredTemperaturoprovVerm(v(k));
    tm=x(k)/2/sqrt(hi*vre);
    tm2=tm^2;
    v(k)=a0*vre*((1+2*tm2)*erfc(tm)-2*tm/sqrt(pi)*exp(-tm2))+b0;
end
te=v;
end

function opr = opredTemperaturoprovVerm(T)
%te0=273.15;
te0=0;
ro=250;
te=arrTem_VVF2()+te0;
cp=arrUdTepl_VVF2();
kocp=MasKoeff(te,cp);
kolam=vychKoefPriblizhLam();
kocp=polyfit(te,cp,2);
cpt=polyval(kocp,T);
lamt=polyval(kolam,T);
opr=lamt/cpt/ro;
end


function [ koe ] = vychKoefPriblizhLam()
te0=273.15;
te=arrTem_VVF2()+te0;
lam=arrKTP_VVF2();
kolam=MasKoeff(te,lam);
kolam=polyfit(te,lam,2);
koe=kolam;
end

function [ m ] = vychKoefPribl(x, yh, yc)
koeh=polyfit(x,yh,1);
koec=polyfit(x,yc,1);
m=[koeh(1)  koeh(2) koec(1) koec(2)];
end

function [ ko ] = MasKoeff(ax,ay)
p=length(ax);
for k=1:3
    for j=1:3
        a(j,k)=0;
    end
    b(k)=0;
end
s4=0;
s3=0; sy=0;
s2=0; s2xy=0;
s1=0; s1xy=0;
for k=1:p
    s4=s4+ax(k)^4;
    s3=s3+ax(k)^3;
    s2=s2+ax(k)^2;
    s1=s1+ax(k);
    s2xy=s2xy+ay(k)*(ax(k)^2);
    s1xy=s1xy+ay(k)*ax(k);
    sy=sy+ay(k);
end
a(1,1)=s4;
a(1,2)=s3;
a(1,3)=s2;
a(2,3)=s1;
b(1)=s2xy;
b(2)=s1xy;
b(3)=sy;
a(2,1)=a(1,2);
a(3,1)=a(1,3);
a(2,2)=a(1,3);
a(3,2)=a(2,3);
a(3,3)=p;
ko=inv(a)*b';
end