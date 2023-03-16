function [ temSr ] = nachPribTemN2()
%Netzsch data obtained used
format long g;
te0=273.15;h=30e-3; l=1e-3; x=0:l:h;
tsr2=opredTsredMaxTemp(1e3,x);
tsr3=opredTsredMaxTemp(8e2,x);
tsr1=opredTsredMaxTemp((585*2+6e2)/3,x);
T1 = [(585*2+6e2)/3   (377+(383*2+384*2)/4+396)/3 ((120*3+119)/4+129+138)/3];
T1=T1+te0;
T2 = [1e3   ((697*3+698)/4+703)/2 ((261*3+262)/4+273)/2];
T2=T2+te0;
T3 = [8e2   548 2e2];
T3=T3+te0;
teex1=NachPribN1(T1(1),T1(3),T1(2),x);
teex2=NachPribN1(T2(1),T2(3),T2(2),x);
teex3=NachPribN1(T3(1),T3(3),T3(2),x);
pl=plot(1e3*x,tsr1,'--k',1e3*x,teex1,'-k',1e3*x,tsr2,'--b',1e3*x,teex2,'-b',1e3*x,tsr3,'--r',1e3*x,teex3,'-r');
set(pl,'LineWidth',2); 
hold on; grid on; 
xlabel({'Координата, мм'}); 
ylabel({'Температура, K'}); 
title({'График зависимости температуры от координаты'});
%legend('tsr', 'teex');
temSr=tsr;
end

function ko = Korni(lc,pc,n)
lc(n+1)=lc(n+1)-pc;
ro=roots(lc);
kor=0; q=1; p=length(ro);
for k=1:p
    if (ro(k)>0)
        if (imag(ro(k))==0)
        kor(q)=ro(k);
        q=q+1;
        end
    end
end
min=1e30; p=length(kor);
for k=1:p
    if (kor(k)<min)
        min=kor(k);
    end
end
ko = min;
end

function [ tem ] = RascRaspTemp(T0,a,h,mkh,mkc,x)
a0=mkh(2); %q
a2=mkh(1)/2/a; %n/(2a)
a3=(mkc(1)-2*a*a2)/h/6/a; %m/(2a)
a1=(mkc(2)-a0-a3*(h^3)-a2*(h^2))/h;
ti=(T0-a0)/(2*a*a2);
lx=length(x);
teta=0;
for k=1:lx
    teta(k)=(6*a*a3*x(k)+2*a*a2)*ti+(a3*x(k)^3+a2*x(k)^2)+a1*x(k)+a0;
end
for k=1:lx
    teta(k)=(6*a*a3*x(k)+2*a*a2)*ti+(a3*x(k)^3+...
        a2*x(k)^2)+a1*x(k)+a0+a4*(x(k)^4+12*a*(x(k)^2)*ti+...
        12*(a*ti)^2)+a5*(x(k)^5+20*a*(x(k)^3)*ti+...
        60*((a*ti)^2)*x(k))+a6*(x(k)^6+30*a*(x(k)^4)*ti+...
        180*((a*x(k)*ti)^2)+120*((a*ti)^3))+a7*(x(k)^7+...
        42*(x(k)^5)*a*ti+420*x(k)*(a*ti*x(k))^2+...
        840*(x(k)^2)*((a*ti)^3)+420*((ti*a)^4))+...
        a8*(x(k)^8+63*a*ti*(x(k)^6)+945*((a*(x(k)^2)*ti)^2)+...
        3780*((a*x(k)*ti)^2)*a*ti+5670*((a*ti)^4));
end
tem = teta;
end

function ov = opreTrebVeli(x,y,a)
n = length(x);
f = 1; p=0;
for k=1:n
    if (x(k)>a)
    if (f > 0)
        p = k;
        f = 0;
    end
    end
end
yn = 0;
if (p > 1)
yn = y(p-1) + (y(p) - y(p-1)) * (a - x(p-1)) / (x(p) - x(p-1));
end
ov = yn;
end

function [ ts ] = opredTsredMaxTemp(Tko,x)
te0=273.15;
Tko=Tko+te0;
ro=250;
h=x(length(x));
tem=arrTem_VVF2();
ktp=arrKTP_VVF2();
udte=arrUdTepl_VVF2();
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
n=1;
ko1h=polyfit(heatt1,temhs1,n);
ko1c=polyfit(heatt1,temcs1,n);
ko2h=polyfit(heatt2,temhs2,n);
ko2c=polyfit(heatt2,temcs2,n);
ko3h=polyfit(heatt3,temhs3,n);
ko3c=polyfit(heatt3,temcs3,n);
ko4h=polyfit(heatt4,temhs4,n);
ko4c=polyfit(heatt4,temcs4,n);
tim(1)=Korni(ko1h,Tko,n);
tim(2)=Korni(ko2h,Tko,n);
tim(3)=Korni(ko3h,Tko,n);
tim(4)=Korni(ko4h,Tko,n);
tim=1*tim'
tco(1)=polyval(ko1c,tim(1));
tco(2)=polyval(ko2c,tim(2));
tco(3)=polyval(ko3c,tim(3));
tco(4)=polyval(ko4c,tim(4));
tco=1*tco'
Ts = 0;
for k=1:4
    Ts(k)=(tco(k)+Tko)/2;
    lambda(k)=opreTrebVeli(tem,ktp,Ts(k));
    cp(k)=opreTrebVeli(tem,ktp,Ts(k));
    atp(k)=lambda(k)/ro/cp(k);
end
te1=RascRaspTemp(Tko,atp(1),h,ko1h,ko1c,x);
te2=RascRaspTemp(Tko,atp(2),h,ko2h,ko2c,x);
te3=RascRaspTemp(Tko,atp(3),h,ko3h,ko3c,x);
te4=RascRaspTemp(Tko,atp(4),h,ko4h,ko4c,x);
p=length(te1);
tsr=0;
for k=1:p
tsr(k)=(te1(k)+te2(k)+te3(k)+te4(k));
end
tsr=tsr/4;
ts=tsr;
end