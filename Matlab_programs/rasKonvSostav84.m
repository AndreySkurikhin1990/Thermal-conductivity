function koTeplop = rasKonvSostav84()
%Фракция 8-4 мм
y0=30e-3;
l=1e-3;
x=0:l:y0;
te0=273.15;
por=0.5575;
g=9.81;
ro=20e1;
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
Th=(temvh(16)+temvh(18))/2;
Tc=(temvc(16)+temvc(18))/2;
tep=(tepv(16)+tepv(18))/2;
Tm=nachPribTem84(Th,tep);
Tm=pribltemkTc(Tm,Th,Tc,x);
Nt=length(Tm);
ktpm=danPoTemTepl(temvs,tepv);
psi1=0;
psi2=0;
ktpr=0;
vtem=5;
udte=arrUdTepl_VVF1();
%cp=udte(vtem);
for k=1:Nt
    Tsre=(Tm(k)+Th)/2;
    beta0=1/(Tsre);
    ktpr=ktpm(1)*Tsre+ktpm(2)
PrVo=vychKoefVoz(1,Tsre);
roVoz=vychKoefVoz(2,Tsre);
muVoz=vychKoefVoz(3,Tsre);
lambVoz=vychKoefVoz(4,Tsre);
cp=vychKoefVoz(5,Tsre);
Gras=beta0*(Th-Tm(k))*g*(deq^3)*(roVoz^2)/(muVoz^2)
psi1(k)=1+2*beta0*(-Th+Tm(k))*g*(roVoz^2)*(por^3)*cp*x(k)/18/(apzeos^2)/muVoz/ktpr;
psi2(k)=1-(16e-5)*PrVo*Gras*((por*x(k))^2)*lambVoz/(deq^2)/ktpr;
end
for k=1:Nt
    if (rem(k-1,5)==0)
    te=Tm(k)
    ktpr=ktpm(1)*Tsre+ktpm(2);
    xk=x(k)
    psi11=psi1(k);
    psi12=psi1(k);
    kot=(psi11-1)*ktpr
    kott=(psi12-1)*ktpr
    psi11=1-psi1(k)
    psi12=1-psi1(k)
    end
end
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
c=1*c
%pl=plot(1e3*x,koTep,'-b',1e3*x,koTe,'-k');
%pl=plot(1e3*x,koTep,'-b');
%set(pl,'LineWidth',2); 
%hold on; 
%grid on; 
%xlabel({'Координата, мм'}); 
%ylabel({'Коэффициент теплопроводности без конвекции, Вт/(м*К)'}); 
%title({'График зависимости температуры от координаты'});
koTeplop=0;
end