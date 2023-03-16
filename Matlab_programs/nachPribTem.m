function [ temSr ] = nachPribTem(T0)
ro=250;
vtem=5;
ktp=arrKTP_VVF1();
lambda=ktp(vtem);
udte=arrUdTepl_VVF1();
cp=udte(vtem);
a=lambda/cp/ro;
h=30e-3;
l=1e-3;
te0=273.15;
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
ko1=vychKoefPrib(heatt1,temhs1,temcs1);
ko2=vychKoefPrib(heatt2,temhs2,temcs2);
ko3=vychKoefPrib(heatt3,temhs3,temcs3);
ko4=vychKoefPrib(heatt4,temhs4,temcs4);
x=0:l:h;
tim(1)=(1e3-ko1(2))/ko1(1);
tim(2)=(1e3-ko2(2))/ko2(1);
tim(3)=(1e3-ko3(2))/ko3(1);
tim(4)=(1e3-ko4(2))/ko4(1);
tco(1)=ko1(3)*tim(1)+ko1(4);
tco(2)=ko2(3)*tim(2)+ko2(4);
tco(3)=ko3(3)*tim(3)+ko3(4);
tco(4)=ko4(3)*tim(4)+ko4(4);
te1=rasTemPol(T0,a,h,ko1,x);
te2=rasTemPol(T0,a,h,ko2,x);
te3=rasTemPol(T0,a,h,ko3,x);
te4=rasTemPol(T0,a,h,ko4,x);
p=length(te1);
tsr=0;
for k=1:p
tsr(k)=(te1(k)+te2(k)+te3(k)+te4(k));
end
tsr=tsr/4;
%pl=plot(1e3*x,te1,'-k',1e3*x,te2,'-b',1e3*x,te3,'-r',1e3*x,te4,'-m');
%set(pl,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Координата, мм'}); 
%ylabel({'Температура, °С'}); 
%title({'График зависимости температуры от координаты'});
%legend('te1', 'te2', 'te3', 'te4');
%disp(te1(p));
%Nr=ceil(y0/y0l);
%Ns=ceil(y0l/l);
%q=1;
%for j=1:Nr
%    for k=1:Ns
%        te(j,k)=tsr(q); q=q+1;
%    end
%end
%for k=1:Nr-1
%    te(k,Ns+1)=te(k+1,1);
%end
%te(Nr,Ns+1)=tsr(p);
%disp(te(1,:));
temSr=tsr;
end