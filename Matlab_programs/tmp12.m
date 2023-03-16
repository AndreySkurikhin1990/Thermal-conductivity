%расчет КП по Росселанду
format long g; teta=0; 
alfs=0; alfs=RasMasKoAbs(); 
npp=0; npp=Kramers_n(); 
te0=273.15;
Tz=1e2:1e-1:11e2;
Tz=Tz+te0;
p=length(Tz); alf=0;
for k=1:p
    %alf(k)=1e6/sredRosSieg(Tz(k));
    alf(k)=knusreddvvsrVer(Tz(k));
    disp(k);
end
Tz=1*Tz'
alf=1*alf'
pl=plot(Tz,alf,'-b');
set(pl,'LineWidth',3); hold on; grid on; 
xlabel({'Температура, К'}); 
ylabel({'Коэффициент поглощения по Росселанду, 1/м'}); 
title({'График зависимости КП по Росселанду от температуры'});