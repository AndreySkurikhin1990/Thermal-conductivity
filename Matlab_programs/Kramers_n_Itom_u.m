function [ PokPrel ] = Kramers_n_Itom_u()
format long g;
Sp=(dlvoItom1()+dlvoItom2())/2;
leSp=length(Sp);
c0=299792458;
enu=1e-3;
alsr=SredGrafItom();
for k=1:leSp
    Sp(k)=1e-2/Sp(k);
    om(k)=2*pi*c0/Sp(k);
end
oms=izmMasChast(om,enu);
nnus=PokazPrelomAl(om,oms,alsr,c0);
w1=6e-1; w2=1e0-w1; 
rover=25e1; rosha=219e1;
v1=w1/rover; v2=w2/rosha;
phi1=v1/(v1+v2); phi2=v2/(v1+v2); 
dob=(1.56-1e0)*phi2+(1.543-1e0)*phi1; %dob=(1.543-1)*0.6+0.4*(1.56-1);
nnu=PokazPrelomPribl(om,oms,nnus)+dob;
%n0=postrGraf(nnu,nnus,om,oms);
PokPrel=nnu;
end

function pg = postrGraf(nom,noms,om,oms)
dlnu=izmMasDlVol(om);
dlnus=izmMasDlVol(oms); 
p=length(dlnu); 
dvko=20e-6;
dlnuo=0; dlnuso=0;
j=1; nnuo=0; 
for k=1:p  
    if (dlnu(k)<dvko)    
      nnuo(j)=nom(k); 
      dlnuo(j)=dlnu(k); 
      j=j+1; 
    end
end
p1=j-1;
p=length(dlnus);
j=1; na=0;
for k=1:p  
    if (dlnus(k)<20e-6)    
        na(j)=noms(k);
        dlnuso(j)=dlnus(k); 
        j=j+1;
    end
end
p2=j-1;
pl=plot(dlnuo(2:p1-1)*1e6,nnuo(2:p1-1),'-b',dlnuso(2:p2-1)*1e6,na(2:p2-1),'-k');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'}); 
ylabel({'Показатель преломления'}); 
title({'График зависимости показателя преломления от длины волны'});
pg=0;
end

function [ dl ] = izmMasDlVol(om)
la=0;
c0=299792458;
for k=1:length(om)
    la(k)=2*pi*c0/om(k);
end
dl=la;
end