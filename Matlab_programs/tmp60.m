%КП по Росселанду для шамота
%Выбор шамота:
%1 - уплотненный (16-20 %)
%2 - низкоплотный (30-33 %)
%3 - повышенноплотный (10-16 %)
%0 - среднеплотный (24-30 %):
%0 - от 20 до 24 %
%1 - от 24 до 30 %
function qlu = tmp60()
format long g; te0=273.15;
alfs=0; alfs=1e6*SredGrafSha();
npp=0; npp=Kramers_n_Sha_uk();
dl=0; dl=dlinyvoln();
Tz=0; Tz=1e2:1e2:12e2; Tz=Tz+te0;
p=length(Tz);
vybsha=1; vystsha=0; salosha=39e-2;
for k=1:p
    npT(k)=usredVelichPlank(dl,npp,npp,Tz(k));
    alsr(k)=usredVelichPlank(dl,alfs,npp,Tz(k));
    alf(k)=sredRosSieg(Tz(k),npp,alfs,dl);
    koeftep(k)=tmp59(salosha,Tz(k),vybsha,vystsha);
end
npT=npT';
alsr=alsr';
alf=alf';
GnpT=0; 
for k=2:p
    GnpT(k-1)=(npT(k)-npT(k-1))/(Tz(k)-Tz(k-1));
end
k=p; GnpT(k)=2*GnpT(k-1)-GnpT(k-2);
sigma=5.67e-8; lamizl=0;
for k=1:p
    t=abs(8*npT(k)*sigma*(Tz(k)^3)/(3*alf(k)));
    t=t*(2*npT(k)+Tz(k)*GnpT(k));
    lamizl(k)=t;
end
t=zapvfile_alpha(alf);
t=zapvfile_lamizl(lamizl);
t=zapvfile_tem(Tz);
lamizl=lamizl'
%-----
%t=postgraf(Tz,lamizl,koeftep,p);
qlu=0;
end

function [ dl ] = dlinyvoln()
Sp=(dlvoSham1()+dlvoSham2())/2;
leSp=length(Sp);
for k=1:leSp
	Sp(k)=1e-2/Sp(k);
end
dl=Sp;
end

function t = zapvfile_alpha(dl)
p=length(dl); 
fid = fopen('Koefficient_pogloscheniya_po_Rosselandu_Shamot.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',dl(k));
end
fclose(fid);
t=0;
end

function t = zapvfile_tem(dl)
p=length(dl); 
fid = fopen('Temperatura_KP_po_Rosselandu_Shamot.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',dl(k));
end
fclose(fid);
t=0;
end

function t = zapvfile_lamizl(dl)
p=length(dl); 
fid = fopen('Koefficient_teploprovodnosti_Luchistyy_po_Rosselandu_Shamot.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',dl(k));
end
fclose(fid);
t=0;
end

function ns = usredVelichPlank(dv,uv,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
for k=1:length(npp)
    vl=c0/npp(k);
    c1=PP*(vl^2);
    c2=PP*vl/PB;
    lambda=dv(k)/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
    Ibn(k)=uv(k)*Ib(k);
dvm(k)=lambda;
end
nc=trapz(dvm,Ibn);
nz=trapz(dvm,Ib);
ns=nc/nz;
end

function t = postgraf(Tz,lamizl,koeftep,p)
pl=plot(Tz(1:p-1),lamizl,'-b',Tz,koeftep,'-k');
legend('Излучательная компонента', 'Общий КТП','location','best');
set(pl,'LineWidth',3); hold on; grid on; 
xlabel({'Температура, К'}); 
ylabel({'Составляющие коэффициента теплопроводности, Вт/(м*К)'}); 
title({'График зависимости составляющих коэффициента теплопроводности от температуры'});
t=0;
end