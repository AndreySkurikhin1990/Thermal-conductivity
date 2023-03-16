%Только для шамота: определение средних значений степени черноты, показателя преломления и его
%квадрата, отражательной способности, коэффициента поглощения
function t = tmp53()
te0=273.15; te=2e1:1e0:2e3; te=te+te0;
npp=Kramers_n_Sha();
dv=(dlvoSham1()+dlvoSham2())/2;
t=ZapisFileOptiodv(dv);
%dv=dlvoVer53101(); 
for k=1:length(dv)
    dv(k)=1e-2/dv(k);
end
%posp=SreZnaPogSpos(te,npp,dv); %поглощательная способность
kp=SreZnaKoefPogl(te,npp,dv); %коэффициент поглощения
%n2s=SreZnaPokPrel2(te,npp,dv); %показатель преломления в квадрате
%ns=SreZnaPokPrel(te,npp,dv); %показатель преломления
t=0;
end
%определение среднего значения степени черноты
function szkp = SreZnaPogSpos(te,npp,dv)
eps=epsillam(npp);
for k =1:length(te)
ab(k)=usrednen(te(k),eps,npp,dv); 
end
ab=abs(ab');
for k=1:length(ab)
    refl(k)=1-ab(k);
end
refl=refl'
%p=plot(te,ab,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Степень черноты (поглощательная способность)'}); 
%title({'График зависимости степени черноты от температуры'});
szkp=0;
end
%определение среднего значения коэффициента поглощения
function szkp = SreZnaKoefPogl(te,npp,dv)
als=1e6*SredGrafSha();
for k =1:length(te)
kp(k)=usrednen(te(k),als,npp,dv); 
end
p=ZapisFileOptio(kp);
kp=abs(kp')
p=plot(te,kp,'-b');
set(p,'LineWidth',2); 
hold on; grid on; 
xlabel({'Температура, К'}); 
ylabel({'Коэффициент поглощения'}); 
title({'График зависимости коэффициента поглощения от температуры'});
szkp=0;
end

function t = ZapisFileOptio(massi)
fid = fopen('Koefficient_pogloscheniya_sha_T.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.12f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiodv(massi)
fid = fopen('Dlina_volny_sha.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

%определение среднего значения квадрата показателя преломления
function szkp = SreZnaPokPrel2(te,npp,dv)
for k=1:length(npp)
    npp2(k)=npp(k)*npp(k);
end
for k =1:length(te)
npp2s(k)=usrednen(te(k),npp2,npp,dv); 
end
npp2s=real(npp2s')
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости квадрата показателя преломления от температуры'});
szkp=0;
end
%определение среднего значения квадрата показателя преломления
function szkp = SreZnaPokPrel(te,npp,dv)
for k =1:length(te)
npps(k)=usrednen(te(k),npp,npp,dv); 
end
npps=real(npps')
%p=plot(te,kp,'-b');
%set(p,'LineWidth',2); 
%hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициент поглощения'}); 
%title({'График зависимости показателя преломления от температуры'});
szkp=0;
end

function uv = usrednen(T,usv,npp,dv)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*c0^2;
c2=PP*c0/PB;
iz0=0; iz1=0; dl=0;
c1t=c1; c2t=c2;
for j=1:length(dv)
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
iz0(j)=c1t/((la^5)*(exp(c2t/(la*T))-1));
iz1(j)=usv(j)*iz0(j); %усредняемая величина
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
uv=chi/zna; %усредненная величина
end

%определение поглощательной способности по формуле Данкла
function [ epsla ] = epsillam(pp)
ep=0;
for k=1:length(pp)
    n=abs(pp(k));
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(n)/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log((n-1)/(n+1))*((n^2-1)^2)/((n^2+1)^3);
ep(k)=abs(eps);
eps=0;
end
epsla=real(ep);
end