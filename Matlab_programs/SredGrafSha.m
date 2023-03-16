function [ alsr ] = SredGrafSha()
dlvo=dlvoSham1(); 
dlvo2=dlvoSham2(); 
Tal=TrkoSham1(); 
Tal2=TrkoSham2(); 
ndlvo=length(dlvo);
mkbr=238.57; 
msh=1.706; 
tol=0.7; 
rokbr=2.75; 
rsh=2.7; 
mkbr=mkbr-msh;
mkbr2=227.973; 
msh2=1.1; 
tol2=0.68; 
rsh2=rsh; 
mkbr2=mkbr2-msh2;
vkbr=mkbr/(1e3*rokbr);
vsh=msh/(1e3*rsh);
vkbr2=mkbr2/(1e3*rokbr);
vsh2=msh2/(1e3*rsh2);
xsh=(vsh/(vsh+vkbr))*tol*1e3;
xsh2=(vsh2/(vsh2+vkbr2))*tol2*1e3;
for k=1:ndlvo      
    dlvo(k)=(1e4/dlvo(k));         
    Tal(k)=-log(Tal(k))/xsh; 
end
ndlvo2=length(dlvo2);
for k=1:ndlvo2 
    dlvo2(k)=(1e4/dlvo2(k));     
    Tal2(k)=-log(Tal2(k))/xsh2;    
end
%alsre=SreZnaPog();
%alsre=postrgr(dlvo,dlvo2,Tal,Tal2);
for k=1:ndlvo
    alsre(k)=(Tal(k)+Tal2(k))/2;
end
%als=zapvfile(1e6*alsre);
alsr=alsre;
end

function pg = postrgr(x1,x2,y1,y2)
p=plot(x1,y1,'-k',x2,y2,'-b');
set(p,'LineWidth',2); 
hold on; grid on; 
xlabel({'Длина волны, мкм'}); 
ylabel({'Коэффициент поглощения, мкм(-1)'}); 
title({'График зависимости КП от длины волны'});
legend('alpha1', 'alpha2','location','best');
pg=0;
end

function szkp = SreZnaPog()
te0=273.15;
te = [350 400 450 500 550 800 1000 1500 2000];
kp = [0.8625 0.8438 0.825 0.8 0.7786 0.6545 0.5409 0.3425 0.2841];
for k =1:length(te)
al(k)=SreSteChe(te(k)); 
end
al=1*al'
p=plot(te,kp,'-k',te,al,'-b');
set(p,'LineWidth',2); 
hold on; grid on; 
xlabel({'Температура, К'}); 
ylabel({'Степень черноты (поглощательная способность)'}); 
title({'График зависимости epsilon от температуры'});
legend('alpha1', 'alpha2','location','best');
szkp=0;
end

function kp = SreSteChe(T)
npp=Kramers_n();
dv=dlvoSham1(); 
for k=1:length(dv)
    dv(k)=1e-2/dv(k);
end
eps=epsillam(npp);
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
iz1(j)=eps(j)*iz0(j);
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
kp=chi/zna;
end

function [ epsla ] = epsillam(pp)
ep=0;
for k=1:length(pp)
    n=pp(k);
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(n)/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log((n-1)/(n+1))*((n^2-1)^2)/((n^2+1)^3);
ep(k)=eps;
end
epsla=ep;
end

function t = zapvfile(dl)
format long g; 
p=length(dl); 
fid = fopen('Koefficient_pogloscheniya_sha.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',dl(k));
end
fclose(fid);
t=0;
end