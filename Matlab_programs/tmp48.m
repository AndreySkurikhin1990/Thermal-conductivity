%определяет показатель поглощения шамота и его средний по Планку КП
function t = tmp48()
format long g;
alsh=preobMasSha_n();
dvsh=dlivoln();
n=length(dvsh);
for k=1:n
    pp(k)=alsh(k)*dvsh(k)/4/pi;
end
tem=1e2:1e1:12e1; 
tem=tem+273.15;
alsr=0;
n=length(tem);
for k=1:n
    if (rem(k,10)==0)
    %    disp(k);
    end
    alsr(k)=alphaPlankSred(tem(k),dvsh,alsh);
end
alsr=alsr'
%tem=tem
alsr=0;
for k=1:n
    if (rem(k,10)==0)
    %    disp(k);
    end
    alsr(k)=UsrePlankShaChas(alsh,dvsh,tem(k));
end
alsr=alsr'
%tem=tem'
t=postGraf(dvsh*1e6,pp);
end
function t = postGraf(dvsh,pp)
nadlvo=3;
kodlvo=10;
no=poisIndex(nadlvo,kodlvo,dvsh);
kmdv=vydPodmas(no,dvsh);
ppn=vydPodmas(no,pp);
n=length(kmdv);
d=(kmdv(n)-kmdv(1));
SrZnPP=trapz(kmdv,ppn)/d %среднее значение ПП в диапазоне от 3 до 10 мкм
p=plot(dvsh,pp,'-b');
set(p,'LineWidth',3); hold on; grid on; 
xlabel({'Длина волны, мкм'}); 
ylabel({'Показатель поглощения'}); 
title({'График зависимости ПП от ДВ'});
t=0;
end
function nsred = alphaPlankSred(tem,dv,alp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=PP*c0^2;
c2=PP*c0/PB;
ct1=c1; ct2=c2; dl=0;
npp=Kramers_n_Sha();
for k=1:length(npp)
    ct1=ct1/(npp(k)^2);
    ct2=ct2/npp(k);
    lambda=dv(k)/npp(k);
    dl(k)=lambda;
    me=(exp(ct2/(lambda*tem))-1);
    me=me*(lambda^5);
    Ib(k)=2*pi*ct1/me;
    Ibc(k)=alp(k)*Ib(k);
ct1=c1;
ct2=c2;
end
nc=trapz(dl,Ibc);
nz=trapz(dl,Ib);
nsred=nc/nz;
end
function alps = UsrePlankShaChas(alpha,dl,tem)
npp=Kramers_n_Sha();
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
nu=0;
for k=1:length(npp)
    nu(k)=c0/dl(k);
end
c1=c0;
for k=1:length(npp)
me=exp((PP*nu(k))/(PB*tem))-1;
c1=c1/npp(k);
Ib(k)=(2*PP*(nu(k)^3)/(c1^2))/me;
Ibc(k)=alpha(k)*Ib(k);
c1=c0;
end
nc=trapz(nu,Ibc);
nz=trapz(nu,Ib);
alps=nc/nz;
end
function [no] = poisIndex(nadlvo,kodlvo,dv)
n=length(dv); fn=1; fk=1; non=1; nok=1;
for k=1:n
    if ((dv(k)>kodlvo) && (fk>0))
        fk=0;
        nok=k;
    end
    if ((dv(k)>nadlvo) && (fn>0))
        fn=0;
        non=k;
    end
end
no=[non         nok];
end
function [kmdv] = vydPodmas(no,x)
non=no(1);
nok=no(2);
q=1;
for k=non:nok
    xk(q)=x(k);
    q=q+1;
end
kmdv=xk;
end