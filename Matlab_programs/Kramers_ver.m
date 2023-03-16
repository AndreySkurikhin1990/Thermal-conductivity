%для моделирования рассеяния излучения
function [ PokPrel ] = Kramers_ver()
format long g;
c0=299792458;
enu=1e-3;
alsr=RasMasKoAbs(); 
Sp=RasshDiapDlinVoln();
leSp=length(Sp);
nu=0; nus=0; hi=0;
for k=1:leSp
    hi(k)=Sp(k)*alsr(k)/4/pi;
    nu(k)=2*pi*c0/Sp(k);
end
nus=izmMasChast(nu,enu); %массив циклических частот a
nnu=PokazPrelomAl(nu,nus,alsr,c0);
%n0=postrGraf(nnu,nnus,nu,c0);
%n0=opredn0(Sp,alsr);
%ns0=opredn0(Sp,als);
nnu=vyborka(nnu,Sp);
%n0=Sravnen(nnu);
PokPrel=nnu';
end

function t = Sravnen(nnu)
nnuk=Kramers_n_uk();
Sp=0; Sp=dliny_voln();
c0=0.543;
nnuk=vyborka(nnuk,Sp)-c0;
nnukr=abs(nnu-nnuk);
mdnu=max(nnukr)
for k=1:length(nnukr)
    if (mdnu==nnukr(k))
        q=k;
        break;
    end
end
q=q'
disp(nnuk(q));
disp(nnu(q));
t=0;
end

function [ dlv ] = dliny_voln()
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k); 
end
dlv = dl;
end

function [ np ] = vyborka(nnu,Sp)
n=length(Sp); shag=1e-1;
nac=1.4; kon=27;
mr=nac:shag:kon;
mr=mr*1e-6;
for k=1:length(mr)
    la=mr(k); f=1; q=2;
    for j=2:n
        if ((Sp(j)>=la) && (f>0))
            q=j; 
            f=0; 
            break;
        end
    end
    if (f<1)
        pr=(nnu(q)-nnu(q-1))/(Sp(q)-Sp(q-1));
        nonu(k)=pr*(la-Sp(q-1))+nnu(q-1);
    end
end
np=nonu;
end

function inte = integpo2ma(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ nus ] = izmMasChast(nu,enu)
nust=0;
lenu=length(nu);
for k=1:lenu-1
    nust(k)=enu*(nu(k+1)-nu(k))+nu(k);
end
nust(lenu)=enu*(nust(lenu-1)-nust(lenu-2))+nust(lenu-1);
nus=nust;
end

function [ pop ] = PokazPrelomAl(nu,nus,ka,vl)
np=0;
p=length(nu);
for k=1:p
    fn=0; q=1;
    for j=1:p-1
        dkpodo=(ka(j+1)-ka(j))/(nu(j+1)-nu(j)); %угловые частоты: nus - a, nu - omega
        podln=((nu(j)+nus(k))/(nu(j)-nus(k)));
        podln=abs(podln);
        fn(q)=dkpodo*log(podln);
        q=q+1;
    end
    fn(p)=2*fn(p-1)-fn(p-2);
    np(k)=1+(vl/pi)*integpo2ma(nu,fn)/2/nu(k);
end
nnu=PokazPrelomPribl(nu,nus,np);
pop=np;
end

function n0 = opredn0(dl,ka)
n0=1+(1/2/pi^2)*integpo2ma(dl,ka);
end

function [ ns ]  = PokazPrelomPribl(nu,nus,na)
nnu=0;
p=length(nu);
for k=2:p
    f=1; q=1;
    for j=1:p
        if ((nus(j)>nu(k)) && (f>0))
            f=0; 
            q=j;
            break;
        end
    end
    de=(na(k)-na(k-1))/(nu(k)-nu(k-1));
    nnu(k)=de*(nus(q)-nu(k))+na(k-1);
end
nnu(p)=(-nus(p)+nu(p))*(na(p)-na(p-1))/(nu(p)-nu(p-1))+na(p);
ns = nnu;
end

function pg = postrGraf(nnu,nnus,nu,vl)
p=length(nu); dvko=27; dvko=dvko*1e-6;
for k=1:p
    nu(k)=nu(k)/2/pi;
    dl(k)=vl/nu(k);
end
dlnuo=0; nnuo=0; j=1; 
for k=1:p  
    if (dl(k)<dvko)
      nnuos(j)=nnus(k);
      nnuo(j)=nnu(k); 
      dlnuo(j)=dl(k); 
      j=j+1; 
    end
end
p=j-1;
pl=plot(dlnuo(2:p-1)*1e6,nnuo(2:p-1),'-b',dlnuo(2:p-1)*1e6,nnuos(2:p-1),'-k');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'}); 
ylabel({'Показатель преломления'}); 
title({'График зависимости показателя преломления от длины волны'});
legend('Пропускание','Абсорбция','location','best')
pg=0;
end