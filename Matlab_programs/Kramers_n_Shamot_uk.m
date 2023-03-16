%определяет показатель преломления шамота
function [ PokPrel ] = Kramers_n_Shamot_uk()
format long g;
c0=299792458;
Sp=dlinyvoln();
leSp=length(Sp);
alsr=1e6*SredGrafSha();
for k=1:leSp
    nu(k)=2*pi*c0/Sp(k);
end
enu=1e-3;
ma=izmMasChast(nu,enu); %определение массива nus на enu от шага nu
na=PokazPrelomAl(nu,ma,alsr,c0); %главная часть - определение ПП по КК
nnu=PokazPrelomPribl(nu,ma,na); %определение второго массива ПП
%als=usredVelich(Sp,alsr,T,nnu)
nnu=vyborka(nnu,Sp);
PokPrel=nnu;
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
        dkpodo=(ka(j+1)-ka(j))/(nu(j+1)-nu(j));
        podln=((nu(j)+nus(k))/(nu(j)-nus(k)));
        podln=abs(podln);
        fn(q)=dkpodo*log(podln);
        q=q+1;
    end
    fn(p)=fn(p-1);
    np(k)=1+(vl/pi)*integpo2ma(nu,fn)/2/nu(k);
end
pop=np;
end

function [ ns ]  = PokazPrelomPribl(nu,nus,na)
nnu=0;
p=length(nu);
for k=1:p-1
    de=(-nus(k)+nu(k))*(na(k+1)-na(k))/(nus(k+1)-nus(k));
    nnu(k)=de+na(k);
end
nnu(p)=(-nus(p)+nu(p))*(na(p)-na(p-1))/(nu(p)-nu(p-1))+na(p);
ns = nnu;
end

function [ dl ] = dlinyvoln()
Sp=(dlvoSham1()+dlvoSham2())/2;
leSp=length(Sp);
for k=1:leSp
	Sp(k)=1e-2/Sp(k);
end
dl=Sp;
end

function [ np ] = vyborka(nnu,Sp)
n=length(Sp); shag=1e-1;
nac=1.4; kon=27;
mr=nac:shag:kon;
mr=mr*1e-6;
for k=1:length(mr)
    la=mr(k); f=1; q=2;
    for j=2:n
        if ((Sp(j)>la) && (f>0))
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

function es = usredVelich(dv,usrvel,T,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*(c0^2);
c2=PP*c0/PB;
iz0=0; iz1=0; dl=0;
c1t=c1; c2t=c2;
for j=1:length(dv)
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
    t=(exp(c2t/(la*T))-1);
iz0(j)=c1t/((la^5)*t); iz0(j)=ProvAdekv(iz0(j));
iz1(j)=usrvel(j)*iz0(j);
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
es=ProvAdekv(chi)/ProvAdekv(zna);
end