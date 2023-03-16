%определяет показатель преломления шамота
function [ PokPrel ] = Kramers_n_Sha_uk()
format long g;
c0=299792458;
Sp=dlinyvoln();
leSp=length(Sp);
alsr=1e6*SredGrafSha();
enu=ZapisFileOptio(alsr);
for k=1:leSp
    nu(k)=2*pi*c0/Sp(k);
end
enu=1e-3;
ma=izmMasChast(nu,enu); %определение массива nus на enu от шага nu
na=PokazPrelomAl(nu,ma,alsr,c0); %главная часть - определение ПП по КК
nnu=PokazPrelomPribl(nu,ma,na); %определение второго массива ПП
nnu=nnu-1+1.56;
%enu=ZapisFile(nnu);
PokPrel=nnu;
end

function t = ZapisFile(massi)
fid = fopen('Pokazatel_prelomleniya_shamota.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptio(massi)
fid = fopen('Koefficient_pogloscheniya_shamota.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
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