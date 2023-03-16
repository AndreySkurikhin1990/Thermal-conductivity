function [ PokPrel ] = Kramers_n_uk()
format long g;
c0=299792458;
enu=1e-3;
Sp=dlinyvoln();
alsr=1e6*SredGraf();
%alsr=1e6*SredGrafKoAbsVer();
%alsr=skleiv2Mas(soglMas(1e6*SredGraf(),1e6*SredGraSt()),1e6*SredGraf());
%alsr=RasMasKoAbs(); 
%n0=ZapisFileOptio(alsr);
%Sp=RasshDiapDlinVoln();
leSp=min([length(alsr),length(Sp)]);
%n0=postrGraf(alsr,alph,Sp,dlinvolny());
%Sp=dlivoln();
%alsr=1e6*SredGraf();
for k=1:leSp
    hi(k)=Sp(k)*alsr(k)/4/pi;
    nu(k)=2*pi*c0/Sp(k);
end
nus=izmMasChast(nu,enu);
nnus=PokazPrelomAl(nu,nus,alsr,c0)-1+1.543;
nnu=PokazPrelomPribl(nu,nus,nnus);
%n0=postrGraf(nnu,nnus,nu,nus);
%hia=PokazPoglPribl(nus,alsr,nu);
%ep=opredEpsilSt(nu/2/pi,nus/2/pi,hi,hia)
%ns=PokazPrelomProv2(nu,nus,hi,hia,nnu,nnus);
%------
%n0=opredn0(Sp,alsr);
%n0=ZapisFile(nnu);
%nnu=vyborka(nnu,Sp)-c0;
PokPrel=nnu;
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

function t = ZapisFileOptio(massi)
fid = fopen('Koefficient_pogloscheniya_ver.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFile(massi)
fid = fopen('Pokazatel_prelomleniya.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function [ oprn ] = opredpokprel(nu,nus,hi,epsSt,leSp)
n=0;
for k=1:leSp
f=0;
for j=1:leSp
    if (nu(j)==nus(k))
        f(j)=0;
    else
        f(j)=(hi(j)*nu(j))/(nu(j)^2-nus(k)^2);
    end
end
n(k)=sqrt(epsSt)+(2/pi)*integpo2ma(nu,f);
end
oprn=n;
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

function [ dl ] = izmMasDlVol(lam)
om=0;
c0=299792458;
lenu=length(lam);
for k=1:lenu
    om(k)=2*pi*c0/lam(k);
end
dl=om;
end

function [ dlv ] = dlinyvoln()
format long g;
dl=0; dl=dlvoVer53101(); 
p=length(dl); 
fid = fopen('Dliny_voln_ver.txt','w');
for k=1:p  
    fprintf(fid,'%0.20f\n',dl(k));
    dl(k)=1e-2/dl(k);
end
fclose(fid);
dlv = dl;
end