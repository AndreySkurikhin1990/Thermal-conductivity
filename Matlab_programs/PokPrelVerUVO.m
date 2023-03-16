function np = PokPrelVerUVO()
veruv=vermFlake50mkm()/50e-6;
veruvdl=vermFlake50mkmDlVoln();
verop=(vermFlake100mkm()+1)/100e-6;
veropdl=vermFlake100mkmDlVoln();
vermuv=vermuvdiap(veruv,veruvdl);
vermop=vermoptdiap(verop,veropdl);
np=postrGrafik(vermuv,vermop);
end

function [ raz ] = razbien(na,ko,nadl,kodl,hdl)
koe=(ko-na)/(kodl-nadl);
sc=ko-koe*kodl;
dli=nadl:hdl:kodl;
ar=0;
for k=1:length(dli)
    ar(k)=koe*dli(k)+sc;
end
raz=ar;
end

function [ uv ] = vermuvdiap(veruv,veruvdl)
p=length(veruv);
hdl=1e-1;
absrb=0;
q=1;
for k=1:p-1
    na=veruv(k);
    ko=veruv(k+1);
    nadl=veruvdl(k);
    kodl=veruvdl(k+1);
    tr=0;
    tr=razbien(na,ko,nadl,kodl,hdl);
    for j=1:length(tr)
        absrb(q)=tr(j);
        if (q==1)
            dlv(q)=veruvdl(1);
        else dlv(q)=dlv(q-1)+hdl;
        end
        q=q+1;
    end
end
tmp=0;
for k=1:length(dlv)
    tmp(1,k)=dlv(k);
    tmp(2,k)=absrb(k);
end
uv=tmp;
end

function [ opt ] = vermoptdiap(verop,veropdl)
p=length(verop);
hdl=1e-1;
absrb=0;
q=1;
for k=1:p-1
    na=verop(k);
    ko=verop(k+1);
    nadl=veropdl(k);
    kodl=veropdl(k+1);
    tr=0;
    tr=razbien(na,ko,nadl,kodl,hdl);
    for j=1:length(tr)
        absrb(q)=tr(j);
        if (q==1)
            dlv(q)=veropdl(1);
        else dlv(q)=dlv(q-1)+hdl;
        end
        q=q+1;
    end
end
tmp=0;
for k=1:length(dlv)
    tmp(1,k)=dlv(k);
    tmp(2,k)=absrb(k);
end
opt=tmp;
end

function n = postrGrafik(vuv,vop)
dluv=vuv(1,:);
dlop=vop(1,:);
abuv=vuv(2,:);
abop=vop(2,:);
q=1;
dlv=0;
absr=0;
for k=1:(length(dluv)+length(dlop))
    if (k<(length(dluv)+1))
    dlv(k)=dluv(k);
    absr(k)=abuv(k);
    else
        dlv(k)=dlop(q);
        absr(k)=abop(q);
        q=q+1;
    end
end
c0=299792458;
ome=0; hi=0; ledl=length(dlv);
for k=1:ledl
    dlv(k)=dlv(k)*1e-9;
    hi(k)=dlv(k)*absr(k)/4/pi;
    ome(k)=2*pi*c0/dlv(k);
end
omes=izmMasChast(ome);
np=PokazPrelomAlfa(ome,omes,absr,c0);
ns=0; dliv=0; q=1;
for k=1:length(np)
    if (dlv(k)>4e-7)
        dliv(q)=dlv(k);
        ns(q)=np(k)-1+1.563;
        q=q+1;
    end
end
pl=plot(dliv*1e9,ns,'-b');
set(pl,'LineWidth',2);
hold on; grid on; 
xlabel({'Длина волны, нм'}); 
ylabel({'Показатель преломления'}); 
title({'График зависимости показателя преломления от длины волны'});
n=0;
end

function [ pop ] = PokazPrelomAlfa(nu,nus,ka,vl)
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
    np(k)=1+(vl/pi)*integpo2mas(nu,fn)/2/nu(k);
end
pop=np;
end

function inte = integpo2mas(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ nus ] = izmMasChast(nu)
eps=1e-3;
nust=0;
lenu=length(nu);
for k=1:lenu-1
    nust(k)=eps*(nu(k+1)-nu(k))+nu(k);
end
nust(lenu)=2*nust(lenu-1)-nust(lenu-2);
nus=nust;
end