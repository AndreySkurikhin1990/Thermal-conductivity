function knus = knusreddvvsrVer(tem,vy)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
npp=Kramers_n();
dv=RasshDiapDlinVoln();
knu=RasMasKoAbs();
s=1;
if (vy==1)
    s=opredKO(tem);
end
knu=s*knu;
c1=PP*c0^2;
c2=PP*c0/PB;
ct1=c1; ct2=c2;
for k=1:length(npp)
    ct1=ct1/(npp(k)^2);
    ct2=ct2/npp(k);
    lambda=dv(k)/npp(k);
    dv(k)=lambda;
    me=(exp(ct2/(lambda*tem))-1);
    Ib(k)=2*pi*ct1/(lambda^5)/me;
    Ibc(k)=knu(k)*Ib(k);
ct1=c1;
ct2=c2;
end
nc=trapz(dv,Ibc);
nz=trapz(dv,Ib);
knus=nc/nz;
end

function w = opredKO(tem)
dko2m=[4.68,4.69,5.65,13.17,20.2,27.81]/1e2; 
Tna=6e2; dTe=2e2; Tko=16e2;
dko2t=Tna:dTe:Tko;
if (tem<Tna)
    dko2=1;
else
for k=1:length(dko2m)
    dko2m(k)=1-dko2m(k);
end
dko2=opredKTPTKTochVer(dko2m, dko2t, tem);
end
wmg=217*1e-3; wal=132*1e-3; wsi=368*1e-3;
koal=rasKoefAlpha(tem,wmg,wsi,wal); 
s=mean(koal);
dko=0.444905; 
dko3=1.3002914; 
w=s*dko*dko2*dko3;
end

function t = opredKTPTKTochVer(ktptks, te, temp)
n=length(te); f=1; p=0;
for k=1:n
if ((te(k)>=temp) && (f>0)) 
        p=k; f=0;
end
end
    if ((p==1) && (f==0)) 
        p=2; 
    end
    if ((f==1) && (p==1)) 
        p=n; f=0;
    end
    if (f==0)
	ko=(ktptks(p)-ktptks(p-1))/(te(p)-te(p-1)); 
    ktp=ktptks(p-1)+ko*(temp-te(p-1));
    else
        ktp=0; 
        ko=0;
    end
    if ((ktp<0) || (ktp>1))
        ktp=1;
    end
t=ktp;
end

function [ ko ] = rasKoefAlpha(te,wmg,wsi,wal)
Tna=3e2; Tko=16e2; delT=1e2; tem=Tna:delT:Tko;
mgo = [1 0.94 0.87 0.801 0.736 0.676 0.635 0.59 ...
    0.567 0.543 0.53 0.525 0.515 0.507];
al2o3 = [1 0.98 0.946 0.898 0.854 0.797 0.753 ...
    0.709 0.676 0.642 0.618 0.598 0.582 0.57];
sio2 = [1 0.984 0.953 0.917 0.854 0.808 0.756 ...
    0.711 0.578 0.523 0.495 0.468 0.448 0.429];
wo=(wmg+wsi+wal);
kumm=vychKoe(tem,[te],mgo);
kuma=vychKoe(tem,[te],al2o3);
kums=vychKoe(tem,[te],sio2);
n=length([te]);
for k=1:n
    koo(k)=kumm(k)*wmg+kuma(k)*wal+kums(k)*wsi;
end
ko=koo/wo;
end

function [ koe ] = vychKoe(tem,te,epra)
n=length(tem);
m=length(te);
for q=1:m
    f=1; j=1;
for k=1:n
    if (te(q)<tem(k))
        if (f>0)
        j=k;
        f=0;
        end
    end
end
if ((j==1) && (f==0)) 
    j=2;
end
if ((j==1) && (f==1))
    j=n;
end
ep(q)=epra(j-1)+(epra(j)-epra(j-1))*(te(q)-tem(j-1))/(tem(j)-tem(j-1));
end
koe=ep;
end