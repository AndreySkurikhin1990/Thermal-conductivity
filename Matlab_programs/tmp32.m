%определяет коэффициенты уменьшения степени черноты с температурой
function t = tmp32()
te0=273.15;
t3=arrTem3VVF1(); t3=te0+t3;
t1=arrTem1VVF1(); t1=te0+t1;
t2=arrTem2VVF1(); t2=te0+t2;
ktp=arrKTP_VVF1(); ts=(t1+t2)/2;
no=5;
ktp1=ktp(no)
Tk=t2(no)
T0=t1(no)
Tsr=ts(no)
Tsred=t3(no)
tol=30e-3
qv=ktp1*(T0-Tk)/tol
npp=Kramers_n();
dl=RasshDiapDlinVoln();
p=length(dl);
Nt=2e1
alfs=RasMasKoAbs();
hksi=tol/Nt; 
ksi=0:hksi:tol; 
Nt=length(ksi);
for k=1:p
    dv(k)=dl(k)/npp(k);
end
mko=koefOprTemp(T0,Tk,Tsred);
at=mko(1);
bt=mko(2);
ct=mko(3);
for k=1:Nt
    temp(k)=at*(ksi(k)^2)+bt*ksi(k)+ct;
end
koal=rasKoefAlpha(temp,npp);
t=0;
end

function [ ko ] = rasKoefAlpha(te,pp)
tem=3e2:1e2:16e2;
mgo = [1 0.94 0.87 0.801 0.736 0.676 0.635 0.59 ...
    0.567 0.543 0.53 0.525 0.515 0.507];
al2o3 = [1 0.98 0.946 0.898 0.854 0.797 0.753 ...
    0.709 0.676 0.642 0.618 0.598 0.582 0.57];
sio2 = [1 0.984 0.953 0.917 0.854 0.808 0.756 ...
    0.711 0.578 0.523 0.495 0.468 0.448 0.429];
wmg=20.44; wsi=36.01; wal=13.84; wo=70.29;
kumm=vychKoe(tem,te,mgo);
kuma=vychKoe(tem,te,al2o3);
kums=vychKoe(tem,te,sio2);
for k=1:length(pp)
    ep(k)=epsilam(pp(k));
end
n=length(te);
for k=1:n
    koo(k)=kumm(k)*wmg/wo+kuma(k)*wal/wo+kums(k)*wsi/wo;
end
ko=koo;
end

function epsi = epsilam(n)
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(abs(n))/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log(abs(n-1)/abs(n+1))*((n^2-1)^2)/((n^2+1)^3);
epsi=eps;
end

function [ koe ] = vychKoe(tem,te,epra)
for q=1:length(te)
    f=1;
for k=2:length(tem)
    if (te(q)<tem(k))
        if (f>0)
        j=k;
        f=0;
        end
    end
end
ep(q)=epra(j-1)+(epra(j)-epra(j-1))*(te(q)-tem(j-1))/(tem(j)-tem(j-1));
end
koe=ep;
end

function [ koe ] = koefOprTemp(T1,T2,T3)
c=T1;
x1=0; x3=15e-3; x2=30e-3;
x11=x2^2; x12=x2; x21=x3^2; x22=x3;
b1=T2-c; b2=T3-c;
de=x11*x22-x12*x21;
de1=b1*x22-b2*x12;
de2=x11*b2-x21*b1;
koe = [de1/de de2/de c];
end