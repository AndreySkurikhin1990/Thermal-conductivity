function t = tmp6()
format long g;
t=RasEpsilon();
end

function v = vychkoef()
te0=273.15;
%ktp=[0.2,0.28,0.68];
%te=[800,1000,1600];
ktp=[0.11,0.2,0.25,0.4]*4184/36e2;
te=[1,6,8,11]*1e2;
te=te+te0;
n=length(te);
ko=polyfit(te,ktp,n-1);
tem=35e1+te0;
%tem=6e2+te0;
v=polyval(ko,tem);
end

function v = vych()
format long g;
e2=0; e3=0; ko=0; N=6; p=2; x=2; h=0.25;
for k=1:N+1
e2(k)=integroexpon(p,x);
ko(k)=x;
x=x+h;
end
for k=1:N+1
    e3(k)=(exp(-ko(k))-ko(k)*e2(k))/2;
end
for k=1:N+1
    e4(k)=(exp(-ko(k))-ko(k)*e3(k))/3;
end
disp(e4);
disp(ko);
v=0;
end

function v = proNorRas()
sig=erfcinv(1e-4);
sig=2*sig;
ami=-3/sqrt(2);
ama=-ami;
P=erf(ama)-erf(ami);
P=P/2;
ama=60;
ami=0;
mo=(ami+ama)/2;
sig=abs(ama-ami)/sig/sqrt(2);
sig=sqrt(2)*sig
z1=(ami-mo)/sqrt(2)/sig;
z2=(ama-mo)/sqrt(2)/sig;
P=(erf(z2)-erf(z1))/2;
P=(1-P)*1e2;
v=P;
end

function p = RasEpsilon()
wmg=217e-3; wal=132e-3; wsi=368e-3;
r=[626,921,626,921,633,923,633,923,628,924,628,924,626,907,626,907,624,900,624,900,625,920,625,920]';
koal=rasKoefAlpha(r,wmg,wsi,wal);
koal=0;
r=[620,891,620,891,627,882,627,882,621,905,621,905,620,891,627,882,621,905];
koal=rasKoefAlpha(r,wmg,wsi,wal)
p=0;
end

function [ ko ] = rasKoefAlpha(te,wmg,wsi,wal)
Tna=3e2; Tko=16e2; delT=1e2; tem=Tna:delT:Tko;
mgo=[1,0.94,0.87,0.801,0.736,0.676,0.635,0.59,0.567,0.543,0.53,0.525,0.515,0.507];
al2o3=[1,0.98,0.946,0.898,0.854,0.797,0.753,0.709,0.676,0.642,0.618,0.598,0.582,0.57];
sio2=[1,0.984,0.953,0.917,0.854,0.808,0.756,0.711,0.578,0.523,0.495,0.468,0.448,0.429];
wo=(wmg+wsi+wal);
kumm=vychKoe(tem,te,mgo);
kuma=vychKoe(tem,te,al2o3);
kums=vychKoe(tem,te,sio2);
n=length(te);
for k=1:n
    koo(k)=kumm(k)*wmg/wo+kuma(k)*wal/wo+kums(k)*wsi/wo;
end
ko=koo';
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