function pr = tmp20()
y0=30e-3; d=y0; te0=273.15; ep=1e-20; temvh=arrTemHigh(); temvh=temvh+te0; temvc=arrTemCold(); temvc=temvc+te0; 
tepo=arrTepPot84(); qv=(1e4/13.85)*tepo; pt=length(temvc); ts=(temvc+temvh)/2;
tepv=koeftepv(qv,y0,temvh,temvc,pt);
w=15; Tk=temvc(w); T0=temvh(w); Tz=T0;
ts1=ts(15);
ts2=ts(16);
koeftep1=tepv(15);
koeftep2=tepv(16);
kktp1=(koeftep2-koeftep1)/(ts2-ts1); kktp2=koeftep1-kktp1*ts1; 
k1=abs(T0-Tk)/d; qv=koeftep1*k1;
koeftep3=kktp1*T0+kktp2;
k1=qv/koeftep3;
tol=30e-3
%Tk=T0-abs(k1*tol)
sig=5.67e-8;
npp=Kramers_n();
dl=RasshDiapDlinVoln();
p=length(dl);
Nt=p-1;
alfs=RasMasKoAbs();
tau=0; 
hksi=tol/Nt; 
ksi=0:hksi:tol; 
k1=(Tk-T0)/tol
k2=T0
for k=1:p
    temp(k)=k1*ksi(k)+k2;
    dv(k)=dl(k)/npp(k);
end
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=2*pi*PP*(c0^2);
C2=PP*c0/PB;
c1t=C1;
c2t=C2;
%1----------
for j=1:p
it1=0; eb1=0;
for k=1:p
    tem=temp(k);
    c1t=c1t/(npp(k)^2);
    c2t=c2t/npp(k);
lam=dl(k)/npp(k);
me=exp(c2t/lam/tem)-1;
eb1(k)=alfs(k)*c1t/(lam^5)/me;
it1(k)=eb1(k)*ProvAd(integroexpon(2,alfs(k)*(tol-ksi(j))));
c1t=C1;
c2t=C2;
end
if (rem(j,1e3)==0)
    disp(j);
end
I1g(j)=inte2ma(dv,it1);
end
I1gg=2*pi*inte2ma(ksi,I1g)
%eb1=0; I1g=0;
%for j=1:p
%   tem=temp(j);
%   c1t=c1t/(npp(j)^2);
%   c2t=c2t/npp(j);
%   lam=dl(j)/npp(j);
%   me=exp(c2t/lam/tem)-1;
%   eb1(j)=alfs(j)*c1t/(lam^5)/me;
%   it1=0; 
%   for k=1:p
%it1(k)=eb1(j)*ProvAd(integroexpon(2,alfs(j)*(tol-ksi(k))));
    %end
%c1t=C1;
%c2t=C2;
%if (rem(j,1e3)==0)
%    disp(j);
%end
%I1g(j)=inte2ma(ksi,it1);
%end
%I1gg=2*pi*inte2ma(dv,I1g)
%2----------
%I1qr=0;
%for j=1:p
%it2=0; eb2=0;
%for k=1:p
    %tem=temp(k);
    %c1t=c1t/(npp(k)^2);
    %c2t=c2t/npp(k);
    %lam=dl(k)/npp(k);
    %me=exp(c2t/lam/tem)-1;
    %eb2(k)=c1t/(lam^5)/me;
    %it2(k)=eb2(k)*alfs(k);
%c1t=C1;
%c2t=C2;
%end
%if (rem(j,1e3)==0)
    %disp(j);
%end
%I1qr(j)=inte2ma(dv,it2);
%end
%I1q=4*pi*inte2ma(ksi,I1qr)
%3----------
%for j=1:p
%it3=0; eb3=0;
%for k=1:p
    %tem=temp(k);
    %c1t=c1t/(npp(k)^2);
    %c2t=c2t/npp(k);
    %lam=dl(k)/npp(k);
    %me=exp(c2t/lam/tem)-1;
    %eb3(k)=c1t/(lam^5)/me;
    %it3(k)=eb3(k)*alfs(k)*ProvAd(integroexpon(2,alfs(k)*ksi(j)));
%c1t=C1;
%c2t=C2;
%end
%I3g(j)=inte2ma(dv,it3);
%if (rem(j,1e3)==0)
    %disp(j);
%end
%end
%I3gg=2*pi*inte2ma(ksi,I3g)
it3=0; eb3=0;
for j=1:p
    tem=temp(j);
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    lam=dl(j)/npp(j);
    me=exp(c2t/lam/tem)-1;
    eb3(j)=c1t/(lam^5)/me;
for k=1:p
    it3(k)=eb3(j)*alfs(j)*ProvAd(integroexpon(2,alfs(j)*ksi(k)));
end
if (rem(j,1e3)==0)
    disp(j);
end
c1t=C1;
c2t=C2;
I3g(j)=inte2ma(ksi,it3);
end
I3gg=2*pi*inte2ma(dv,I3g)
%----------------------
end

function inte = inte2ma(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function mk = ProvAd(m)
ep=1e-40;
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
if (abs(m)<ep)  
    m=0;  
end; 
mk=real(abs(m));
end