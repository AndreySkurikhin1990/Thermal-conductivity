function knus = knusreddv(knu,npp,dv,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c10=PP*c0^2;
c20=PP*c0/PB;
c1=c10;
c2=c20;
for k=1:length(npp)
    c1=c1/(npp(k)^2);
    c2=c2/npp(k);
    lambda=dv(k)/npp(k);
    dl(k)=lambda;
    me=(exp(c2/lambda/tem)-1);
    Ib(k)=2*pi*c1/(lambda^5)/me;
Ibc(k)=knu(k)*Ib(k);
c1=c10;
c2=c20;
end
nc=trapz(dl,Ibc);
nz=trapz(dl,Ib);
knus=nc/nz;
end