function knus = knusreddvvsr(knu,npp,dv,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=PP*c0^2;
c2=PP*c0/PB;
ct1=c1; ct2=c2; dl=0;
for k=1:length(npp)
    ct1=ct1/(npp(k)^2);
    ct2=ct2/npp(k);
    lambda=dv(k)/npp(k);
    dl(k)=lambda;
    Ib(k)=2*pi*ct1/((lambda^5)*(exp(ct2/(lambda*tem))-1));
    Ibc(k)=knu(k)*Ib(k);
ct1=c1;
ct2=c2;
end
nc=trapz(dl,Ibc);
nz=trapz(dl,Ib);
knus=nc/nz;
end