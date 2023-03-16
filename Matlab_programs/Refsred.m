function refs = Refsred(ronu,npp,dl,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=c0;
for k=1:length(npp)
nu(k)=c0/dl(k);
me=(PP*nu(k))/(PB*tem);
me=exp(me)-1;
c1=c1/npp(k);
Ib(k)=(2*PP*(nu(k)^3)/(c1^2))/me;
Ibc(k)=ronu(k)*Ib(k);
c1=c0;
end
nc=trapz(nu,Ibc);
nz=trapz(nu,Ib);
refs=nc/nz;
end