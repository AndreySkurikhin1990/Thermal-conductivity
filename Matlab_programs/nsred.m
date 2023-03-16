function ns = nsred(nu,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
for k=1:length(npp)
me=exp((PP*nu(k))/(PB*tem))-1;
Ib(k)=(2*PP*(nu(k)^3)/(c0^2))/me;
Ibn(k)=npp(k)*Ib(k);
end
nc=trapz(nu,Ibn);
ns=nc/trapz(nu,Ib);
end