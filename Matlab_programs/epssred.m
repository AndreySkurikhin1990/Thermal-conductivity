function es = epssred(dv,eps,T,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*c0^2;
c2=PP*c0/PB;
iz0=0; iz1=0; dl=0;
c1t=c1; c2t=c2;
for j=1:length(dv)
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
iz0(j)=ProvAdekv(c1t/((la^5)*(exp(c2t/(la*T))-1)));
iz1(j)=eps(j)*iz0(j);
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
es=ProvAdekv(chi)/ProvAdekv(zna);
end