function ref = kootsre(dv,refl,T)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*c0^2;
c2=PP*c0/PB;
I0=0; ref1=0;
for j=1:length(dv)
    la=dv(j);
I0(j)=c1/((la^5)*(exp(c2/(la*T))-1));
ref1(j)=refl(j)*I0(j);
end
ref=trapz(dv,ref1)/trapz(dv,I0);
end