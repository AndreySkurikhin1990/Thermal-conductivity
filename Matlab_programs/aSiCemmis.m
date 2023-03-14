function [ asic ] = aSiCemmis(dv)
c0=299792458;
nu=c0/dv;
eps0=6.7;
nupi=4.327e13;
gai=1.428e11;
nui=2.38e13;
epss=eps0+(nui^2-nu^2)*nupi^2/((nui^2-nu^2)^2+(nu*gai)^2);
epsss=gai*nu*(nupi^2)/((nui^2-nu^2)^2+(nu*gai)^2);
n=0.5*(epss+sqrt(epss^2+epsss^2)); n=sqrt(n);
k=0.5*(-epss+sqrt(epss^2+epsss^2)); k=sqrt(k);
asic(1)=emmet(n,k);
asic(2)=n;
asic(3)=k;
end