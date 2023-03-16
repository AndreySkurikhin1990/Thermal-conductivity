function rs = sredRosSieg(tem,npp,alfs,dl)
te0=273.15; 
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4));
np2=nsredPlank2(dl,npp,tem);
Ibc=0; Ibz=0;
for k=1:length(npp)
c1m(k)=C1/(npp(k)^2);
c2m(k)=C2/npp(k);
t=(pi/2)*(c1m(k)*c2m(k))/sig;
dlv(k)=dl(k)/npp(k);
chi=exp(c2m(k)/tem/dlv(k));
zna=(chi-1)^2;
Ibz(k)=t*chi/zna/(dlv(k)^6)/(tem^5)/np2; 
Ibc(k)=Ibz(k)/alfs(k);
end
chasc=trapz(dlv,Ibc);
chasz=trapz(dlv,Ibz);
dlsvprfo2=chasc/chasz;
rs=1/dlsvprfo2;
end

function ns = nsredPlank2(dv,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
for k=1:length(npp)
    vl=c0;
    vl=vl/npp(k);
    c1=PP*(vl^2);
    c2=PP*vl/PB;
    %lambda=dv(k);
    lambda=dv(k)/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
    %Ib(k)=Ib(k)/(npp(k)^2);    
Ibn(k)=(npp(k)^2)*Ib(k);
dvm(k)=lambda;
end
nc=trapz(dvm,Ibn);
nz=trapz(dvm,Ib);
ns=nc/nz;
end