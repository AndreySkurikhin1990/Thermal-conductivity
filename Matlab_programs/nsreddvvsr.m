function ns = nsreddvvsr(dv,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
vl=c0;
for k=1:length(npp)
    vl=vl/npp(k);
    c1=PP*vl^2;
    c2=PP*vl/PB;
    %lambda=dv(k);
    lambda=dv(k)/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
    %Ib(k)=Ib(k)/(npp(k)^2);    
Ibn(k)=(npp(k)^2)*Ib(k);
vl=c0;
end
nc=trapz(dv,Ibn);
ns=abs(real(nc/trapz(dv,Ib)));
end