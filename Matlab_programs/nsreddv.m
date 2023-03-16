function ns = nsreddv(dv,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=PP*c0^2;
c2=PP*c0/PB;
c01=c1;
c02=c2;
for k=1:length(npp)
    lambda=dv(k);
    lambda=dv(k)/npp(k);
    c1=c1/(npp(k)^2);
    c2=c2/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
    %Ib(k)=Ib(k)/(npp(k)^2);    
Ibn(k)=npp(k)*Ib(k);
c1=c01;
c2=c02;
end
nc=trapz(dv,Ibn);
ns=abs(real(nc/trapz(dv,Ib)));
end