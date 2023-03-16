function knus = knusreddvvsrSha(tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=PP*c0^2; 
ct1=c1;
c2=PP*c0/PB; 
ct2=c2;
dv=dliVol();
knuSha=preobMasSha_n(); 
npp=Kramers_n_Sha();
for k=1:length(npp)
    ct1=ct1/(npp(k)^2);
    ct2=ct2/npp(k);
    lambda=dv(k)/npp(k);
    dl(k)=lambda;
    Ib(k)=2*pi*ct1/((lambda^5)*(exp(ct2/(lambda*tem))-1));
    Ibc(k)=knuSha(k)*Ib(k);
ct1=c1;
ct2=c2;
end
nc=trapz(dl,Ibc);
nz=trapz(dl,Ib);
knus=nc/nz;
end

function [ dv ] = dliVol()
dl=dlvoVer53101();
p=length(dl);
for k=1:p
    dl(k)=1e-2/dl(k);
end
dv=dl;
end