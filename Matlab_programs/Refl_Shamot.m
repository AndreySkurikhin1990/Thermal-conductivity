function kootr = Refl_Shamot(tem)
npp=Kramers_n_Sha_uk();
dl=dliVol();
knuSha=preobMasSha_n(); 
alfs=1e6*SredGraf();
nsh=Kramers_n_Sha(); 
p=length(npp);
for k=1:p  
    hiSha(k)=abs(dl(k)*knuSha(k))/(4*pi);
    hi(k)=abs(dl(k)*alfs(k))/(4*pi);
    npp(k)=abs(npp(k)); nsh(k)=abs(nsh(k)); 
    nsha(k)=sqrt((nsh(k)^2+hiSha(k)^2)/(npp(k)^2+hi(k)^2)); nsha(k)=abs(nsha(k));
    ronu(k)=(5e-1)+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2);
    ronu(k)=ronu(k)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); 
    podln=(nsha(k)-1)/(nsha(k)+1);
    podln=abs(podln);
    ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));
    ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log(podln)/((nsha(k)^2+1)^3);
    ronu(k)=abs(ronu(k));
    if (ronu(k)>1)
        ronu(k)=1;
    end
end
Refl=knusreddv(ronu,nsha,dl,tem)
Refl=Refsred(ronu,nsha,dl,tem)
Refl=Reflsreddvvsr(ronu,nsha,dl,tem)
kootr=Refl;
end

function [ dv ] = dliVol()
dl=dlvoVer53101();
p=length(dl);
for k=1:p
    dl(k)=1e-2/dl(k);
end
dv=dl;
end

function knus = Reflsreddvvsr(ronu,nsha,dl,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=PP*c0^2;
c2=PP*c0/PB;
ct1=c1; ct2=c2;
for k=1:length(nsha)
    ct1=ct1/(nsha(k)^2);
    ct2=ct2/nsha(k);
    lambda=dl(k)/nsha(k);
    dv(k)=lambda;
    me=exp(ct2/lambda/tem)-1;
    Ib(k)=2*pi*ct1/(lambda^5)/me;
    Ibc(k)=ronu(k)*Ib(k);
ct1=c1;
ct2=c2;
end
nc=trapz(dv,Ibc);
nz=trapz(dv,Ib);
knus=nc/nz;
end