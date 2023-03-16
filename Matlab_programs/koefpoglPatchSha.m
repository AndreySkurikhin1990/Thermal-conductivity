function kpPatch = koefpoglPatchSha(s,tem)
npp=0; npp=Kramers_n_Sha(); 
dl=0; dl=dlvoVer53101(); p=length(dl);
for k=1:p  
    dl(k)=1e-2/dl(k); 
end;
alfs=0; alfs=preobMasSha_n();
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=2*pi*PP*(c0^2);
C2=PP*c0/PB;
c1t=C1;
c2t=C2;
dlv=dl;
%pi=3.1415926535897932;
%sig=2*C1*(pi^5)/(15*(C2^4));
for k=1:p
c1t=c1t/(npp(k)^2);
c2t=c2t/npp(k);
lam=dl(k)/npp(k);
dlv(k)=lam;
me=(exp(c2t/lam/tem)-1);
eb(k)=c1t/(lam^5)/me;
ne=ProvAdekv(integroexpon(2,alfs(k)*s));
chi(k)=alfs(k)*eb(k)*ne;
zna(k)=eb(k)*ne;
c1t=C1;
c2t=C2;
end
chis=trapz(dlv,chi);
znam=trapz(dlv,zna);
kpPatch=chis/znam;
end