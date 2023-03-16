te0=273.15; npp=0; npp=Kramers_n(); 
dl=0; dl=dlvoVer53101(); p=length(dl);
for k=1:p  
    dl(k)=1e-2/dl(k); 
end;
alfs=0; alfs=1e6*SredGraf; 
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4))
Cons=(pi/2)*(C1*C2)/sig;
temp=2e2:1e2:1e3; temp=temp+te0; 
tem=0; Ib=0; Ibc=0; Ibz=0; chasc=0; chacz=0; 
alf1=0; alf2=0; me=0; ote=0; chi=0; zna=0;
for n=1:length(temp)
    tem=temp(n);
for k=1:length(npp)
dlv(k)=dl(k)/npp(k);
dlv(k)=dl(k);
%me=exp(C2/(dlv(k)*tem))-1;Ib(k)=2*pi*C1/(dlv(k)^5)/me;ote=(sig/Ib(k))^(1/4);
eb=sig*(tem^4);
ote=(sig/eb)^(1/4);
chi=exp(C2*ote/dlv(k));
zna=(chi-1)^2;
Ibz(k)=Cons*(ote^5)*chi/zna/(dlv(k)^6); 
Ibc(k)=Ibz(k)/alfs(k);
end
tem=tem-te0
chasc=trapz(dlv,Ibc)*1e6
chasz=trapz(dlv,Ibz)*1e6
dlsvprfo1=sprko(tem+te0)
dlsvprfo2=chasc*1e6/chasz
end