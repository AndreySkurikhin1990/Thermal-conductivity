%считает среднее значение коэффициента отражения вермикулита
function t = tmp50()
npp=0; npp=Kramers_n();
dl=RasshDiapDlinVoln();
p=length(npp);
for k=1:p
    npp(k)=abs(npp(k));
    ronu(k)=0.5+((npp(k)-1)*(3*npp(k)+1))/(6*(npp(k)+1)^2);
    ronu(k)=ronu(k)-(2*(npp(k)^3)*(npp(k)^2+2*npp(k)-1))/((npp(k)^2+1)*(npp(k)^4-1)); 
    podln=(npp(k)-1)/(npp(k)+1);
    podln=abs(podln);
    ronu(k)=ronu(k)+(8*(npp(k)^4)*((npp(k)^4)+1)*log(npp(k)))/((npp(k)^2+1)*((npp(k)^4-1)^2));
    ronu(k)=ronu(k)+(npp(k)^2)*((npp(k)^2-1)^2)*log(podln)/((npp(k)^2+1)^3);
end;
%kootsre(dv,refl,T)
tem=1e2:1e2:12e2; tem=tem+273.15;
for k=1:length(tem)
Ref(k)=Refsred(ronu,npp,dl,tem(k));
end
Ref=Ref'
tem=tem'
t=0;
end