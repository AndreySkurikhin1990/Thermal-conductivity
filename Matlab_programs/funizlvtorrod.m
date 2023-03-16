function fivr = funizlvtorrod(npp,lam,T)
%format long g;
pp=6.6260755e-34;
pb=1.380658e-23;
c0=299792458;
c1=2*pi*pp*c0^2;
c2=pp*c0/pb;
sig=(c1*pi^4)/(15*c2^4);
lt=0;
lt=lam*T;
N=1e4;
    ltm=0;
    c1m=0;
    c2m=0;
    F0ltm=0;
for j=1:length(lt)
ltm(j)=lt(j)/npp(j);
c1m(j)=c1/(npp(j)^2);
c2m(j)=c2/npp(j);
    F0ltm(j)=c1m(j)*(T^4)/((ltm(j)^5)*(exp(c2m(j)/ltm(j))-1));
end
F0=trapz(ltm,F0ltm);
for t=1:length(lt)
    ltm=0;
    c1m=0;
    c2m=0;
    F0ltm=0;
    if (t==1) 
        p=2; else p=t; end;
for j=1:p
    %s=0;
    %dze=c2/lt(j);
%for k=1:N
    %s=s+(npp(j)^2)*exp(-k*dze)*(dze^3+3*dze^2/k+6*dze/k^2+6/k^3)/k;
%end
%F0lt(j)=(15/(pi^4))*s/(npp(j)^2);
ltm(j)=lt(j)/npp(j);
c1m(j)=c1/(npp(j)^2);
c2m(j)=c2/npp(j);
F0ltm(j)=c1m(j)*(T^4)/((ltm(j)^5)*(exp(c2m(j)/ltm(j))-1));
end
F0lt(t)=trapz(ltm,F0ltm)/F0;
end
fz0lt=0;
for j=1:length(F0lt)
ltm(j)=lt(j)/npp(j);
c1m(j)=c1/(npp(j)^2);
c2m(j)=c2/npp(j);
    fz0lt(j)=F0lt(j)+(ltm(j)/4)*(c1m(j)*(T^4)/((ltm(j)^5)*(exp(c2m(j)/ltm(j))-1)))/F0;
end
%for k=1:length(lt)
    %np=1; 
%    np=npp(k);
%    dze=c2/(lt(k)/np); fslt=(15/((pi^4)*c2))*(dze^5/(exp(dze)-1));
%fz0lt(k)=F0lt(k)+(lt(k)/4/np)*fslt;
%end
%Ibb=0;
%for k=1:length(lam)
%    lambda=lam(k);
%    Ib(k)=c1/((lambda^5)*(exp(c2/(lambda*T))-1));
%end
%I=trapz(lam,Ib);
%disp(I/(sig*(T^4)));
fivr=fz0lt;
end