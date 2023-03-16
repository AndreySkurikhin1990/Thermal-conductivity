function intmu = inteMuSerSt(knu,rost,tol,x)
epnu=0;
for k=1:length(rost)
    epnu(k)=1-rost(k);
end;
Nt=1e3;
mu=0:1/Nt:1;
f0=0;
for k=2:length(mu)
    beta(k)=(1-(rost(k)^2)*exp(-2*knu(k)*tol/mu(k)))^(-1);
    f0(k)=exp(-knu(k)*x/mu(k))-rost(k)*exp(-knu(k)*(2*tol-x)/mu(k));
    f0(k)=epnu(k)*beta(k)*f0(k)*mu(k);
end;
f0(1)=0;
if (x==0) 
    f0(1)=1;
end
intmu=trapz(mu,f0);
end