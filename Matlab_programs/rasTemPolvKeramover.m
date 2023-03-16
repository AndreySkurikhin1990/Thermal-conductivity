function [ temvy ] = rasTemPolvKeramover(tau0,koeftep,R1,R2,np,kopove,Tna,Tkon,tol,Q)
sig=5.67e-8;
R1
R2
tau0
koeftep
np
kopove
Tna
Tkon
tol
Q
Nt=5e1-1; htau=tau0/Nt; tau=0:htau:tau0; taus=tau;
temp=nachPribTemN(Tna,Tkon,tol/Nt,tol,np,koeftep,kopove)
for ite=1:4
for k=1:length(tau) 
rkp=2*(np^2)*sig/kopove/koeftep; inttau=0;
for j=1:length(taus)
inc1=ProvAdekv(integroexpon(3,taus(j)));
inc2=ProvAdekv(integroexpon(3,abs(tau(k)-taus(j))));
inc3=I3tautausUpr(0,taus(j),R1,R2,tau0);
inc4=I3tautausUpr(tau(k),taus(j),R1,R2,tau0);
inttau(j)=(inc1-inc2+inc3-inc4)*(temp(j)^4);
inttau(j)=ProvAd(inttau(j), 0, 0, 0);
end
inta=trapz(taus,inttau);
inta0=I1tauUpr(0,R1,R2,tau0)-I1tauUpr(tau(k),R1,R2,tau0);
inta=inta+(Tna^4)*inta0;
inta0=I2tauUpr(0,R1,R2,tau0)-I2tauUpr(tau(k),R1,R2,tau0);
inta=inta-(Tkon^4)*inta0;
inta=rkp*inta;
inta=inta+Tna-Q*tau(k)/koeftep/kopove;
temp(k)=ProvAd(inta, Tna, Tkon, 1);
temp2=overrecarrto2darr(ite,temp,temp,length(temp));
disp(k);
end
end
temvy=temp2;
end

function mk = ProvAd(m, Tn, Tk, no)
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end
if (no==1)
if (m>Tk)  
    m=Tk;  
end
if (m<Tn)  
    m=Tn;  
end
end
mk=abs(m);
end