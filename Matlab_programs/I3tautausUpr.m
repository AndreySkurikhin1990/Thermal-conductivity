function per3 = I3tautausUpr(tau,taus,R1,R2,tau0)
I=0; k=0;
k=k+1; I(k)=ProvAdekv(integroexpon(3,tau0+taus));
k=k+1; I(k)=R1*R2*ProvAdekv(integroexpon(3,2*tau0+tau-taus));
k=k+1; I(k)=R1*R2*ProvAdekv(integroexpon(3,2*tau0-tau+taus));
k=k+1; I(k)=R2*ProvAdekv(integroexpon(3,2*tau0-tau-taus));
k=k+1; I(k)=R1*R2*ProvAdekv(integroexpon(3,2*tau0+tau+taus));
k=k+1; I(k)=((R1*R2)^2)*ProvAdekv(integroexpon(3,4*tau0+tau-taus));
k=k+1; I(k)=((R1*R2)^2)*ProvAdekv(integroexpon(3,4*tau0-tau+taus));
k=k+1; I(k)=R1*(R2^2)*ProvAdekv(integroexpon(3,4*tau0-tau-taus));
k=k+1; I(k)=(R1*R2)*ProvAdekv(integroexpon(3,4*tau0+tau+taus));
k=k+1; I(k)=((R1*R2)^2)*ProvAdekv(integroexpon(3,6*tau0+tau-taus));
k=k+1; I(k)=((R1*R2)^2)*ProvAdekv(integroexpon(3,6*tau0-tau+taus));
k=k+1; I(k)=(R1*(R2^2))*ProvAdekv(integroexpon(3,6*tau0-tau-taus));
k=k+1; I(k)=((R1*R2)^3)*ProvAdekv(integroexpon(3,6*tau0+tau+taus));
k=k+1; I(k)=((R1*R2)^4)*ProvAdekv(integroexpon(3,8*tau0+tau-taus));
k=k+1; I(k)=((R1*R2)^4)*ProvAdekv(integroexpon(3,8*tau0-tau+taus));
k=k+1; I(k)=R2*((R1*R2)^4)*ProvAdekv(integroexpon(3,8*tau0-tau-taus));
I0=0;
for k=1:length(I)
I0=I0+I(k);
end
per3=I0;
end