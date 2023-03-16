function per1 = I1tauUpr(tau,R1,R2,tau0)
I=0;
I(1)=ProvAdekv(integroexpon(4,tau));
I(2)=R2*ProvAdekv(integroexpon(4,(2*tau0-tau)));
I(3)=R1*R2*ProvAdekv(integroexpon(4,(2*tau0+tau)));
I(4)=R1*(R2^2)*ProvAdekv(integroexpon(4,(4*tau0-tau)));
I(5)=((R1*R2)^2)*ProvAdekv(integroexpon(4,(4*tau0+tau)));
I(6)=((R1*R2)^2)*R2*ProvAdekv(integroexpon(4,(6*tau0-tau)));
I(7)=((R1*R2)^3)*ProvAdekv(integroexpon(4,(6*tau0+tau)));
I(8)=((R1*R2)^3)*R2*ProvAdekv(integroexpon(4,(8*tau0-tau)));
I0=0;
for k=1:length(I)
    I0=I0+I(k);
end
per1=I0*(1-R1);
end