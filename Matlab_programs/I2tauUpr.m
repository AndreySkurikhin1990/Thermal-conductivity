function per2 = I2tauUpr(tau,R1,R2,tau0)
I=0;
I(1)=ProvAdekv(integroexpon(4,tau0-tau));
I(2)=R1*ProvAdekv(integroexpon(4,tau0+tau));
I(3)=R1*R2*ProvAdekv(integroexpon(4,3*tau0-tau));
I(4)=R2*(R1^2)*ProvAdekv(integroexpon(4,3*tau0+tau));
I(5)=((R1*R2)^2)*ProvAdekv(integroexpon(4,5*tau0-tau));
I(6)=((R1*R2)^2)*R1*ProvAdekv(integroexpon(4,5*tau0+tau));
I(7)=((R1*R2)^3)*ProvAdekv(integroexpon(4,7*tau0-tau));
I(8)=((R1*R2)^3)*R1*ProvAdekv(integroexpon(4,7*tau0+tau));
I0=0;
for k=1:length(I)
I0=I0+I(k);
end
per2=I0*(1-R2);
end