function inte1 = I1tau(tau,R1,R2,tau0)
hmu=1e-2; mu0=1e-2; muk=1;
mu=mu0:hmu:muk;
for k=1:length(mu)
I1=exp(-tau/mu(k));
I2=R2*exp(-(2*tau0-tau)/mu(k));
I0(k)=(I1+I2)*beta(tau0,mu(k),R1,R2)*(1-R1)*(mu(k)^2);
end
inte1=trapz(mu,I0);
end

function beta = rasbeta(tau0,mu,R1,R2)
beta=(1-R1*R2*exp(-2*tau0/mu))^-1;
end