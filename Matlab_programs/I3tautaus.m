function inte3 = I3tautaus(taus,tau,R1,R2,tau0)
hmu=1e-2; mu0=1e-2; muk=1;
mu=mu0:hmu:muk;
for k=1:length(mu)
I1=exp(-(tau+taus)/mu(k));
I2=R1*R2*exp(-(2*tau0+tau-taus)/mu(k));
I3=R1*R2*exp(-(2*tau0-tau+taus)/mu(k));
I4=R2*exp(-(2*tau0-tau-taus)/mu(k));
I0(k)=(I1+I2+I3+I4)*beta(tau0,mu(k),R1,R2)*mu(k);
end
inte3=trapz(mu,I0);
end

function beta = rasbeta(tau0,mu,R1,R2)
beta=(1-R1*R2*exp(-2*tau0/mu))^-1;
end