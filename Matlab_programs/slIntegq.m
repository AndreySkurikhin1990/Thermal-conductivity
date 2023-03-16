function inte = slIntegq(alp,L,r1,r2,a1,a2,x,Tw1,Tw2,h)
sig=5.668e-8;
E3=integroexpon(3,alp*L);
tm=1-r1*r2*4*(E3^2);
tr=integroexpon(3,alp*x);
if (x==0) 
    tr=0.5; end
tu=integroexpon(3,alp*(L-x));
if (x==L) 
    tu=0.5; end
psi1=(2*tr-4*E3*tu*r2)/tm;
psi2=(2*tu-4*E3*tr*r1)/tm;
qr=psi1*a1*sig*Tw1^4-psi2*a2*sig*Tw2^4;
ksi=0; ko=0; pod=0; q=1;
while (ksi<=x)
tm=integroexpon(2,alp*ksi);
if (ksi==0) 
    tm=1; end
tr=integroexpon(2,alp*(x-ksi));
if (ksi==x) 
    tr=1; end
tu=integroexpon(2,alp*(L-ksi));
if (ksi==L) 
    tu=1; end
fib1=tr+r1*psi1*tm-r2*psi2*tu;
pod(q)=fib1*temp(q)^4;
ko(q)=ksi;
ksi=ksi+h;
q=q+1;
end
pod=(sig*alp*2)*pod;
integ1=trapz(ko,pod);
ksi=x; ko=0; pod=0; q=1;
while (ksi<=L)
tm=integroexpon(2,alp*ksi);
if (ksi==0) 
    tm=1; end
tr=integroexpon(2,alp*(ksi-x));
if (ksi==x) 
    tr=1; end
tu=integroexpon(2,alp*(L-ksi));
if (ksi==L) 
    tu=1; end
fib2=tr+r2*psi2*tu-r1*psi1*tm;
pod(q)=fib2*temp(q)^4;
ko(q)=ksi;
ksi=ksi+h;
q=q+1;
end
pod=(sig*alp*2)*pod;
integ2=trapz(ko,pod);
qr=integ1-integ2+qr;
inte=qr;
end