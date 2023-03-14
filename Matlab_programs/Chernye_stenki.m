%ng=(1.545+1.581+1.538+1.542+1.534+1.553)/6; np=(1.525+1.561+1.51+1.514+1.516+1.536)/6; nm=ng; n=(nm+ng+np)/3;
n=Kramers_n();
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=PP*c0^2;
c2=PP*c0/PB;
sig=2*c1*(pi^5)/(15*(c2^4));
tol=3e4; t0=273.15; T0=563+t0; TH=141+t0; Nt=2e1;
koo=0:tol/Nt:tol; koo=(1e-6)*koo;
Sp=dlvoVer53101; leSp=length(Sp); knu=(1e6)*SredGraf;
for  k=1:Nt+1 
    tem(k)=((T0^4+(TH^4-T0^4)*(k-1)/Nt)^0.25); end;
for k=1:leSp    
Sp(k)=(1e-2)/Sp(k); end;
for p=1:Nt+1
x=koo(p);
    for k=1:leSp
lambda=Sp(k);
Ib0=c1/((lambda^5)*(exp(c2/(lambda*T0))-1));
IbH=c1/((lambda^5)*(exp(c2/(lambda*TH))-1));
g=knu(k)*x; 
if (g==0) 
    E3=0.5; else 
    E3=exp(-g)/g;
if (E3<1e-12) 
    E3=0; else 
    E3=integroexpon(3,g); end; end;
f1(k)=(n(k)^2)*Ib0*E3;
g=knu(k)*(tol-x); 
if (g==0) E3=0.5; else E3=exp(-g)/g;
if (E3<1e-12) E3=0; else E3=integroexpon(3,g); end; end;
f2(k)=IbH*E3*n(k)^2;
ksi=0; E2=0; ff3=0; Ibb=0; jnu=0; t=1;
for j=1:p-1
    ksi(j)=koo(j);
    Ibb=c1/((lambda^5)*(exp(c2/(lambda*tem(j)))-1));
    jnu(j)=knu(k)*(n(k)^2)*Ibb;
    g=knu(k)*(x-ksi(j)); 
    if (g==0) 
        E2(j)=1; else 
        E2(j)=exp(-g)/g;
    if (E2(j)<1e-12) 
        E2(j)=0; else 
        E2(j)=integroexpon(2,g); end; end;
    ff3(t)=jnu(j)*E2(j);
    t=t+1;
end
m=0;
for j=2:t-1
    m=m+(ksi(j)-ksi(j-1))*(ff3(j-1)+ff3(j))/2;
end
    f3(k)=m;
ksis=0; E2=0; ff4=0; jnu=0; t=1;
    for j=p:Nt
    ksis(t)=koo(j);
    Ibb=c1/((lambda^5)*(exp(c2/(lambda*tem(j)))-1));
    jnu(j)=knu(k)*(n(k)^2)*Ibb;
    g=knu(k)*(ksis(t)-x); 
    if (g==0) 
        E2(j)=1; else 
        E2(j)=exp(-g)/g;
    if (E2(j)<1e-12) 
        E2(j)=0; else 
        E2(j)=integroexpon(2,g); end; end;
    ff4(t)=jnu(j)*E2(j);
    t=t+1;
    end
m=0;
for j=2:t-1
    m=m+(ksis(j)-ksis(j-1))*(ff4(j-1)+ff4(j))/2;
end
    f4(k)=m;
    f0(k)=2*pi*(f1(k)-f2(k)+f3(k)-f4(k));
    end
    E(p)=2*pi*trapz(Sp,f0);
end
koo=(1e3)*koo;
q=plot(koo,E,'-b');
set(q,'LineWidth',3);
xlabel('Координата, мм');
ylabel('Плотность теплового потока, Вт/м2');
title('График E(x)');
hold on; grid on;