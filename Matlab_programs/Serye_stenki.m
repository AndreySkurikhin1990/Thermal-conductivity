n=Kramers_n();
PP=6.6260755e-34; PB=1.380658e-23; c0=299792458; 
c1=2*1e24*pi*PP*c0^2; c2=PP*1e6*c0/PB; 
tol=3e4; t0=273.15; 
T0=927+t0; 
TH=369+t0; 
Nt=20; 
Nit=1;
koo=0:tol/Nt:tol; 
Sp=dlvoVer53101; 
leSp=length(Sp); 
knuSha=preobMasSha_n(); 
knu=1e6*SredGraf();
for  k=1:Nt+1 
    tem(k)=((T0^4+(TH^4-T0^4)*(k-1)/Nt)^0.25); 
end
%nsh=Kramers_n_Sha(); 
nsh=preobMasSha_n();
nve=Kramers_n();
for k=1:leSp
    Sp(k)=1e4/Sp(k);
    hiSha(k)=Sp(k)*knuSha(k)/(4*pi); 
    hi(k)=Sp(k)*knu(k)/(4*pi);
    nsha(k)=((nsh(k)^2+hiSha(k)^2)/(nve(k)^2+hi(k)^2))^(1/2);
ronu(k)=(1/2)+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2);
ronu(k)=ronu(k)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1));
ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));
ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log((nsha(k)-1)/(nsha(k)+1))/((nsha(k)^2+1)^3);
end;
for q=1:Nit
for p=1:Nt+1
x=koo(p);
    for k=1:leSp
lambda=Sp(k);
Ib0=c1/((lambda^5)*(exp(c2/(lambda*T0))-1));
IbH=c1/((lambda^5)*(exp(c2/(lambda*TH))-1));
f1(k)=(n(k)^2)*Ib0*inteMuSerSt(knuSha,ronu,tol,x);
f2(k)=IbH*(n(k)^2)*inteMuSerSt2(knuSha,ronu,tol,x);
ksi=0; E2=0; ff3=0; Ibb=0; jnu=0; t=1;
for j=1:p-1
    ksi(j)=koo(j);
    Ibb=c1/((lambda^5)*(exp(c2/(lambda*tem(j)))-1));
    jnu(j)=knu(k)*(n(k)^2)*Ibb;
    g=knu(k)*(x-ksi(j)); 
    if (g==0) 
        E2(j)=1; 
    else E2(j)=exp(-g)/g;
    if (E2(j)<1e-12) 
        E2(j)=0; 
    else
        E2(j)=integroexpon(2,g);
    end; 
    end;
    ff3(t)=jnu(j)*E2(j);
    t=t+1;
end
m=0; 
    if (p>1)
    for q=1:(length(ksi)-1) 
        m=m+(ksi(q+1)-ksi(q))*(ff3(q)+ff3(q+1))/2;
    end; 
    end;
        f3(k)=m;
ksis=0; E2=0; ff4=0; jnu=0; t=1;
    for j=p:Nt+1
    ksis(t)=koo(j);
    Ibb=c1/((lambda^5)*(exp(c2/(lambda*tem(j)))-1));
    jnu(j)=knu(k)*(n(k)^2)*Ibb;
    g=knu(k)*(ksis(t)-x); 
    if (g==0) 
        E2(j)=1; 
    else E2(j)=exp(-g)/g;
    if (E2(j)<1e-12) 
        E2(j)=0; 
    else E2(j)=integroexpon(2,g); 
    end; 
    end;
    ff4(t)=jnu(j)*E2(j);
    t=t+1;
end
m=0;
if (p<Nt)
    for j=1:(length(ff4)-1)
        m=m+(ksis(j+1)-ksis(j))*(ff4(j)+ff4(j+1))/2;
    end; end;
    f4(k)=m;
ff5=0; jnu=0;
    for j=1:Nt+1
    Ibb=c1/((lambda^5)*(exp(c2/(lambda*tem(j)))-1));
    jnu(j)=knu(k)*(n(k)^2)*Ibb;
    ff5(j)=jnu(j)*inteMuSerSt3(knuSha,ronu,tol,x,koo(j));
    end
    f5(k)=trapz(koo,ff5);
    f0(k)=f1(k)-f2(k)+f3(k)-f4(k)+f5(k);
    end
    E(p)=2*pi*trapz(Sp,f0)*1e-6;
    disp(p);
end
%Et=0;
%for j=1:length(E)-1
%    Et(j)=(E(j)+E(j+1))/2;
%end
%for j=1:length(Et)
%    E(j)=Et(j);
%end
end
b=plot(koo,E,'-b');
set(b,'LineWidth',3);
xlabel('Координата, мкм');
ylabel('Плотность теплового потока, Вт/м2');
title('График E(x)');