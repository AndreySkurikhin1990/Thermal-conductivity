%аппроксимация температурного поля, проверка закона СБ
ep=0; lambda=0; Sp=0;  knu=0;
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=PP*c0^2;
c2=PP*c0/PB;
sig=2*c1*(pi^5)/(15*(c2^4));
n=Kramers_n();
tol=3e4; t0=273.15; T0=563+t0; TH=141+t0; Nt=3e1; 
Sp=dlvoVer53101; leSp=length(Sp); knu=1e6*SredGraf; nu=0;
for  k=1:Nt+1 
    tem(k)=((T0^4+(TH^4-T0^4)*(k-1)/Nt)^0.25); 
end
for k=1:leSp
Sp(k)=1e-2/Sp(k);
end
for p=1:Nt+1
T=tem(p);
ept=0; Ib=0;
    for k=1:leSp
        nu(k)=c0/Sp(k);
me=exp((PP*nu(k))/(PB*T))-1;
Ib(k)=(2*PP*(nu(k)^3)/(c0^2))/me;
ept(k)=Ib(k)*knu(k);
    end
    otn(p)=trapz(Sp,Ib);
    ep(p)=trapz(Sp,ept)/otn(p);
    otn(p)=otn(p)/(sig*(T^4));
end
dl=1e-2:1e-3:5e3; p=length(dl); T=T0; Ib=0; dl=1e-6*dl; nuf=0;
for k=1:p
    nuf(k)=c0/dl(k);
    %if (k<10) disp(nuf(k)); end;
me=exp((PP*nuf(k))/(PB*T))-1;
Ib(k)=2*PP*(nuf(k)^3)/(c0^2);
Ib(k)=Ib(k)/me;
end
disp(trapz(nuf,Ib)/(sig*(T^4)));
%q=plot(tem-t0,otn,'-b');set(q,'LineWidth',3);xlabel('Температура, °С');ylabel('Отношение');title('График epsilon(t)');hold on; grid on;