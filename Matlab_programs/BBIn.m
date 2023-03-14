function frwl = BBIn(tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*c0^2;
C2=PP*c0/PB;
sig=2*C1*pi^5/(15*C2^4);
te0=273.15;
T=tem+te0;
koe=1e-6;
Sp=dlvoVer5311();
r=length(Sp);
for k=1:r 
    Sp(k)=1e-2/Sp(k); end;
als=SredGraf;
alsr=trapz(als,Sp);
alsr=alsr/(Sp(r)-Sp(1));
for k=1:r
    si(k)=2*pi*C1/((Sp(k)^5)*(exp(C2/(Sp(k)*T))-1));
    e(k)=si(k)*als(k);
    dlv(k)=Sp(k);
end
I=alsr*sig*(T^4);
for k=r:-1:1
    dll=0; ee=0;
    for j=1:k 
        dll(j)=Sp(j); 
        ee(j)=e(j); end;
    F=trapz(dll,ee);
    F=F/I;
        if (F<99e-2) 
            break; end;
end;
ndll=k;
%disp(r);
%disp(ndll);
%disp(F);
for k=2:ndll
    dlm=0; ec=0;
    for j=1:k 
        dlm(j)=Sp(j); 
        ec(j)=e(j);
    end
    F=trapz(dlm,ec);
    F=F/I;
    if (F>1e-2) 
            break; end;
end
%disp(F);
ndlm=k;
%disp(ndlm);
%j=1; for k=ndlm:ndll    dlr(j)=dlv(k);     j=j+1; end;
dlr(1)=ndlm;
dlr(2)=ndll;
%p=plot(Sp/koe,e,'-b');
%set(p,'LineWidth',3);
%hold on; grid on;
%xlabel({'Длина волны, мкм'});
%ylabel({'Сила излучения, Вт/(м^2*мкм)'});
%title({'График зависимости спектральной силы излучения от длины волны'});
frwl=dlr;
end