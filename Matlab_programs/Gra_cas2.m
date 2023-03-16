dlvo=[3000 2969.6 2947.8 2908.7 2847.8 2826.1 2782.6 2721.7 2691.3 2639.1 2500 2442.5 2436.7 2419.5 2400 2382.2 2322.9 2293.1 2250 2229.9 2206.9 2175.3 2125 2065.7 2062 2031.6 2005.6 2000 1972.4 1958.6 1936.2 1922.4 1905.2 1889.7 1862];
Tal=[40 50.8 60 64.6 66.7 71.5 74.4 76.8 80 81.2 81 80.8 80 79 77.4 78.3 79.4 77.4 72.9 75.5 65.7 60 49.5 44.5 44.2 47.4 48.6 48.2 40 30.1 20 9.3 7.5 5.2 4.6];
ndlvo=length(dlvo);
dob=201;
for k=1:ndlvo 
    dlvo(k)=(0.01/dlvo(k))*1e6;
    Tal(k)=-log(0.01*Tal(k))/dob; 
end
ns=1.616;
lamb=0.24;
pi=3.1415926535897932;
sig=5.668e-8;
tC1=550;
tC2=500;
teta(1)=tC1+273.15;
teta(2)=tC2+273.15;
ndlvo=length(dlvo);
tetas=(teta(1)+teta(2))/2;
lam=2898/tetas;
for q=2:ndlvo-1 
if dlvo(q)>lam break;
end
end
ko=0; la=0; hi=0; te=0; koR=0; koe=0;
la(1)=dlvo(q-1); la(2)=dlvo(q); al(1)=Tal(q-1); al(2)=Tal(q); 
koe=(al(2)-al(1))/(la(2)-la(1));
alp=koe*(lam-la(1))+al(1);
hi=(alp*lam)/(4*pi);
ro=((ns-1)^2+hi^2)/((ns+1)^2+hi^2);
alp=alp*1e6;
m=(2*alp*sig*(teta(1)^4-teta(2)^4)/(lamb*(teta(1)-teta(2)))+alp^2)^0.5;
d=1e-3;
H1=2*sig*(teta(1)^4-teta(2)^4);
H1=H1+alp*lamb*(teta(1)-teta(2));
H2=m*lamb*(teta(1)-teta(2))*(1+ro)*(1+exp(-m*d));
H2=H2/((1-ro)*(1-exp(-m*d)));
H1=H1+H2;
H2=2*(1+ro)/(1-ro);
H2=H2+alp*d;
H3=2*(1+ro)*alp^2;
H3=H3/((m^2)*(1-ro));
H2=H2-H3;
H3=d*(1+ro)*(1+exp(-m*d))*alp^2;
H3=H3/(m*(1-ro)*(1-exp(-m*d)));
H2=H2+H3;
H=H1/H2;
A=alp*lamb*teta(1)+2*sig*teta(1)^4+0.5*alp*H*d-(lamb*(teta(1)-teta(2))*m^2)/(2*alp);
B=0.5*(teta(1)-teta(2)-(H*d*alp^2)/(lamb*m^2))/(1-exp(-m*d));
Nx=1e4;
tetax=0;
koo=0:d/Nx:d;
for k=1:length(koo)
tetax(k)=B*(exp(-m*koo(k))-exp(-m*(d-koo(k))))-(koo(k)*H*alp^2-alp*A-6*alp*sig*tetas^4)/(lamb*m^2);
end
koo=koo*1e6;
tetax=tetax-273.15;
p=plot(koo,tetax,'-b');
set(p,'LineWidth',3);
hold on;
grid on;
xlabel({'Координата, мкм'});
ylabel({'Температура, град С'});
title({'График зависимости температуры от координаты'});