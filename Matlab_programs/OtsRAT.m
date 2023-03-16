dlvo=[3000 2969.6 2947.8 2908.7 2847.8 2826.1 2782.6 2721.7 2691.3 2639.1 2500 2442.5 2436.7 2419.5 2400 2382.2 2322.9 2293.1 2250 2229.9 2206.9 2175.3 2125 2065.7 2062 2031.6 2005.6 2000 1972.4 1958.6 1936.2 1922.4 1905.2 1889.7 1862];
Tal=[40 50.8 60 64.6 66.7 71.5 74.4 76.8 80 81.2 81 80.8 80 79 77.4 78.3 79.4 77.4 72.9 75.5 65.7 60 49.5 44.5 44.2 47.4 48.6 48.2 40 30.1 20 9.3 7.5 5.2 4.6];
ndlvo=length(dlvo);
dob=201;
for k=1:ndlvo 
    dlvo(k)=(0.01/dlvo(k))*1e6;
    Tal(k)=-log(0.01*Tal(k))/dob; 
end
temp=350:50:500;
temp=temp+273.15;
npp=1.616;
na=length(temp);
lam=0;
pi=3.1415926535897932;
lv=500;
m=1;
maxn=201;
koPo=0;
ts=0;
Tr=0;
for n=3:2:maxn
tol=lv/n;
alp=0;
    for k=1:na 
        lam(k)=2898/temp(k);
        for q=2:ndlvo-1 
            if dlvo(q)>lam(k) 
            break;
        end
        end
ko=0; la=0; hi=0; te=0; koR=0; koe=0;
la(1)=dlvo(q-1); la(2)=dlvo(q); al(1)=Tal(q-1); al(2)=Tal(q); 
koe=(al(2)-al(1))/(la(2)-la(1)); 
        alp=koe*(lam(k)-la(1))+al(1);
        te=exp(-alp*tol);
        koPoV(m,k)=(1-te)*100;
end
        ts(m)=tol;
        m=m+1;
end
temp=0;
temp=300:5:500;
temp=temp+273.15;
na=length(temp);
alp=0;
for k=1:na 
        lam(k)=2898/temp(k);
        for q=2:ndlvo-1 
            if dlvo(q)>lam(k) 
            break;
        end
        end
ko=0; la=0; hi=0; te=0; koR=0; koe=0;
la(1)=dlvo(q-1); la(2)=dlvo(q); al(1)=Tal(q-1); al(2)=Tal(q); 
koe=(al(2)-al(1))/(la(2)-la(1)); 
        alp=koe*(lam(k)-la(1))+al(1);
        hi=(alp*lam(k))/(4*pi);
        koR=((npp-1)^2+hi^2)/((npp+1)^2+hi^2);
        Tr(k)=100*(1-koR);
end
subplot(2,1,1);
p1=plot(ts,koPoV(:,1)','-m',ts,koPoV(:,2)','-r',ts,koPoV(:,3)','-c',ts,koPoV(:,4)','-b');
set(p1,'LineWidth',2);
grid on; 
xlabel({'Толщина слоя, мкм'}); 
ylabel({'Коэффициент поглощения, %'}); 
title({'Зависимость коэффициента поглощения в объеме от толщины слоя'});
subplot(2,1,2);
temp=temp-273.15;
p2=plot(temp,Tr,'-b');
set(p2,'LineWidth',2);
grid on; 
xlabel({'Температура, град С'}); 
ylabel({'Коэффициент прохождения, %'}); 
title({'Зависимость коэффициента прохождения на границе от температуры'});