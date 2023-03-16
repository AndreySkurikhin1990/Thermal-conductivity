%расчет конвективной составл€ющей при линеаризации пол€ температур
format long g;
Tna=1e3;
Tko=235;
lamb=0.245;
y0=30e-3;
Nt=3e3;
h=y0/Nt;
x=0:h:y0;
konak=(Tko-Tna)/y0;
for k=1:length(x)
    tem(k)=konak*x(k)+Tna;
end
lambda=rasKonvSostav(tem,lamb,Tna,Tko,x)