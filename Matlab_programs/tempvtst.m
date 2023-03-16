function tem = tempvtst(T0,Th,len,N)
temv = [20 50 100 250 500]; 
te0=273.15;
temv = temv+te0; 
tepv = [0.045 0.047 0.051 0.065 0.096];
tepv = tepv*4184/3600; 
kotepv = koef(temv,tepv,length(tepv))'; 
Tz=(T0+Th)/2;
koeftep = abs(polyval(kotepv,Tz));
qv=abs(koeftep*(T0-Th)/(len*N));
Tz=T0;
for k=2:N
    koeftep = abs(polyval(kotepv,Tz));    
    Tz=Tz-qv*len/koeftep;
    tem(k)=Tz;
end
tem(1)=T0;
tem(N+1)=Th;
end