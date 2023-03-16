te0=273.15;
temvh=arrTemHigh(); 
temvh=temvh+te0; 
temvc=arrTemCold(); 
temvc=temvc+te0; 
tepo=arrTepPot84(); 
qv=(1/13.85e-4)*tepo; 
dop=[0.223543486340722
         0.408363698723139
         0.158234992027462
         0.412127160055975
         0.203844587055728
         0.360912423331915
         0.192931527902446
         0.359313615049683
         0.111783569987879
         0.238073543308064
         0.111301844365152
         0.253030375559658
         0.0946023670975657
         0.209362879056982
         0.171677558280794
         0.394798682703568
         0.206320819761141
         0.439209997762441
         0.104426255919279
         0.245152276440127
         0.111355275987299
         0.185364175025205
         0.235996993861958
         0.484263395535222];
disp(length(dop));
q=1; r=1; ktko=0;
tevho=0; temco=0;
temvs=(temvh+temvc)/2;
for k=1:length(dop)
    if (rem(k,2)==0)
        ktkoh(q)=dop(k);
        tevho(q)=temvs(k);
        q=q+1;
    else
        ktkoc(r)=dop(k);
        tevco(r)=temvs(k);
        r=r+1;
    end
end
p=length(tevho); tem1=0; tem2=0; tepv1=0; tepv2=0;
for k=1:p
tem1=tem1+tevho(k); tepv1=tepv1+ktkoh(k)*tevho(k);
end
p=length(tevco); 
for k=1:p
tem2=tem2+tevco(k); tepv2=tepv2+ktkoc(k)*tevco(k);
end
tepv1=tepv1/tem1
tepv2=tepv2/tem2
tem1=tem1/length(tevho)-te0
tem2=tem2/length(tevco)-te0
pl=plot([tem2 tem1],[tepv2 tepv1],'-b');
set(pl,'LineWidth',2); 
hold on; grid on; 
xlabel({'Температура, °С'}); 
ylabel({'Коэффициент кондуктивной теплопроводности, Вт/(м*К)'}); 
title({'График зависимости кондуктивного КТП от температуры'});