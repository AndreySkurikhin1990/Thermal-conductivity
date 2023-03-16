%определяет КП вермикулита (ранние измерения)
dlv=dlvoVer5311();
al1=AbsKoVer5311();
al2=AbsKoVer5312();
al3=AbsKoVer5313();
tr1=TrkoVer5311();
tr2=TrkoVer5312();
tr3=TrkoVer5313();
p=length(al1);
mkbr=250.239; mv=0.283; tol=0.73; rokbr=2.75; rov=0.4; 
mkbr2=249.740; mv2=0.464; tol2=0.72; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); 
vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
mkbr3=249.223; mv3=0.812; tol3=0.72; rov3=0.4; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); 
xv=(vv/(vv+vkbr))*tol*1e3; 
xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3;
t1=0; t2=0; t3=0;
for k=1:p
    t1(k)=1-al1(k)-tr1(k);
    t2(k)=1-al2(k)-tr2(k);
    t3(k)=1-al3(k)-tr3(k);
end
alp1=0; alp2=0; alp3=0;
for k=1:p
    dlv(k)=1e4/dlv(k);
    alp1(k)=-log(1-al1(k))/xv;
    alp2(k)=-log(1-al2(k))/xv2;
    alp3(k)=-log(1-al3(k))/xv3;
end;
for k=1:p
    tr1(k)=-log(tr1(k))/xv;
    tr2(k)=-log(tr2(k))/xv2;
    tr3(k)=-log(tr3(k))/xv3;
end;
disp(t2);
p=plot(dlv,tr1,'-b',dlv,tr2,'-k',dlv,tr3,'-r',dlv,alp1,'-m',dlv,alp2,'-g',dlv,alp3,'-c');
set(p,'LineWidth',1); 
hold on; 
grid on; 
xlabel({'Длина волны, мкм'}); 
ylabel({'Коэффициент поглощения, мкм^-1'}); 
title({'График зависимости температуры от координаты'});