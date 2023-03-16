function [ alsre ] = IntegrDlVolVer55;
dlvo=dlvoVer5551(); dlvo2=dlvoVer5552(); 
Tal=TrkoVer5551(); Tal2=TrkoVer5552(); 
ndlvo=length(dlvo);  ndlvo2=length(dlvo2); 
mkbr=249.913; mv=0.315; tol=0.74; rokbr=2.75; rov=0.4; 
mkbr2=249.607; mv2=0.473; tol2=0.74; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);

dlvo3=dlvoVer5553(); dlvo4=dlvoVer5561(); 
Tal3=TrkoVer5553(); Tal4=TrkoVer5561(); 
ndlvo3=length(dlvo3);  ndlvo4=length(dlvo4); 
mkbr3=249.218; mv3=0.709; tol3=0.72; 
mkbr4=249.929; mv4=0.293; tol4=0.72; rov3=0.4; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;

dlvo5=dlvoVer5562(); dlvo6=dlvoVer5563(); 
Tal5=TrkoVer5562(); Tal6=TrkoVer5563(); 
ndlvo5=length(dlvo5);  ndlvo6=length(dlvo6); 
mkbr5=249.695; mv5=0.528; tol5=0.71; rokbr=2.75; rov5=0.4; 
mkbr6=249.306; mv6=0.83; tol6=0.7; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);

dlvo7=dlvoVer5571(); dlvo8=dlvoVer5572(); 
Tal7=TrkoVer5571(); Tal8=TrkoVer5572(); 
ndlvo7=length(dlvo7);  ndlvo8=length(dlvo8); 
mkbr7=250.405; mv7=0.27; tol7=0.78; rov7=0.4;
mkbr8=249.625; mv8=0.493; tol8=0.73; rov8=rov7;
vkbr7=mkbr7/(1e3*rokbr); vv7=mv7/(1e3*rov7); vkbr8=mkbr8/(1e3*rokbr); vv8=mv4/(1e3*rov8);

dlvo9=dlvoVer5573(); Tal9=TrkoVer5573(); ndlvo9=length(dlvo9); 
mkbr9=249.348; mv9=0.764; tol9=0.76; rov9=0.4;
vkbr9=mkbr9/(1e3*rokbr); vv9=mv9/(1e3*rov9); vkbr9=mkbr9/(1e3*rokbr); vv9=mv9/(1e3*rov9);

xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; xv7=(vv7/(vv7+vkbr7))*tol7*1e3; xv8=(vv8/(vv8+vkbr8))*tol8*1e3;
xv9=(vv9/(vv9+vkbr9))*tol9*1e3;

for k=1:ndlvo      dlvo(k)=(1e4/dlvo(k)); Tal(k)=-log(Tal(k))/xv; end; 
for k=1:ndlvo2     dlvo2(k)=(1e4/dlvo2(k));     Tal2(k)=-log(Tal2(k))/xv2;    end;  
for k=1:ndlvo3     dlvo3(k)=(1e4/dlvo3(k));     Tal3(k)=-log(Tal3(k))/xv3; end; 
for k=1:ndlvo4     dlvo4(k)=(1e4/dlvo4(k));      Tal4(k)=-log(Tal4(k))/xv4; end;
for k=1:ndlvo5     dlvo5(k)=(1e4/dlvo5(k));      Tal5(k)=-log(Tal5(k))/xv5; end;
for k=1:ndlvo6     dlvo6(k)=(1e4/dlvo6(k));      Tal6(k)=-log(Tal6(k))/xv6; end;
for k=1:ndlvo7     dlvo7(k)=(1e4/dlvo7(k));      Tal7(k)=-log(Tal7(k))/xv7; end;
for k=1:ndlvo8     dlvo8(k)=(1e4/dlvo8(k));      Tal8(k)=-log(Tal8(k))/xv8; end;
for k=1:ndlvo9     dlvo9(k)=(1e4/dlvo9(k));      Tal9(k)=-log(Tal9(k))/xv9; end;

for k=1:ndlvo te1(k)=2898/dlvo(k)-t0; end; 
for k=1:ndlvo2 te2(k)=2898/dlvo2(k)-t0; end; 
for k=1:ndlvo3 te3(k)=2898/dlvo3(k)-t0; end; 
for k=1:ndlvo4 te4(k)=2898/dlvo4(k)-t0; end;
for k=1:ndlvo5 te5(k)=2898/dlvo5(k)-t0; end;
for k=1:ndlvo6 te6(k)=2898/dlvo6(k)-t0; end;
for k=1:ndlvo7 te7(k)=2898/dlvo7(k)-t0; end;
for k=1:ndlvo8 te8(k)=2898/dlvo8(k)-t0; end;
for k=1:ndlvo9 te9(k)=2898/dlvo9(k)-t0; end;

j=1; for k=1:ndlvo if (te1(k)>0) if (te1(k)<1200) tem1(j)=te1(k); al1(j)=Tal(k); dlv1(j)=dlvo(k); j=j+1;  end; end; end; 
j=1; for k=1:ndlvo2 if (te2(k)>0) if (te2(k)<1200) tem2(j)=te2(k); al2(j)=Tal2(k);  dlv2(j)=dlvo(k); j=j+1;  end; end; end;  
j=1; for k=1:ndlvo3 if (te3(k)>0) if (te3(k)<1200) tem3(j)=te3(k); al3(j)=Tal3(k);  dlv3(j)=dlvo3(k); j=j+1;  end; end; end;  
j=1; for k=1:ndlvo4 if (te4(k)>0) if (te4(k)<1200) tem4(j)=te4(k); al4(j)=Tal4(k);  dlv4(j)=dlvo4(k); j=j+1;  end; end; end; 
j=1; for k=1:ndlvo5 if (te5(k)>0) if (te5(k)<1200) tem5(j)=te5(k); al5(j)=Tal5(k);  dlv5(j)=dlvo5(k); j=j+1;  end; end; end; 
j=1; for k=1:ndlvo6 if (te6(k)>0) if (te6(k)<1200) tem6(j)=te6(k); al6(j)=Tal6(k);  dlv6(j)=dlvo6(k); j=j+1;  end; end; end; 
j=1; for k=1:ndlvo7 if (te7(k)>0) if (te7(k)<1200) tem7(j)=te7(k); al7(j)=Tal7(k);  dlv7(j)=dlvo7(k); j=j+1;  end; end; end; 
j=1; for k=1:ndlvo8 if (te8(k)>0) if (te8(k)<1200) tem8(j)=te8(k); al8(j)=Tal8(k);  dlv8(j)=dlvo8(k); j=j+1;  end; end; end; 
j=1; for k=1:ndlvo9 if (te9(k)>0) if (te9(k)<1200) tem9(j)=te9(k); al9(j)=Tal9(k);  dlv9(j)=dlvo9(k); j=j+1;  end; end; end; 

j=1; for k=1:ndlvo if (dlvo(k)>3) if (dlvo(k)<10) dll(j)=dlvo(k); all(j)=Tal(k); j=j+1; end; end; end;
j=1; for k=1:ndlvo2 if (dlvo2(k)>3) if (dlvo2(k)<10) dll2(j)=dlvo2(k); all2(j)=Tal2(k); j=j+1; end; end; end;
j=1; for k=1:ndlvo3 if (dlvo3(k)>3) if (dlvo3(k)<10) dll3(j)=dlvo3(k); all3(j)=Tal3(k); j=j+1; end; end; end;
j=1; for k=1:ndlvo4 if (dlvo4(k)>3) if (dlvo4(k)<10) dll4(j)=dlvo4(k); all4(j)=Tal4(k); j=j+1; end; end; end;
j=1; for k=1:ndlvo5 if (dlvo5(k)>3) if (dlvo5(k)<10) dll5(j)=dlvo5(k); all5(j)=Tal5(k); j=j+1; end; end; end;
j=1; for k=1:ndlvo6 if (dlvo6(k)>3) if (dlvo6(k)<10) dll6(j)=dlvo6(k); all6(j)=Tal6(k); j=j+1; end; end; end;
j=1; for k=1:ndlvo7 if (dlvo7(k)>3) if (dlvo7(k)<10) dll7(j)=dlvo7(k); all7(j)=Tal7(k); j=j+1; end; end; end;
j=1; for k=1:ndlvo8 if (dlvo8(k)>3) if (dlvo8(k)<10) dll8(j)=dlvo8(k); all8(j)=Tal8(k); j=j+1; end; end; end;
j=1; for k=1:ndlvo9 if (dlvo9(k)>3) if (dlvo9(k)<10) dll9(j)=dlvo9(k); all9(j)=Tal9(k); j=j+1; end; end; end;

ndlvo=length(dlvo); ndlvo2=length(dlvo2); ndlvo3=length(dlvo3); ndlvo4=length(dlvo4);
als=trapz(dlvo,Tal)/(dlvo(ndlvo)-dlvo(1))
als2=trapz(dlvo2,Tal2)/(dlvo2(ndlvo2)-dlvo2(1))
als3=trapz(dlvo3,Tal3)/(dlvo3(ndlvo3)-dlvo3(1))
als4=trapz(dlvo4,Tal4)/(dlvo4(ndlvo4)-dlvo4(1))
ndlvo5=length(dlvo5); ndlvo6=length(dlvo6); ndlvo7=length(dlvo7); ndlvo8=length(dlvo8); ndlvo9=length(dlvo9);
als5=trapz(dlvo5,Tal5)/(dlvo5(ndlvo5)-dlvo5(1))
als6=trapz(dlvo6,Tal6)/(dlvo6(ndlvo6)-dlvo6(1)) 
als7=trapz(dlvo7,Tal7)/(dlvo7(ndlvo7)-dlvo7(1))
als8=trapz(dlvo8,Tal8)/(dlvo8(ndlvo8)-dlvo8(1))
als9=trapz(dlvo9,Tal9)/(dlvo9(ndlvo9)-dlvo9(1))

ndlvo=length(dlv1); ndlvo2=length(dlv2); ndlvo3=length(dlv3); ndlvo4=length(dlv4);
ndlvo5=length(dlv5); ndlvo6=length(dlv6); ndlvo7=length(dlv7); ndlvo8=length(dlv8); ndlvo9=length(dlv9);

alsr(1)=trapz(dlv1,al1)/(dlv1(length(dlv1))-dlv1(1)); alsr(2)=trapz(dlv2,al2)/(dlv2(length(dlv2))-dlv2(1)); 
alsr(3)=trapz(dlv3,al3)/(dlv3(length(dlv3))-dlv3(1)); alsr(4)=trapz(dlv4,al4)/(dlv4(length(dlv4))-dlv4(1));
alsr(5)=trapz(dlv5,al5)/(dlv5(length(dlv5))-dlv5(1)); alsr(6)=trapz(dlv6,al6)/(dlv6(length(dlv6))-dlv6(1));
alsr(7)=trapz(dlv7,al7)/(dlv7(length(dlv7))-dlv7(1)); alsr(8)=trapz(dlv8,al8)/(dlv8(length(dlv8))-dlv8(1));
alsr(9)=trapz(dlv9,al9)/(dlv9(length(dlv9))-dlv9(1));

alsr(10)=dlv1(1); alsr(11)=dlv1(length(dlv1)); alsr(12)=dlv2(1); alsr(13)=dlv2(length(dlv2)); 
alsr(14)=dlv3(1); alsr(15)=dlv3(length(dlv3)); alsr(16)=dlv4(1); alsr(17)=dlv4(length(dlv4));
alsr(18)=dlv5(1); alsr(19)=dlv5(length(dlv5)); alsr(19)=dlv6(1); alsr(20)=dlv6(length(dlv6));
alsr(21)=dlv7(1); alsr(22)=dlv7(length(dlv7)); alsr(23)=dlv8(1); alsr(24)=dlv8(length(dlv8)); 
alsr(25)=dlv9(1); alsr(26)=dlv9(length(dlv9));

%p=plot(dlvo,Tal,'-b',dlvo2,Tal2,'-g',dlvo3,Tal3,'-m');set(p,'LineWidth',2);  hold on; grid on;
%legend(p,'0.13 %','0.19 %','0.28 %');  set(gca,'FontName','TimesET'); set(gca,'FontSize',10);
%xlabel({'Длина волны, мкм'}); ylabel({'Коэффициент поглощения, мкм-1'}); title({'График зависимости коэффициента поглощения от длины волны'}); figure;
%p1=plot(dlvo4,Tal4,'-b',dlvo5,Tal5,'-r',dlvo6,Tal6,'-c');set(p1,'LineWidth',2);  hold on; grid on;
%legend(p1,'0.12 %','0.21 %','0.33 %');  set(gca,'FontName','TimesET'); set(gca,'FontSize',10);
%xlabel({'Длина волны, мкм'}); ylabel({'Коэффициент поглощения, мкм-1'}); title({'График зависимости коэффициента поглощения от длины волны'}); figure;
%p2=plot(dlvo7,Tal7,'-m',dlvo8,Tal8,'-g',dlvo9,Tal9,'-k');
%set(p2,'LineWidth',2); hold on; grid on; xlabel({'Длина волны, мкм'}); ylabel({'Коэффициент поглощения, мкм-1'});  title({'График зависимости коэффициента поглощения от длины волны'});
%legend(p2,'0.11 %','0.2 %','0.31 %');  set(gca,'FontName','TimesET'); set(gca,'FontSize',10);

p=plot(dll,all,'-b',dll2,all2,'-g',dll3,all3,'-m'); set(p,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'}); ylabel({'Коэффициент поглощения, мкм-1'}); title({'График зависимости коэффициента поглощения от длины волны'});
legend(p,'0.13 %','0.19 %','0.28 %');  set(gca,'FontName','TimesET'); set(gca,'FontSize',10); figure;
p1=plot(dll4,all4,'-r',dll5,all5,'-c', dll6,all6,'-b'); set(p1,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'}); ylabel({'Коэффициент поглощения, мкм-1'}); title({'График зависимости коэффициента поглощения от длины волны'});
legend(p1,'0.12 %','0.21 %','0.33 %');  set(gca,'FontName','TimesET'); set(gca,'FontSize',10); figure;
p2=plot(dll7,all7,'-g',dll8,all8,'-m',dll9,all9,'-r');set(p2,'LineWidth',2); hold on; grid on;
xlabel({'Длина волны, мкм'}); ylabel({'Коэффициент поглощения, мкм-1'}); title({'График зависимости коэффициента поглощения от длины волны'});
legend(p2,'0.11 %', '0.2 %','0.31 %');  set(gca,'FontName','TimesET'); set(gca,'FontSize',10);

ndll=length(dll); ndll2=length(dll2); ndll3=length(dll3); ndll4=length(dll4); 
ndll5=length(dll5); ndll6=length(dll6); ndll7=length(dll7); ndll8=length(dll8); ndll9=length(dll9); 

als=trapz(dll,all)/(dll(ndll)-dll(1))
als2=trapz(dll2,all2)/(dll2(ndll2)-dll2(1))
als3=trapz(dll3,all3)/(dll3(ndll3)-dll3(1))
als4=trapz(dll4,all4)/(dll4(ndll4)-dll4(1))
als5=trapz(dll5,all5)/(dll5(ndll5)-dll5(1))
als6=trapz(dll6,all6)/(dll6(ndll6)-dll6(1))
als7=trapz(dll7,all7)/(dll7(ndll7)-dll7(1))
als8=trapz(dll8,all8)/(dll8(ndll8)-dll8(1))
als9=trapz(dll9,all9)/(dll9(ndll9)-dll9(1))

als=1/((als+als2+als4+als3+als5+als6+als7+als8+als9)/9)
alsre=alsr;
IntegrDlVol=alsre;
end