function [ alsr ] = SredGrafItom()
format long g;
dlvo=dlvoItom1(); 
dlvo2=dlvoItom2(); 
Tal=TrkoItom1(); 
Tal2=TrkoItom2(); 
ndlvo=length(dlvo);
mkbr=230.078; 
mitom=0.952; 
tol=0.64; 
rokbr=2.75; 
ritom=0.44; 
mkbr=mkbr-mitom;
mkbr2=231.006; 
mitom2=1.189; 
tol2=0.64; 
ritom2=ritom; 
mkbr2=mkbr2-mitom2;
vkbr=mkbr/(1e3*rokbr);
vitom=mitom/(1e3*ritom);
vkbr2=mkbr2/(1e3*rokbr);
vitom2=mitom2/(1e3*ritom2);
xsh=(vitom/(vitom+vkbr))*tol*1e3;
xsh2=(vitom2/(vitom2+vkbr2))*tol2*1e3;
for k=1:ndlvo      
    dlvo(k)=(1e4/dlvo(k));         
    Tal(k)=-log(Tal(k))/xsh; 
end
ndlvo2=length(dlvo2);
for k=1:ndlvo2 
    dlvo2(k)=(1e4/dlvo2(k));     
    Tal2(k)=-log(Tal2(k))/xsh2;    
end
%alsre=postrgr(dlvo,dlvo2,Tal,Tal2);
for k=1:(ndlvo+ndlvo2)/2
    alsre(k)=(Tal(k)+Tal2(k))/2;
end
alsre=1e6*alsre;
k=ZapisFile(alsre);
alsr=alsre;
end

function pg = postrgr(x1,x2,y1,y2)
p=plot(x1,y1,'-k',x2,y2,'-b');
set(p,'LineWidth',2); 
hold on; grid on; 
xlabel({'Длина волны, мкм'}); 
ylabel({'Коэффициент поглощения, 1/мкм'}); 
title({'График зависимости КП от длины волны'});
legend('Образец 1', 'Образец 2','location','best');
pg=0;
end

function t = ZapisFile(massi)
fid = fopen('Koefficient_pogloscheniya_itom.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.12f\n',massi(k));
end
fclose(fid);
t=0;
end