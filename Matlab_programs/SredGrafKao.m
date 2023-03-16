function [ alsr ] = SredGrafKao()
format long g;
dlvo=dlvokao1(); 
dlvo2=dlvokao2(); 
Tal=Trkokaol1(); 
Tal2=Trkokaol2(); 
ndlvo=length(dlvo);
mkbr=234.132; 
mkao=0.928; 
tol=0.68; 
rokbr=2.75; 
rkao=2.6; 
mkbr=mkbr-mkao;
mkbr2=234.323; 
mkao2=1.04; 
tol2=0.65; 
rkao2=rkao; 
mkbr2=mkbr2-mkao2;
vkbr=mkbr/(1e3*rokbr);
vkao=mkao/(1e3*rkao);
vkbr2=mkbr2/(1e3*rokbr);
vkao2=mkao2/(1e3*rkao2);
xsh=(vkao/(vkao+vkbr))*tol*1e3;
xsh2=(vkao2/(vkao2+vkbr2))*tol2*1e3;
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
for k=1:ndlvo
    alsre(k)=(Tal(k)+Tal2(k))/2;
end
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