function alsre = InteDlVoShamot;
%dlvo=[3000 2969.6 2947.8 2908.7 2847.8 2826.1 2782.6 2721.7 2691.3 2639.1 2500 2442.5 2436.7 2419.5 2400 2382.2 2322.9 2293.1 2250 2229.9 2206.9 2175.3 2125 2065.7 2062 2031.6 2005.6 2000 1972.4 1958.6 1936.2 1922.4 1905.2 1889.7 1862];
%Tal=[40 50.8 60 64.6 66.7 71.5 74.4 76.8 80 81.2 81 80.8 80 79 77.4 78.3 79.4 77.4 72.9 75.5 65.7 60 49.5 44.5 44.2 47.4 48.6 48.2 40 30.1 20 9.3 7.5 5.2 4.6];
%dob=201*(1e-6); 
dlvo=dlvoSham1(); dlvo2=dlvoSham2(); 
Tal=TrkoSham1(); Tal2=TrkoSham2(); 
ndlvo=length(dlvo);  ndlvo2=length(dlvo2); t0=273.15;
mkbr=238.57; msh=1.706; tol=0.7; rokbr=2.75; rsh=2.7;
mkbr2=227.973; msh2=1.1; tol2=0.68; rsh2=rsh;
vkbr=mkbr/(1e3*rokbr); 
vsh=msh/(1e3*rsh); 
vkbr2=mkbr2/(1e3*rokbr); 
vsh2=msh2/(1e3*rsh2);
xsh=(vsh/(vsh+vkbr))*tol*1e3; 
xsh2=(vsh2/(vsh2+vkbr2))*tol2*1e3;
for k=1:ndlvo      
    dlvo(k)=(1e4/dlvo(k));         
    Tal(k)=-log(Tal(k))/xsh; 
end 
for k=1:ndlvo2     
    dlvo2(k)=(1e4/dlvo2(k));     
    Tal2(k)=-log(Tal2(k))/xsh2;    
end
for k=1:ndlvo 
    te1(k)=2898/dlvo(k)-t0; 
end
for k=1:ndlvo2 
    te2(k)=2898/dlvo2(k)-t0;
end
j=1; 
for k=1:ndlvo 
    if ((te1(k)>0) && (te1(k)<1000) )
            tem1(j)=te1(k); 
            al1(j)=Tal(k); 
            dlv1(j)=dlvo(k); 
            j=j+1;  
    end
end
j=1; 
for k=1:ndlvo 
    if ((te2(k)>0) && (te2(k)<1e3)) 
        tem2(j)=te2(k);
        al2(j)=Tal2(k);
        dlv2(j)=dlvo(k); 
        j=j+1;  
    end 
end
als=trapz(dlvo,Tal)/(dlvo(ndlvo)-dlvo(1)); 
als2=trapz(dlvo2,Tal2)/(dlvo2(ndlvo2)-dlvo2(1)); 
alsr(1)=trapz(dlv1,al1)/(dlv1(length(dlv1))-dlv1(1)); 
alsr(2)=trapz(dlv2,al2)/(dlv2(length(dlv2))-dlv2(1));
p=plot(dlvo,Tal,'-b',dlvo2,Tal2,'-g');
set(p,'LineWidth',2); hold on; grid on;
xlabel({'Äëèíà âîëíû, ìêì'}); 
ylabel({'Êîıôôèöèåíò ïîãëîùåíèÿ, ìêì-1'}); 
title({'Ãğàôèê çàâèñèìîñòè êîıôôèöèåíòà ïîãëîùåíèÿ îò äëèíû âîëíû'});
alsre=(alsr(1)+alsr(2))/2;
IntegrDlVol=alsre;
end