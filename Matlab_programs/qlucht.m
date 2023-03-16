format long g; teta=0; te0=273.15; T0=1000+te0; tol=3e4; 
y0=1e3; l=2e1; Ntt=ceil(tol/l); npp=Kramers_n(); dl=dlvoVer53101; p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k); 
end;
temvh=[487 870 585 1000 603 1000 603 1000 1000];
temvc=[214 435 109 235 108.8 238 101 199 228];
tepo=[3.0671 3.9805 3.9596 8.6377 3.5615 7.123 2.3003 5.3674 8.1705];
qv = (1e4/13.85)*tepo;
temvh=temvh+te0;
temvc=temvc+te0;
r=length(temvh);
for k=1:r
    kotep(k)=qv*tol*1e-6/(temvc(k)-temvh(k));
end
temvs=(1/2)*(temvc+temvh);
kktp1=(tepv(2)-tepv(1))/(temvs(2)-temvs(1));
kktp2=tepv(2)-kktp1*temvs(2);
xa = 0:l:tol; 
teta = nachPribTem(T0,l,19e-4/3600,tol,y0)
Tzh=temvh(2); 
Tzc=temvh(1);
p=length(xa);
for k=1:p-1
    koeftep = abs(kktp1*Tzh+kktp2); 
    laob(k)=koeftep;
    teta(k)=Tzc;
    Tzh=Tzh-qv(2)*abs(xk-xn)*1e-6/koefteph; 
    Tzc=Tzc-qv(1)*abs(xk-xn)*1e-6/koeftepc;
end
tetah(1)=temvh(2); 
tetah(p)=temvc(2); 
tetac(1)=temvh(1); 
tetac(p)=temvc(1); 
sigma=5.668e-8;
qx1=0; 
qx2=0; 
laizc=0;
laizh=0;
alfc=1/2e2; alfh=1/2e2; npc=1.6; nph=1.6;
for k=2:p-1
    xn=xa(k-1)*1e-6; 
    xk=xa(k)*1e-6;
    grc=(tetac(k)-tetac(k-1))/abs(xk-xn);
    grh=(tetah(k)-tetah(k-1))/abs(xk-xn);
    Tzc=tetac(k);
    Tzh=tetah(k);
    %alfc=1/sprko(Tzc);
    %alfh=1/sprko(Tzh);
    %npc=real(nsreddv(dl,npp,Tzc));
    %nph=real(nsreddv(dl,npp,Tzh));
    t=abs(16*(npc^2)*sigma*(Tzc^3)/(3*alfc*1e6));
    laizc(k)=t;
    t=t*grc;
    qx1(k)=t;
    t=abs(16*(nph^2)*sigma*(Tzh^3)/(3*alfh*1e6));
    laizh(k)=t;
    t=t*grh;
    qx2(k)=t;
    %disp(k);
end
tetac=tetac-te0;
tetah=tetah-te0;
disp(tetac);
disp(tetah);
%for k=1:length(qx) if (rem(k-2,2)==0) disp(qx(k)); end; end;
%for k=1:length(qx) if (rem(k-2,2)==0) disp(xa(k)*1e-6); end; end;
p=plot(tetac(2:p-1),laizc(2:p-1),'-b',tetah(2:p-1),laizh(2:p-1),'-k',tetac(2:p-1),laobc(2:p-1),'-m',tetah(2:p-1),laobh(2:p-1),'-c');
%p=plot(tetac(2:p-1),tetah(2:p-1),'-k');
set(p,'LineWidth',1); hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Плотность излучательного теплового потока, кВт/м2'}); 
%title({'График зависимости составляющих коэффициентов теплопроводности от координаты x'});
%ylabel({'Коэффициент одной из теплопроводностей, Вт/(м*К)'}); 
%legend('Лучистая 356 °C','Лучистая 618 °C', 'Общая 356 °C','Общая 618 °C');