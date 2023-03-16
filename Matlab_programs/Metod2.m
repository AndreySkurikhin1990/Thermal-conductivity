l=5e-2;
lv=5e-4;
t1=65+273.15;
t2=445+273.15;
eps=(0.6e-8)*4184/3600; 
sig=5.668e-8;
epsil=eps/sig;
temv = [20 50 100 250 500];
temv = temv+273.15;
tepv = [0.045 0.047 0.051 0.065 0.096];
tepv=tepv*4184/3600;
kotepv = koef(temv,tepv,length(tepv))';
tv=350:50:550;
tv=tv+273.15;
tepvo=[0.0491 0.0521 0.0548 0.0574 0.0598];
kotepvo=koef(tv,tepvo,length(tepvo))';
for k=1:5 
    polj(k)=kotepv(k)/(6-k);
end
polj(6)=0;
jt1=polyval(polj,t1)/l;
jv=polyval(polj,t2)/l-jt1;
tepitom=0.23611111111;
vf=0.2;
n1=round((l/lv)*(vf^(1/3)));
lav=(l-lv*n1)/(n1-1);
maxn=149;
lamb=0;
nsl=0;
for dgv=1:4 
    m=1;
    for kl=3:2:maxn
    lamb(dgv,m)=jv*lv/RasTem(kl, dgv, kotepvo, jv, t2, l, lv, eps, sig);
    nsl(dgv,m)=kl;
    m=m+1;
    end
end
p=plot(nsl(1,:),lamb(1,:),'-m',nsl(2,:),lamb(2,:),'-r',nsl(3,:),lamb(3,:),'-c',nsl(4,:),lamb(4,:),'-b');
set(p,'LineWidth',3);
hold on;
grid on;
xlabel({'Число слоев'});
ylabel({'Коэффициент теплопроводности, Вт/(м*К)'});
title({'Зависимость коэффициента теплопроводности от количества слоев'});