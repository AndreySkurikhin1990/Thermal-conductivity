function [lam] = tmp35()
%рассчитывает теплопроводность плоскопараллельной засыпки при любых
%температурах путем аппроксимации из данных предыдущих измерений
%строит график изменения общей плотности теплового потока от координаты
format long g; y0=(3e4)*(1e-6); te0=273.15; ep=1e-20; 
tem0=0; tem0=arrTem_VVF2; tem0=tem0+te0; 
ktp0=0; ktp0=arrKTP_VVF2; ketp0=polyfit(tem0,ktp0,2);
%temvh=0; temvh=arrTemHigh207(); temvh=temvh+te0;
%temvc=0; temvc=arrTemCold207(); temvc=temvc+te0;
%tepo=0; tepo=arrTepPot207(); qv=0; qv=(1e4/13.85)*tepo; 
%pt=length(temvc); ts=0; ts=sredtemp(temvh,temvc,pt); 
%tepv=0; tepv=koeftepv(qv,y0,temvh,temvc,pt);
%temvs=(1/2)*(temvc+temvh);
%kc=0; kc=danPoTemC2071(temvs,temvc);
%kh=0; kh=danPoTemH2071(temvs,temvh);
%kt=0; kt=danPoTemTepl2071(temvs,tepv);
tepv=arrKTP_VVF1();
tepvs1=(tepv(1)+tepv(2)+tepv(3))/3
tepvs2=tepv(4)
tepvs3=(tepv(5)+tepv(6))/2
temvh=0; temvh=arrTem1VVF1(); temvh=temvh+te0;
temvc=0; temvc=arrTem2VVF1(); temvc=temvc+te0;
temv3=0; temv3=arrTem3VVF1(); temv3=temv3+te0;
temvh1=(temvh(1)+temvh(2)+temvh(3))/3
temvc1=(temvc(1)+temvc(2)+temvc(3))/3
temv31=(temv3(1)+temv3(2)+temv3(3))/3
temvh2=temvh(4)
temvc2=temvc(4)
temv32=temv3(4)
temvh3=(temvh(5)+temvh(6))/2
temvc3=(temvc(5)+temvc(6))/2
temv33=(temv3(5)+temv3(6))/2
temvs=(temvc+temvh)/2;
temvss1=(temvs(1)+temvs(2)+temvs(3))/3;
temvss2=temvs(4);
temvss3=(temvs(5)+temvs(6))/2;
qob1=tepvs1*(temvh1-temvc1)/y0;
qob2=tepvs2*(temvh2-temvc2)/y0;
qob3=tepvs3*(temvh3-temvc3)/y0;
qob=[qob1,qob2,qob3]
tem1=[temvss1,temvss2,temvss3]
teps1=[tepvs1,tepvs2,tepvs3]
ketp1=polyfit(tem1,teps1,2);
ktexi600=polyfit([0,15,30]*1e-3,[temvh1,temv31,temvc1],2);
ktexi800=polyfit([0,15,30]*1e-3,[temvh2,temv32,temvc2],2);
ktexi1000=polyfit([0,15,30]*1e-3,[temvh3,temv33,temvc3],2);
xi=0:1e-2:30; xi=xi*1e-3; q=1; xk=0;
for k=1:length(xi)
    texi1000(k)=polyval(ktexi1000,xi(k));
    ktpxi1000(k)=polyval(ketp1,texi1000(k));
    texi800(k)=polyval(ktexi800,xi(k));
    ktpxi800(k)=polyval(ketp1,texi800(k));
    texi600(k)=polyval(ktexi600,xi(k));
    ktpxi600(k)=polyval(ketp1,texi600(k));
    if (k>1)
        dtexi=(texi1000(k)-texi1000(k-1))/(xi(k)-xi(k-1));
        qobxi1000(q)=abs(ktpxi1000(k)*dtexi);
        xk(q)=xi(k); tk1000(q)=texi1000(k);
        dtexi=(texi800(k)-texi800(k-1))/(xi(k)-xi(k-1));
        qobxi800(q)=abs(ktpxi800(k)*dtexi);
        tk800(q)=texi800(k);
        dtexi=(texi600(k)-texi600(k-1))/(xi(k)-xi(k-1));
        qobxi600(q)=abs(ktpxi600(k)*dtexi);
        xk(q)=xi(k); tk600(q)=texi600(k);
        ktpxim600(q)=ktpxi600(k); 
        ktpxim800(q)=ktpxi800(k);
        ktpxim1000(q)=ktpxi1000(k);
        q=q+1;
    end
end
ktps600=0; ktps800=0; ktps1000=0; n=length(ktpxi600);
for k=1:n
ktps600=ktps600+ktpxi600(k);
ktps800=ktps800+ktpxi800(k);
ktps1000=ktps1000+ktpxi1000(k);
end
ktps600=ktps600/n; ktps800=ktps800/n; ktps1000=ktps1000/n;
ktps=[ktps600,ktps800,ktps1000]
kqobxi1000=polyfit(tk1000,qobxi1000,2);
qobxi1000=qobxi1000';
qobxi800=qobxi800';
qobxi600=qobxi600';
tk1000=tk1000';
tk800=tk800';
tk600=tk600';
%pl=plot(tk1000,qobxi1000,'-b',tk800,qobxi800,'-k',tk600,qobxi600,'-r');
%pl=plot(xk,tk1000,'-b',xk,tk800,'-k',xk,tk600,'-r');
pl=plot(tk1000,ktpxim1000,'-b',tk800,ktpxim800,'-b',tk600,ktpxim600,'-b');
set(pl,'LineWidth',3); hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Плотность теплового потока, Вт/м2'}); 
%title({'График зависимости ПТП от температуры'});
%xlabel({'Координата, м'}); 
%ylabel({'Температура, К'}); 
%title({'График зависимости температуры от координаты'});
%legend('907 К','773 К','633 К','Location','Best');
xlabel({'Температура, К'}); 
ylabel({'КТП, Вт/(м*К)'}); 
title({'График зависимости КТП от температуры'});
dite=0; dite=3e1:1:1e3; dite=dite+te0; le=length(dite); proc=0;
for k=1:le 
    ktp1(k)=polyval(ketp1,dite(k));
    ktp2(k)=polyval(ketp0,dite(k));
    proc(k)=(ktp1(k)-ktp2(k))/ktp1(k);
end
temvh=0; temvh=arrTemHigh(); temvh=temvh+te0; 
temvc=0; temvc=arrTemCold(); temvc=temvc+te0; 
tepo=0; tepo=arrTepPot84(); qv=0; qv=(1e4/13.85)*tepo; 
pt=length(temvc); ts=0; ts=sredtemp(temvh,temvc,pt); 
tepv=0; tepv=koeftepv(qv,y0,temvh,temvc,pt);
temvs=(1/2)*(temvc+temvh);
kt=0; kt=danPoTemTepl(temvs,tepv);
kc=0; kc=danPoTemC(temvs,temvc);  
kh=0; kh=danPoTemH(temvs,temvh); 
temvhs1=(temvh(16)+temvh(18))/2;
temvcs1=(temvc(16)+temvc(18))/2;
temvss1=(temvhs1+temvcs1)/2
tepls1=(tepv(16)+tepv(18))/2
qob=tepls1*(temvhs1-temvcs1)/y0
temvhs2=(temvh(15)+temvh(17))/2;
temvcs2=(temvc(15)+temvc(17))/2;
temvss2=(temvhs2+temvcs2)/2
tepls2=(tepv(15)+tepv(17))/2
tep2=tepls2*(temvhs2-temvcs2)/y0
koe1=(tepls2-tepls1)/(temvss2-temvss1);
koe2=tepls2-temvss2*koe1;
tepl=(temvc(15)+temvc(17))/2; tepl=(temvc(16)+temvc(18))/2;
tepl=(temvh(15)+temvh(17))/2; tepl=(temvh(16)+temvh(18))/2;
tepl=(tepv(15)+tepv(17))/2; tepl=(tepv(16)+tepv(18))/2;
tem3=0; tem3=[301,367,466,615,765,914,1064,1264]; 
ktp3=0; ktp4=0; kt=0; kt=[koe1,koe2]; proc3=0; %kproc=polyfit(dite,proc,2);
for k=1:length(tem3) 
    ktp3(k)=polyval(kt,tem3(k));
    %ktp4(k)=(1-polyval(kproc,tem3(k)))*ktp3(k);
    proc3(k)=fitti(dite,proc,tem3(k));
    ktp4(k)=(1-proc3(k))*ktp3(k);
end
tem3=tem3';
ktp3=ktp3';
ktp4=ktp4';
tesr=PoKo(ktp3,tem3,ktp2,tem3,temvss2,temvss1);
temvhs3=temvhs2+(temvhs1-temvhs2)*(tesr-temvss2)/(temvss1-temvss2);
temvcs3=temvcs2+(temvcs1-temvcs2)*(tesr-temvss2)/(temvss1-temvss2);
qobs=fitti(tem3,ktp3,tesr)*(temvhs3-temvcs3)/y0;
ektp=[0.034,0.05,0.077,0.128,0.199,0.294,0.423,0.669]; %расчетный эффективный КТП при параллельной засыпке
ete=[301,367,466,615,765,914,1064,1264];
for k=1:7 
    tte(k)=ete(k);
end
te0=273.15;
tktp=[1.44,1.66,2.16,3.19,4.45,6.14,8.15];
vte=arrTempAir()+te0;
ktpvo=koefTeploprovAir();
kktpvo=polyfit(vte,ktpvo,2);
kektp=polyfit(ete,ektp,2);
ktktp=polyfit(tte,tktp,2);
%p=plot(dite,ktp1,'-b',dite,ktp2,'-k',tem3,ktp3,'-r',tem3,ktp4,'-m');
%p=plot(tem3,ktp3,'-r',tem3,ktp4,'-m');
%p=plot(dite,proc,'-r',tem3,proc3,'-m');
%set(p,'LineWidth',1); hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Коэффициенты теплопроводности, Вт/(м*К)'}); 
%title({'График зависимости КТП от температуры'});
%legend('Стационарный метод', 'TCT','Location','Best');
te1=1000+te0; te2=(242+294)/2+te0;
%te1=temvhs3; te2=temvcs3;
ko0=(te2-te1)/y0; 
%deltaT=qob*y0/opredTeploprovVozd((te1+te2)/2)
ko1=te1;
h=12e-6;
hk=0;
por84=0.5575;
hv=h*por84/(1-por84);
ks=floor(1e-3/(hv+h));
tol=1e-3; 
Ns=floor((y0-hv)/(hv+h));
deltaT=0; xko=0; xna=0; tko=0; tesr=0; tett=0;
for q=1:30
    x=0; tex=0;
x(1)=hk+h/2; 
ht=h+hv;
tex(1)=polyval([ko0,ko1], x(1));
for k=2:ks
    x(k)=x(k-1)+ht;
    tex(k)=polyval([ko0,ko1],x(k));
end
xko(q)=x(ks);
xna(q)=x(1);
tko(q)=tex(ks);
dex=x(ks)-x(1);
det=tex(ks)-tex(1);
dtx(q)=abs(det/dex);
tna(q)=tko(q)-dtx(q)*dex;
tesr(q)=(tex(1)+tex(ks))/2;
lamvo=opredTeploprovVozd(tesr(q));
lampp(q)=fitti(tem3,ktp4,tesr(q));
tett(q)=opredktptvch(lampp(q),x,h,hv,qob,tna(q),xna(q),length(x));
hk=hk+1e-3;
deltaT=deltaT+qob*1e-3/lamvo;
end
tesr=tesr';
lampp=lampp';
deltaT=deltaT';
tett=tett';
%hk=RaschQ(tna(1),tem3,ktp3,[ko0 ko1],0:1e-5:30e-3,qob);
%hk=RasTemObraz(qob,[ko0,ko1],0:1e-4:30e-3,te1,ktp4,tem3,1e-4)
dlvo=1.3:1e-1:27;
for k=1:length(dlvo)
    nkbr(k)=OprPokPreKBr(dlvo(k));
end
nkbr=nkbr';
lam=[0];
end

function opn = OprPokPreKBr(la)
lam=[1      2       10      20      25]; 
pp=[1.487   1.48        1.52        1.48        1.453]; 
ko=polyfit(lam,pp,length(pp)-1); 
opn=polyval(ko,la);
end

function ktptt = opredktptvch(ektp,x,h,hv,qob,tn,xn,Ns)
tsre=0;
sdtma=0;
slam=0;
dto=0;
for k=1:Ns
dtk=qob*abs(x(k)-xn)/ektp;
dto=dto+dtk;
tk=tn-dtk;
tsre(k)=(tk+tn)/2;
lamvo=opredTeploprovVozd(tsre(k));
lam(k)=PoisKorn(lamvo,tk,qob,h,hv,tn);
slam=slam+lam(k)*tsre(k);
sdtma=sdtma+qob*hv/opredTeploprovVozd(tsre(k));
tn=tk;
xn=x(k);
end
tsre=tsre';
slam=slam/dto;
ktptt=slam;
end

function fi = fitti(arx,ary,xk)
f=1; q=length(arx); n=q;
for k=1:n
    if (f>0)
        if (arx(k)>xk)
            f=0;
            q=k;
            break;
        end
        if (arx(k)==xk)
            f=0;
            q=k;
            break;
        end
    end
end
if (f>0)
    ko=(ary(n)-ary(n-1))/(arx(n)-arx(n-1));
    ka=ary(n)+ko*(xk-arx(n));
end
if (f==0)
if (q==1)
ka=ary(1);
else
ko=(ary(q)-ary(q-1))/(arx(q)-arx(q-1));
ka=ko*(xk-arx(q-1))+ary(q-1);
end
end
fi=ka;
end

function kor = PoKo(ktp1,te1,ktp2,te2,ta,tb)
eps=1e-6;
Nkit=1e2;
dlin=abs(ta-tb);
h=0;
tc=(ta+tb)/2;
while (dlin>eps)
tc=(ta+tb)/2;
ktpa2=fitti(te2,ktp2,ta);
ktpb2=fitti(te2,ktp2,tb);
ktpc2=fitti(te2,ktp2,tc);
ktpa1=fitti(te1,ktp1,ta);
ktpb1=fitti(te1,ktp1,tb);
ktpc1=fitti(te1,ktp1,tc);
fc=ktpc2-ktpc1;
fa=ktpa2-ktpa1;
fb=ktpb2-ktpb1;
if ((fa*fc)<0)
if ((fb*fc)>0)
    tb=tc; 
end
end
if ((fa*fc)>0) 
if ((fb*fc)<0) 
    ta=tc; 
end
end
dlin=abs(ta-tb);
h=h+1;
if (h>Nkit) 
    break; 
end
end
kor=tc;
end

function kor = PoisKorn(lamvo,tk,qob,h,hv,tn)
lama=1e6;
lamb=1e-6;
ep=1e-8;
dlin=abs(lama-lamb);
Nkit=1e2;
k=0;
while (dlin>ep)
    lamc=(lama+lamb)/2;
    tka=tn-qob*(abs(h/lama)+abs(hv/lamvo));
    tkb=tn-qob*(abs(h/lamb)+abs(hv/lamvo));
    tkc=tn-qob*(abs(h/lamc)+abs(hv/lamvo));
    fa=tk-tka;
    fb=tk-tkb;
    fc=tk-tkc;
    if ((fc*fb)>0) 
        if ((fa*fc)<0) 
            lamb=lamc; 
        end
    end
    if ((fc*fa)>0) 
        if ((fb*fc)<0) 
            lama=lamc; 
        end
    end
    k=k+1;
    if (k>Nkit) 
        break; 
    end
    dlin=abs(lama-lamb);
end
kor=lamc;
end

function pr = RaschQ(tn,tem,ktp,koe,x,qob)
n=length(x);
qobs=0;
h=abs(x(1)-x(2)); temp(1)=tn;
for k=2:n
    lampp=fitti(tem,ktp,temp(k-1));
    dte(k)=qob*h/lampp;
    temp(k)=temp(k-1)-dte(k);
end
qva(1)=qob; qobs(1)=qob; 
sq=0; sx=0;
for k=2:n
    dx=abs(x(k)-x(k-1));
    lampp=fitti(tem,ktp,temp(k-1));
dt=abs(polyval(koe,x(k))-polyval(koe,x(k-1)));
qobs(k)=qob;
qva(k)=PoisKornQ(dt,lampp,dx);
sq=sq+qva(k)*dx;
sx=sx+dx;
temp(k)=tn-abs(dt);
tn=temp(k);
end
q=1; qm=0;
for k=1:n
    if (rem(k,200)==0)
       qm(q)=qva(k);
       q=q+1;
    end
end
qm=qm'
tol=abs(sx);
sq=sq/tol
qob'
qva=qva';
temp=temp';
%p=plot(x,qobs,'-b',x,qva,'-k');
%set(p,'LineWidth',1); hold on; grid on; 
%ylabel({'Плотность теплового потока, Вт/(м2)'}); 
%xlabel({'Координата, м'}); 
%title({'График зависимости ПТП от координаты q=q(x)'});
%legend('Метод 1', 'Метод 2','Location','Best');
pr=0;
end

function kor = PoisKornQ(dt1,lam,dx)
qoba=1e6;
qobb=1e-6;
ep=1e-8;
dlin=abs(qoba-qobb);
Nkit=1e2;
k=0;
while (dlin>ep)
    qobc=(qoba+qobb)/2;
    dtec=qobc*dx/lam;
    dtea=qoba*dx/lam;
    dteb=qobb*dx/lam;
    fa=dt1-dtea;
    fb=dt1-dteb;
    fc=dt1-dtec;
    if ((fc*fb)>0) 
        if ((fa*fc)<0) 
            qobb=qobc; 
        end
    end
    if ((fc*fa)>0) 
        if ((fb*fc)<0) 
            qoba=qobc; 
        end
    end
    k=k+1;
    if (k>Nkit) 
        break; 
    end
    dlin=abs(qoba-qobb);
end
kor=qobc;
end

function qt = RasTemObraz(qob,koe,x,tn,lam,tem,h)
n=length(x);
tep(1)=tn;
for k=2:n
    te(k)=polyval(koe,x(k));
    lamp=fitti(tem,lam,tep(k-1));
    tep(k)=te(k-1)-qob*h/lamp;
end
p=plot(x(4:n),te(4:n),'-b',x(4:n),tep(4:n),'-k');
set(p,'LineWidth',2); hold on; grid on; 
ylabel({'Температура, К'}); 
xlabel({'Координата, м'}); 
title({'График зависимости температурыот координаты T=T(x)'});
legend('Метод 1', 'Метод 2','Location','Best');
qt=0;
end