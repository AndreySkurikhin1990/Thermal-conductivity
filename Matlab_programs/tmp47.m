%Промежуточные расчеты
function [ t ] = tmp47()
t=nacha();
end
function [ tm ] = nacha()
format long g; %определяет коэффициент рассеяния для таблетки KBr+вермикулит
t=2*PoisKorn();
ti=RasFra150_200(t);
tm=ti;
end
function [ pk ] = PoisDn(vyfr,vyko,xv,vv,mo,si,mi,ma)
TrRef=SredGrafRass(vyfr,vyko,xv); dlv=1e6*dlivoln(); numele=1e0; na=1e0; 
xve=opredEffPar(xv,vyfr,vyko); vve=opredEffPar(vv,vyfr,vyko); dnre=0; q=1;
for w=na:na+numele-1
    dlvo=dlv(w); TR=TrRef(w); n1=OprPokPreKBr(dlvo); dna=-1e0; dnb=2e0; ep=1e-6; Nit=1e0; k=0; ra=abs(dna-dnb); w=w'
while ((ra>ep) && (k<Nit))
    dnc=(dna+dnb)/2;
    na2=n1+dna;
    if (dna==0)
        fa=1;
    else
    fa=opreKoefProp(vve,mo,si,mi,ma,dlvo,n1,na2); 
    end
    fa=fa-TR;
    nb2=n1+dnb;
    if (dnb==0)
        fb=1;
    else
    fb=opreKoefProp(vve,mo,si,mi,ma,dlvo,n1,nb2); 
    end
    fb=fb-TR;
    nc2=n1+dnc;
    if (dnc==0)
        fc=1;
    else
    fc=opreKoefProp(vve,mo,si,mi,ma,dlvo,n1,nc2); 
    end
    fc=fc-TR;
    ro=vybPoKo(fa,fb,fc,dna,dnb,dnc); 
    dna=ro(1); dnb=ro(2); dnc=ro(3); 
    k=k+1;
    ra=abs(dna-dnb);
end
dnre(q)=dnc; q=q+1;
end
pk=dnre;
end
%Данная функция должна возвращать вероятность попадания при выстрелах, mo - матем. ожид-е, si - СКО, mi, ma - мин. и макс. разм.
function mdn = opreKoefProp(vvk,mo,si,mi,ma,dlvo,n1,n2)
Nr=1e0; raot=13e3; blis=0; ocp=0; sepoo=0; rosr=0; %raot - размер отверстия, blis - число удачных попаданий
for g=1:Nr %усреднение по реализациям
suv=0; k=1; nk=2e0;
while ((suv<vvk) && (k<nk))
hk=mo+si*randn(); hk=proverka(hk,mi,ma); %определение длины направляющей
tetani=pi*rand(); %направление излучения - угол между плоскостью основания и осью Z (вверх)
phini=pi*(-1/2+rand()); %угол между проекцией оси цилиндра на горизонтальную ось и осью X - коллинеарно направлению излучения (ось X вниз)
nalu=[0    0   -1]; nalusk=PovorotyIzNovKIskh(nalu,tetani,phini); %направление излучения в старых координатах
ocp=ocp+1; %общее число прицеливаний
    dk=mo+si*randn(); dk=proverka(dk,mi,ma); rk=dk/2; %определение диаметра цилиндра
    sob=opredPloshOsnBokPov(rk,tetani,phini,hk);
    suv=suv+pi*hk*(rk^2); %проецирование на вертикальную ось, перпендикулярную направлению излучения
    %se=sob(1); sb=sob(2); rva=sob(3); rvb=sob(4); hv=sob(5); phima=sob(6); phimi=sob(7); 
    %sep=se+sb; sepoo=sepoo+sep; %сечение рассеяния
    %vys=VysLucha(se,sep,tetani,hk,rk,n1,n2,dlvo,raot,phini,rva,rvb,hv,nalusk,phima,phimi); 
    %blis=blis+vys; %число благополучных исходов, 1 - успешно, 0 - нет
k=k+1;
end
end
w=blis/ocp; rosr=rosr/Nr;
p=sepoo/Nr/(pi*(raot^2)/4);
mdn=p*w*rosr+(1-p); %пропускание за счет рассеяния
end
%выстрел луча по цилиндру - в плоскости угла падения
function vy = VysLucha(se,spo,teta,hk,rk,n1,n2,dlv,razotv,phi,va,vb,hv,nisk,phima,phimi)
s=spo*rand(); po=0;
if (s<se)
    p=1;
    ktpnk=PoisKTPN(va,vb); %координаты точки падения луча на основание в вертикальной плоскости
    ktpsk=PovorotyIzNovKIskh(ktpnk,teta,phi); %в старых координатах
    naluskoc=[nisk(1)         nisk(2)         0]; %направление излучения на основании цилиндра
    phiup=acos(nisk(3)); %угол падения
    phiups=PoiskTetaShtrNach(phiup,n1,n2);  %phiups - угол преломления, phiup - угол падения
    v=VykhodIzCylOsn(naluskoc,ktpsk,teta,phi,hk,rk); %ищем, откуда будет выходить луч, в плоскости самого основания
    if (v==2)
    phiupss=PoiskTetaShtrNach(pi/2-phiup,n1,n2);  %выходит через боковую поверхность
    po=PopNePop(dlv,abs(phiup-phiupss),razotv); %попадание в диапазон малых углов
    else po=1;  %через основание - никакого ухода луча нет
    end
else p=2; nu=pi/2-teta;
    ktpnk=PoisKTPNBok(hv,rk,phima,phimi,phi,teta); %координаты точки падения луча на основание в вертикальной плоскости
    ktpsk=PovorotyIzNovKIskh(ktpnk,teta,phi); %в старых координатах
    phitp=atan(ktpsk(2)/ktpsk(1)); 
    nalu=[0  0   -1]; nalusk=PovorotyIzNovKIskh(nalu,teta,phi); 
    phip=atan(nalu(2)/nalu(1));
    if (((phitp==0) || (phitp==pi)) && ((phip==0) || (phip==pi/2) || (phip==-pi/2)))
    nus=PoiskTetaShtrNach(abs(nu),n1,n2);
    else
        nus=PoiskTetaShtrNach(abs(phip),n1,n2); %угол преломления в плоскости сечения, параллельно основанию
        phipz=acos(nalusk(3)/rk); %угол между направлением излучения и образующей в касательной плоскости, параллельной оси Z
    end
    v=VykhodIzCylBokPov(n1,n2,-hk*rand(),2*rk,hk,nus,phip,phipz,phitp); %1 - выход через основание, 2 - через бок. пов-ть
    if (v>1)
        po=PopNePop(dlv,2*abs(phip-nus),razotv); %попадание в диапазон малых углов - выход через бок. пов-ть
    else
        nuss=PoiskVykhUglaSSKonech(dn,abs(nus));  %выходит через основание
        po=PopNePop(dlv,abs(phipn-nuss),razotv); %попадание в диапазон малых углов
    end
end
vy=po;
end
function [ kup ] = PoisKTPN(av,bv)
ry=bv*(-1+2*rand());
pn=asin(ry/bv);
pk=sign(pn)*(pi-pn);
phi=pn+(pk-pn)*rand();
rx=av*cos(phi);
kup=[rx     ry      0];
end
function [ kup ] = PoisKTPNBok(hv,r,phima,phimi,phi,teta)
    kup=[0  0   0];
end
%падает и выходит через основание или боковую стенку
function pr10 = PraCha10(phi,phis,n1,n2)
pr10=n1*sin(phi)-n2*sin(phis);
end
%показатель преломления бромида калия
function opn = OprPokPreKBr(la)
lam=[1      2       10      20      25]; 
pp=[1.487   1.48        1.52        1.48        1.453]; 
ko=polyfit(lam,pp,length(pp)-1); 
opn=polyval(ko,la);
end
%поиск угла преломления, phi>0
function pk = PoiskTetaShtrNach(phi,n1,n2)
a=0; b=pi/2; ep=1e-7; Nit=1e2; k=0; ra=abs(a-b);
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; 
    fa=PraCha10(phi,a,n1,n2);
    fb=PraCha10(phi,b,n1,n2); 
    fc=PraCha10(phi,c,n1,n2);
    ro=vybPoKo(fa,fb,fc,a,b,c); 
    a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
pk=c;
end
%------------------пока не трогать--------------
function pd = PopNePop(lam,deltatetass,raot)
p=0; deko=lam/raot;
if (deltatetass>deko)
    p=0;
else p=1;
end
pd=p;
end
function v = VykhodIzCylBokPov(n1,n2,z,dk,hk,hi,phip,phipz,phitp)
vyh=1; dn=n2-n1;
if (((phitp==pi/2) || (phitp==-pi/2)) && (phip==0) && (phipz==0))
    vyh=1;
end
if (((phitp==pi/2) || (phitp==3*pi/2)) && (phip==0) && (phipz==0))
    vyh=1;
end
if (phip==0)
beta=atan(abs((hk-z)/dk));
if (abs(phipz)>beta)
    vyh=1;
else vyh=2;
end
if ((dn<0) && (abs(phip)>beta))
        vyh=1; %1 - выход через основание, 2 - через бок. пов-ть
        else vyh=2;
end
end
if (hi==pi/2)
    vyh=2;
end
dlpuos=dk*abs(cos(hi));
dlpuob=abs(hk-z); 
if (dlpuob<dlpuos)
    vyh=1;
else vyh=2;
end
v=vyh;
end
function v = VykhodIzCylOsn(nalu,tpisk,teta,phi,h,r)
    dlo=PoiskDlPutiLuchaOsn(r,nalu,tpisk);
    if (dlo<h)
        vyh=2; %выход из боковой поверхности
    else vyh=1; %выход из основания
    end
    vyh=OprTochVykhKritSluch(phi,teta,nalu,tpisk,h,r); %в плоскости основания
v=vyh;
end
function [ v ] = vybPoKo(fa,fb,fc,a,b,c)
if ((fc*fb)>0) 
        if ((fa*fc)<0) 
            b=c; 
        end
end
    if ((fc*fa)>0) 
        if ((fb*fc)<0) 
            a=c; 
        end
    end
    v=[a        b       c];
end
function [ alsr ] = SredGrafRass(vyfr,vyko,xv)
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
%vyfr - выбор фракции: 1 - 60, 2 - 60-100, 3 - 100-150, 4 - 150-200; vyko - выбор массовой доли: 1 - 0,1 %, 2 - 0,2 %, 3 - 0,3 %
if (vyfr==4)
Tal=privedkEdiPropus(TrkoVer5681()); Tal2=privedkEdiPropus(TrkoVer5682()); Tal3=privedkEdiPropus(TrkoVer5683()); 
Tal4=privedkEdiPropus(TrkoVer5691()); Tal5=privedkEdiPropus(TrkoVer5692()); Tal6=privedkEdiPropus(TrkoVer5693());
alsre=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv,vyko);
end
alsr=alsre;
end
function [ prked ] = privedkEdiPropus(ar)
arr=ar;
p=length(arr);
for k=1:p
    if (arr(k)>1)
        arr(k)=1;
    end
    if (arr(k)<0)
        arr(k)=0;
    end
end
prked=arr;
end
function [ al ] = opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv,vyko)
nTal=length(Tal);
for k=1:nTal 
    Ta(k)=-log(Tal(k))/xv(1); Ta2(k)=-log(Tal2(k))/xv(2); Ta3(k)=-log(Tal3(k))/xv(3); Ta4(k)=-log(Tal4(k))/xv(4); Ta5(k)=-log(Tal5(k))/xv(5); Ta6(k)=-log(Tal6(k))/xv(6); 
end
switch (vyko)
    case 1
alsr=opredAlphSred2(Tal,Tal4);
    case 2
alsr=opredAlphSred2(Tal2,Tal5);
    case 3
alsr=opredAlphSred2(Tal3,Tal6);
end
al=alsr;
end
function [ al ] =opredAlphSred2(Tal,Tal2)
alsr=0; nTal=length(Tal);
for k=1:nTal 
    alsr(k)=Tal(k)+Tal2(k);
end
al=alsr/2;
end
function pk = PoisKorn()
ro=0; a=1e-5; b=1e5; ep=1e-6; Nit=1e2; k=0; ffv=1e-4; ra=abs(a-b); %ffv = 0,01 % - такая доля частиц не попадает в диапазон фракции
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; 
    fa=erfc(a)-ffv; 
    fb=erfc(b)-ffv; 
    fc=erfc(c)-ffv; 
    ro=vybPoKo(fa,fb,fc,a,b,c); 
    a=ro(1); b=ro(2); c=ro(3); k=k+1;
    ra=abs(a-b);
end
pk=c;
end
function [ ras150_200 ] = RasFra150_200(di)
mkbr=[250.882 249.590 249.213 250.299 249.441 249.365];
mv=[0.320 0.533 0.849 0.223 0.502 0.797]; 
tol=[0.76 0.72 0.69 0.73 0.73 0.73]; 
rokbr=2.75; rov=0.56; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); vv(k)=mv(k)/(1e3*rov); xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vv(k)=vv(k)*((1e4)^3); vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
end
mini=15e1; maxi=2e2; mo=(maxi+mini)/2; %МО
si=(mini-maxi)/di; %СКО
disp('Фракция 150-200 мкм');
depp=PoisDn(4,1,xv,vv,mo,si,mini,maxi);
ras150_200=depp; 
end
function p = proverka(hk,mi,ma)
hp=hk;
if (hk<mi)
    hp=mi;
end
    if (hk>ma)
        hp=ma;
    end
    p=hp;
end
function oep = opredEffPar(mapa,vyfr,vyko)
switch (vyfr)
    case {1,4}
        ne=2; 
        switch (vyko)
            case 1
                no = [1 4];
            case 2
                no = [2 5];
            case 3
                no = [3 6];
        end
    case {2,3}
        ne=3;
        switch (vyko)
            case 1
                no = [1 4 7];
            case 2
                no = [2 5 8];
            case 3
                no = [3 6 9];
        end
end
s=0;
for k=1:ne
    s=s+mapa(no(k));
end
oep=s/ne;
end
%phin - угол направления падения, xp, yp - координаты падения на круге, ищем в плоскости основания
function dlot = PoiskDlPutiLuchaOsnCyl(rok,nalu,tpisk)
ko=nalu(2)/nalu(1);
if ((isinf(ko)) || (isnan(ko)))
    dlo=0;
else
xtope(1)=0; ytope(1)=0; xtope(2)=0; ytope(2)=0;
dlo=0; ep=1e-6; Nit=1e2; 
yp=tpisk(2); xp=tpisk(1); sc=yp-ko*xp;
b=(2*ko*sc); a=(ko^2+1); c=(sc^2-rok^2); diskr=(b^2)-4*a*c;
if (diskr<0)
    rvp(1)=0; rvp(2)=0; rvv(1)=0; rvv(2)=0;
else
    diskr=sqrt(diskr);
    xtope(1)=(-b-diskr)/2/a; 
    ytope(1)=ko*xtope(1)+sc; 
    rvp=[xtope(1)       ytope(1)        0];
    xtope(2)=(-b+diskr)/2/a; 
    ytope(2)=ko*xtope(2)+sc; 
    rvv=[xtope(2)      ytope(2)        0];
end
nalu=nalu';
rvp=rvp';
rvv=rvv';
if ((nalu*rvp')>0)
    dxyz=[rvp(1)-xp       rvp(2)-yp        0]; 
    dlo=DlinaVectora(dxyz);
else
if ((nalu*rvv')>0)
    dxyz=[rvv(1)-xp      rvv(2)-yp        0]; 
    dlo=DlinaVectora(dxyz);
else
    dlo=0;
end
end
end
dlot=dlo;
end
function [ nisk ] = PovorotyIzNovKIskhCyl(vnk,teta,phi)
ksi=vnk(1); eta=vnk(2); dzeta=vnk(3);
xs=ksi; ys=eta*cos(phi)-dzeta*sin(phi); zs=eta*sin(phi)+dzeta*cos(phi);
x=xs*sin(teta)+zs*cos(teta); y=ys; z=-xs*cos(teta)+zs*sin(teta);
nisk=[x     y   z];
end
function [ vnk ] = PovorotyIzIskhKNovCyl(vik,teta,phi)
x=vik(1); y=vik(2); z=vik(3);
xs=x*sin(teta)-z*cos(teta); ys=y; zs=x*cos(teta)+z*sin(teta);
ksi=xs; eta=ys*cos(phi)+zs*sin(phi); dzeta=-ys*sin(phi)+zs*cos(phi);
vnk=[ksi     eta   dzeta];
end
function dlve = DlinaVectora(vect)
s=0;
for k=1:length(vect)
    s=s+(vect(k))^2;
end
dlve=sqrt(s);
end
function [ vnnk ] = PovorotyIzIskhKNovBok(vnsk,teta)
x=vnsk(1); y=vnsk(2); z=vnsk(3);
xs=x*cos(teta)+z*sin(teta); zs=-x*sin(teta)+z*cos(teta); ys=y;
vnnk=[xs    ys  zs];
end
function [ pobp ] = opredPloshOsnBokPov(r,teta,phi,h)
b=0; nu=0:1e-2:2*pi;
pc=[0  0   -h];
pc=PovorotyIzIskhKNov(pc,teta,phi); 
rpc=DlinaVectora(pc); %перенос центра
for k=1:length(nu)
    nuk=nu(k); vxyz=[r*cos(nuk)  r*sin(nuk)   0];
    vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); 
    x0(k)=vxyz(1); y0(k)=vxyz(2); vxyz(3)=0; 
    b(k)=DlinaVectora(vxyz);
end
p=plot(x0,y0,':b');
set(p,'LineWidth',2); hold on; grid on; 
xlabel({'x'}); 
ylabel({'y'}); 
title({'y(x)'});
ma=max(b); phima=PoiskUglaMax(nu,b,ma,length(nu),phi,teta,r)
vxyz=[r*cos(phima)  r*sin(phima)   0]; 
vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); %находим полуоси эллипса - проекции основания на вертикальной плокости
vxyz(3)=0; 
ma=DlinaVectora(vxyz); %большая полуось
phimi=phima+pi/2; 
vxyz=[r*cos(phimi)  r*sin(phimi)   0];
vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); 
vxyz(3)=0;
mi=DlinaVectora(vxyz); %малая полуось
so=pi*ma*mi; sbp=rpc*ma; 
pobp=[so   sbp     ma  mi  rpc     phima   phimi];
end
function [ vy ] = VybABC(fa,fb,fc,a,b,c)
if (fc>fa) 
        if (fc>fb)
            if (fa>fb)
                b=c;
            else
                a=c;
            end
        else
        a=c;    
        end
else
        b=c;
end
vy=[a   b   c];
end
function pk = PoiskUglaMax(nu,ba,m,n,phi,teta,ak)
diap=PoiskDiap(nu,ba,m,n);
a=diap(1); b=diap(2); c=(a+b)/2;
ep=1e-6; Nit=1e2; k=0; ra=abs(a-b); 
while ((ra>ep) && (k<Nit))
    c=(a+b)/2;
    xya=PovorotyIzIskhKNov([ak*cos(a)   ak*sin(a)    0],teta,phi); xya(3)=0; fa=DlinaVectora(xya);
    xyb=PovorotyIzIskhKNov([ak*cos(b)   ak*sin(b)    0],teta,phi); xyb(3)=0; fb=DlinaVectora(xyb);
    xyc=PovorotyIzIskhKNov([ak*cos(c)   ak*sin(c)    0],teta,phi); xyc(3)=0; fc=DlinaVectora(xyc);
    ro=VybABC(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3);
ra=abs(a-b); k=k+1;
end
pk=c;
end
function [ diap] = PoiskDiap(phi,ar,m,n)
f=1; l=1; p=1;
for k=2:(n-1)
    if (ar(k)==m)
        if (f>0)
        l=k-1;
        p=k+1;
        f=0;
        end
    end
end
ug(1)=phi(l); ug(2)=phi(p);
a=min(ug); b=max(ug);
diap=[a     b];
end
function vykh = OprTochVykhKritSluch(phi,teta,nalu,tpisk,h,r) %в плоскости основания
phin=atan(nalu(2)/nalu(1)); phitp=atan(tpisk(2)/tpisk(1)); ep=1e-6;
vyh=1; ro=0;
if ((phi==0) && (teta==0))
vyh=1;
end
if ((phin==0) || (phin==pi) || (phin==2*pi) || (phin==-pi) || (phin==-2*pi) || (phin==pi/2) || (phin==3*pi/2) || (phin==-pi/2) || (phin==-3*pi/2))
    if ((phitp==0) || (phitp==pi) || (phitp==2*pi) || (phitp==-pi) || (phitp==-2*pi) || (phitp==pi/2) || (phitp==3*pi/2) || (phitp==-pi/2) || (phitp==-3*pi/2))
        if (abs(nalu*tpisk')>ep)
            if (abs(tpisk(2))>ep)
                if ((nalu*tpisk')>0)
            dy=r-tpisk(2); dx=tpisk(1); ro=sqrt(dx^2+dy^2);
                else
                dy=r-tpisk(2); dx=tpisk(1); ro=sqrt(dx^2+dy^2);
                end
            else
                if ((nalu*tpisk')>0)
                dx=r-tpisk(1); dy=tpisk(2); ro=sqrt(dx^2+dy^2);
                else
                dx=r-tpisk(1); dy=tpisk(2); ro=sqrt(dx^2+dy^2);
                end
            end
        else
            x=tpisk(1); ro=sqrt(r^2-x^2);
        end
    end
end
if (ro<h)
    vyh=2; %выход из боковой поверхности
    else vyh=1; %выход из основания
end
vykh=vyh;
end
%проверка - промежуточные расчеты
function t = tmp47_1()
format long g; 
maxz=2e2;
minz=16e1;
a=minz+(maxz-minz)*rand(); 
h=minz+(maxz-minz)*rand();
teta=pi*rand(); 
phi=pi*(-1/2+rand());
ktp=PoisKTPNBokPodProg(a,phi,teta,h);
t=0;
end
function [ te ] = tempor()
nu=0:1e-1:2*pi; n=length(nu); pc=PovorotyIzIskhKNov([0  0   -h],teta,phi); rpc=sqrt(((pc(1))^2)+((pc(2))^2));
for k=1:n
    b(k)=0; nuk=nu(k); vxyz=PovorotyIzIskhKNov([a*cos(nuk)  a*sin(nuk)   0],teta,phi); x=vxyz(1); y=vxyz(2); b(k)=sqrt((x^2)+(y^2));
end
ma=max(b); mi=min(b); phima=PoiskUglaMax(nu,b,ma,n,phi,teta,a);
vxyz=PovorotyIzIskhKNov([a*cos(phima)  a*sin(phima)   0],teta,phi);
vxyz=PoisKTPNBok(h,a,phima,phima+pi/2,phi,teta);
ak=1.8; phitp=3*pi/4; nri(1)=-1; nri(2)=-1; tpisk(1)=ak*cos(phitp); tpisk(2)=ak*sin(phitp); plo=PoiskDlPutiLuchaOsn(a,nri,tpisk);
phiup=pi/6; dn=1e-1;
phis=PoiskTetaShtrNach(phiup,dn);
phiss=PoiskTetaShtrNach(pi/2-phis,-dn);
nalusk=PovorotyIzNovKIskh([0    0   -1],teta,phi); %направление излучения в старых координатах
vys=VysLucha(teta,phi,nalusk,ma,mi);
te=[0];
end
function phixy = PrivkOdnoZnach(ug,phinul)
phix1=real(acos(ug(1)));
phix2=real(asin(ug(2)));
if (phix1>(pi/2))
    if (phix2<0)
        phix1=pi-phix1; %x<0, y<0, 3 четверть
        phix1=pi+phix1;
        phix2=abs(phix2);
        phix2=pi+phix2;
    else
        phix2=pi-phix2; %x<0, y<0, 2 четверть
    end
else
    if (phix2<0)
        phix1=-phix1+2*pi; %x>0, y<0, 4 четверть
        phix2=phix2+2*pi;
    end
end
phix1=phix1';
phix2=phix2';
phixy=(phix1+phix2)/2+phinul;
end
function [ ktpl ] = PoisKTPNBokPodProg(r,phi,teta,h)
rpc=[0  0   -h]; rpc=PovorotyIzIskhKNov(rpc,teta,phi);
phima=PoiskUglaMax(phi,teta,r); vxyz=[r*cos(phima)     r*sin(phima)    0]; vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); vxyz(3)=0; vxyz=Pereobozn(vxyz); 
vxyz=PlusOdinPovorot(vxyz); k=1; phim(k)=phima; xym(k)=vxyz(2); vxyz=vxyz';
phima=phima+pi; vxyz=[r*cos(phima)      r*sin(phima)    0]; vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); vxyz(3)=0; vxyz=Pereobozn(vxyz); 
vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phima; xym(k)=vxyz(2); vxyz=vxyz';
phimi=phima-pi/2; vxyz=[r*cos(phimi)      r*sin(phimi)    0]; vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); vxyz(3)=0; vxyz=Pereobozn(vxyz); 
vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phimi; xym(k)=vxyz(2); vxyz=vxyz';
phimi=phimi+pi; vxyz=[r*cos(phimi)      r*sin(phimi)    0]; vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); vxyz(3)=0; vxyz=Pereobozn(vxyz); 
vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phimi; xym(k)=vxyz(2); vxyz=vxyz'; xym=xym';
if (xym(1)>xym(2))
        minx=xym(2); phimix=phim(2); maxx=xym(1); phimax=phim(1);
else
        minx=xym(1); phimix=phim(1); maxx=xym(2); phimax=phim(2); 
end
if (xym(3)>xym(4))
        miny=xym(4); maxy=xym(3); phimay=phim(3); phimiy=phim(4);
else
        miny=xym(3); maxy=xym(4); phimay=phim(4); phimiy=phim(3);
end
    nxy=1e2; hnu=(phimax-phimix)/nxy; numa=phimix:hnu:phimax;
    xyz=[r*cos(phimix)       r*sin(phimix)        0]; xyz=PovorotyIzIskhKNov(xyz,teta,phi); xyz(3)=0; xymi=xyz;
    xyz=[r*cos(phimax)       r*sin(phimax)        0]; xyz=PovorotyIzIskhKNov(xyz,teta,phi); xyz(3)=0; xyma=xyz;
    rxy=(xymi+xyma)/2; phis=(phimax+phimix)/2; rxys=[r*cos(phis)       r*sin(phis)        0]; rxys=PovorotyIzIskhKNov(rxys,teta,phi); rxys(3)=0;
    if (rxys(1)>rxy(1))
        phimax=phimix; phimix=phimix-pi; hnu=(phimax-phimix)/nxy;
    end
        xyz=[r*cos(phimax)       r*sin(phimax)        0]; xyzmas=xyz'; xyz=[r*cos(phimix)       r*sin(phimix)        0]; xyzmis=xyz';
        raxys=xyzmas-xyzmis; draxys=DlinaVectora(raxys); psiras=PrivkOdnoZnach([real(raxys(1)/draxys)         real(raxys(2)/draxys)],0);
    xyz=[r*cos(phimax)       r*sin(phimax)        0]; xyz=PovorotyIzIskhKNov(xyz,teta,phi); xyz(3)=0; xyz=Pereobozn(xyz); xyzma=xyz'; dxyzma=DlinaVectora(xyzma);
    psima=PrivkOdnoZnach([real(xyzma(1)/dxyzma)         real(xyzma(2)/dxyzma)],0); xyzmi=-xyzma; psimi=psima-pi;
    raxy=xyzma-xyzmi; draxy=DlinaVectora(raxy); psira=PrivkOdnoZnach([real(raxy(1)/draxy)         real(raxy(2)/draxy)],0);
    maxzx=(abs(maxx)+abs(minx))/2; maxzy=(abs(maxy)+abs(miny))/2; hpsi=(psima-psimi)/nxy;
    spsis=(psimi+psima)/2; naspsis=[cos(spsis)      sin(spsis)      0]; phis=(phimax+phimix)/2; 
    nasphis=[cos(phis)      sin(phis)      0]; ugpsisphis=acos((naspsis*nasphis')/DlinaVectora(naspsis)/DlinaVectora(nasphis));
    ugpsi=psimi:hpsi:psima; %выбор угла в основании
    rpcp=rpc*rand(); rpcp=PovorotyIzNovKIskh(rpcp,teta,phi);
    for k=1:length(ugpsi)
%прямое преобразование из старых координат в новые   
    psik=ugpsi(k)-psira;
    xyz=[maxzx*cos(psik)      maxzy*sin(psik)       0];
    xyz=DopolnitPovorot(xyz,psira);
    if ((ugpsisphis>pi/2) || (phimax==0))
        xyz=DopolnitPovorot(xyz,pi);
    end
    xs(k)=xyz(1); ys(k)=xyz(2); %новые координаты
    xyz=[r*cos(numa(k))       r*sin(numa(k))        0]; 
    if (sign(hpsi)==-sign(hnu))
        xyz=DopolnitPovorot(xyz,pi);
    end
    xv(k)=xyz(1); yv(k)=xyz(2); %старые координаты
    xyz=PovorotyIzIskhKNov(xyz,teta,phi); xyz(3)=0; xyz=Pereobozn(xyz); xn(k)=xyz(1); yn(k)=xyz(2); %новые координаты
%обратное преобразование из новых координат в старые
    if (sign(hpsi)==-sign(hnu))
        xyz=DopolnitPovorot(xyz,pi);
    end
xyz=DopolnitPovorot(xyz,-psira);
nuk=PrivkOdnoZnach([real(xyz(1)/maxzx)         real(xyz(2)/maxzy)],psira+(phimi-psimi));
        xyz=[-max([maxzx  maxzy])*sin(nuk)      max([maxzx  maxzy])*cos(nuk)       0];
        if ((ugpsisphis>pi/2) || (phimax==0))
            xyz=DopolnitPovorot(xyz,pi);
        end
        xw(k)=xyz(1);
        yw(k)=xyz(2); %старые координаты
    end
    %nuk=PrivkOdnoZnach([real(xw/r)         real(yw/r)],0)
        phimix
        phimax
        psimi
        psima
        hpsi
        hnu
        ugpsisphis
    pl=plot(xw,yw,':r',xv,yv,':b'); 
%pl=plot(xs,ys,':r',xn,yn,':b'); %новые координаты
set(pl,'LineWidth',4); grid on; 
xlabel({'x'}); ylabel({'y'}); title({'y(x)'});
ktpl=rpcp;
end
function [ mpv ] = Pereobozn(axy)
x=axy(1); y=axy(2); z=axy(3);
ksi=y; eta=-x; dzeta=z;
mpv=[ksi    eta         dzeta];
end
function [ mpv ] = PereoboznObrat(axy)
ksi=axy(1); eta=axy(2); dzeta=axy(3);
x=-eta; y=ksi; z=dzeta;
mpv=[x    y   z];
end
function [ mpv ] = DopolnitPovorot(axy,phi)
x=axy(1); y=axy(2); z=axy(3);
ksi=x*cos(phi)-y*sin(phi); eta=x*sin(phi)+y*cos(phi); dzeta=z;
mpv=[ksi    eta         dzeta];
end
function [ mpv ] = PlusOdinPovorot(axy)
x=axy(1); y=axy(2); phiksi=x/y; phiksi=atan(phiksi);
ksi=x*cos(phiksi)-y*sin(phiksi); eta=x*sin(phiksi)+y*cos(phiksi);
mpv=[ksi    eta];
end
function pk = PoiskUglaMax(phi,teta,ak)
bk=0; nu=0:1e-2:2*pi; n=length(nu);
for k=1:n
    nuk=nu(k); vxyz=[ak*cos(nuk)        ak*sin(nuk)         0];
    vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); 
    vxyz=Pereobozn(vxyz);
    x0(k)=vxyz(1); y0(k)=vxyz(2); vxyz(3)=0; 
    bk(k)=DlinaVectora(vxyz);
end
ma=max(bk);
diap=PoiskDiap(nu,bk,ma,n);
a=diap(1); b=diap(2); c=(a+b)/2;
ep=1e-6; Nit=1e2; k=0; ra=abs(a-b); 
while ((ra>ep) && (k<Nit))
    c=(a+b)/2;
    xya=PovorotyIzIskhKNov([ak*cos(a)   ak*sin(a)    0],teta,phi); xya(3)=0; fa=DlinaVectora(xya);
    xyb=PovorotyIzIskhKNov([ak*cos(b)   ak*sin(b)    0],teta,phi); xyb(3)=0; fb=DlinaVectora(xyb);
    xyc=PovorotyIzIskhKNov([ak*cos(c)   ak*sin(c)    0],teta,phi); xyc(3)=0; fc=DlinaVectora(xyc);
    ro=VybABC(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3);
ra=abs(a-b); k=k+1;
end
pk=c;
end
function [ diap] = PoiskDiap(phi,ar,m,n)
f=1; l=1; p=1;
for k=2:(n-1)
    if (ar(k)==m)
        if (f>0)
        l=k-1;
        p=k+1;
        f=0;
        end
    end
end
ug(1)=phi(l); ug(2)=phi(p);
a=min(ug); b=max(ug);
diap=[a     b];
end
function dlot = PoiskDlPutiLuchaOsnPP(rok,nalu,tpisk)
phin=nalu(2)/nalu(1); phin=atan(phin);
ko=tan(phin); dlo=0; ep=1e-6; Nit=1e2;
yp=tpisk(2); xp=tpisk(1); sc=yp-ko*xp;
b=(2*ko*sc); a=(ko^2+1); c=(sc^2-rok^2);
diskr=b^2-4*a*c;
if (diskr<0)
    xtope(1)=0;
    xtope(2)=0;
    ytope(1)=0;
    ytope(2)=0;
else
    xtope(1)=(-b-sqrt(diskr))/2/a;
    xtope(2)=(-b+sqrt(diskr))/2/a;
    ytope(1)=ko*xtope(1)+sc;
    rvp(1)=xtope(1); rvp(2)=ytope(1);
    ytope(2)=ko*xtope(2)+sc;
    rvv(1)=xtope(2); rvv(2)=ytope(2);
end
if (nalu*rvp'>0)
    vy=1;
else
if (nalu*rvv'>0)
    vy=2;
else
    vy=0;
end
end
phi=0:1e-3:2*pi; n=length(phi);
for k=1:n
    phik=phi(k);
    xk(k)=rok*cos(phik);
    yk(k)=rok*sin(phik);
end
xl=-2*a:1e-3:2*a; n=length(xl);
for k=1:n
    yl(k)=ko*xl(k)+sc;
end
%pl=plot(xk,yk,'-k',xl,yl,'-b');
%set(pl,'LineWidth',2); hold on; grid on; 
%xlabel({'x'}); 
%ylabel({'y'}); 
%title({'График y(x)'});
%legend('y(x)', 'z(x)', 'z(y)','location','best');
dlot=0;
end
function [ v ] = vybPoKo(fa,fb,fc,a,b,c)
if ((fc*fb)>0) 
        if ((fa*fc)<0) 
            b=c; 
        end
end
    if ((fc*fa)>0) 
        if ((fb*fc)<0) 
            a=c; 
        end
    end
    v=[a b c];
end
function [ nisk ] = PovorotyIzNovKIskhPP(vnk,teta,phi)
ksi=vnk(1); eta=vnk(2); dzeta=vnk(3);
xs=ksi; ys=eta*cos(phi)-dzeta*sin(phi); zs=eta*sin(phi)+dzeta*cos(phi);
x=xs*sin(teta)+zs*cos(teta); y=ys; z=-xs*cos(teta)+zs*sin(teta);
nisk=[x     y   z];
end
function [ vnk ] = PovorotyIzIskhKNovPP(vik,teta,phi)
x=vik(1); y=vik(2); z=vik(3);
xs=x*sin(teta)-z*cos(teta); ys=y; zs=x*cos(teta)+z*sin(teta);
ksi=xs; eta=ys*cos(phi)+zs*sin(phi); dzeta=-ys*sin(phi)+zs*cos(phi);
vnk=[ksi     eta   dzeta];
end
function dlve = DlinaVectora(vect)
s=0;
for k=1:length(vect)
    s=s+(vect(k))^2;
end
dlve=sqrt(s);
end
function vy = VysLucha(teta,phi,nisk,va,vb)
    ktpnk=PoisKTPN(va,vb); %координаты точки падения луча на основание в вертикальной плоскости
    ktpsk=PovorotyIzNovKIskh(ktpnk,teta,phi); %в старых координатах
    phiup=nisk(3); phiup=acos(phiup); phiup=(phiup+acos(-sin(teta)*cos(phi)))/2;
    vy=0;
end
function [ kup ] = PoisKTPN(av,bv)
ry=bv*(-1+2*rand());
pn=asin(ry/bv);
pk=sign(pn)*(pi-pn);
phi=pn+(pk-pn)*rand();
rx=av*cos(phi);
kup=[rx     ry  0];
end
function pk = PoiskTetaShtrNach(phi,dnr)
a=0; b=pi/2; ep=1e-6; Nit=4e1; k=0; ra=abs(a-b);
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; fa=PraCha10(phi,a)-dnr; fb=PraCha10(phi,b)-dnr; fc=PraCha10(phi,c)-dnr;
    ro=vybPoKo(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
pk=c;
end
function pr10 = PraCha10(phi,phis)
pr10=sin(phi)/sin(phis)-1; %dn/n1
end
function [ kup ] = PoisKTPNBok(hv,r,phima,phimi,phi,teta)
ep=1e-6; vxyz=PovorotyIzIskhKNov([r*cos(phima)  r*sin(phima)   0],teta,phi); vxyz=PlusOdinPovorot(vxyz); k=1; phim(k)=phima; xym(k)=vxyz(2); 
phima=phima+pi; vxyz=PovorotyIzIskhKNov([r*cos(phima)  r*sin(phima)   0],teta,phi); vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phima; xym(k)=vxyz(2); 
vxyz=PovorotyIzIskhKNov([r*cos(phimi)  r*sin(phimi)   0],teta,phi); vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phimi; xym(k)=vxyz(2); 
phimi=phimi+pi; vxyz=PovorotyIzIskhKNov([r*cos(phimi)  r*sin(phimi)   0],teta,phi); vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phimi; xym(k)=vxyz(2); 
if (xym(1)>xym(2))
        minx=xym(2); phimix=phim(2); maxx=xym(1); phimax=phim(1);
else
        minx=xym(1); phimix=phim(1); maxx=xym(2); phimax=phim(2); 
end
if (xym(3)>xym(4))
        miny=xym(4); phimiy=phim(4); maxy=xym(3); phimay=phim(3); phimiy=(phimiy+phimay)/2;
else
        miny=xym(3); phimiy=phim(3); maxy=xym(4); phimay=phim(4); phimiy=(phimiy+phimay)/2;
end
f=1; q=1;
phimiy=180/pi*phimiy'
phimay=180/pi*phimay'
phimix=180/pi*phimix'
phimax=180/pi*phimax'
phimi=min([phimiy phimix])
phima=max([phimay phimax])
phix=0:1e-4:2*pi;
for k=1:length(phix)
    phimia(k)=phimi; phimaa(k)=phima;
ksitp=minx+(maxx-minx)*rand(); etatp=miny*(1-rand()); 
vxy=PovorotyIzNovKIskh([ksitp   etatp   0],teta,phi); x=vxy(1); y=vxy(2); alpha=atan(y/x); 
if (alpha<0)
alpha=alpha+pi; 
end
al(k)=alpha*180/pi;
if ((alpha<phimax) && (alpha>phimix) && (alpha<phimay) && (alpha>phimiy) && (f>0))
    q=q+1;
else
end
end
pl=plot(phix,phimia,':b',phix,phimaa,':k',phix,al,':m');
set(pl,'LineWidth',1); hold on; grid on; 
hold on; grid on; 
xlabel({'x'}); 
ylabel({'y'}); 
title({'y(x)'});
q=q'
    zk=hv*rand(); %попали по бок. пов-ти (в образующую)
    if (teta>pi/2)
    phitp=pi*(-1+2*rand())/2; %точка падения луча на верхнюю часть по отношению к OZ
    else
    phitp=pi*(1+2*rand())/2; %точка падения луча на нижнюю часть, OZ вверх
    end
    kup=[0  0   0];
end