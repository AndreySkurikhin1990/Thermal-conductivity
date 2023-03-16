%определение коэффициента теплопроводности вермикулита иными способами
function [ srlam ] = opredKoeTepSrInSp(vyb, po, tem, te0, laefm)
for k=1:length(tem)
ts=tem(k)-te0;
lavo=opredTeploprovVozd(ts);
laef=laefm(k);
switch (vyb)
    case 2
        lam = ReshNeyaUravMaxKoTepl(laef,lavo,po); %максимальный КТП - слои идут параллельно тепловому потоку
    case 3
        lam = ReshNeyaUravMinKoTepl(laef,lavo,po); %минимальный КТП - слои идут перпендикулярно тепловому потоку %не подходит
    case 4
        lam = ReshNeyaUravKoTeplLiechtenecker(laef, lavo, po); 
    case {5,6}
        la = opredRdeltaSrVzaPro(po,laef,lavo); 
        if (vyb<6)
            lam = la(1); 
        else lam = la(2); 
        end
    case 7
        lam = KoTeplFrickEllipsoidChasti(laef, lavo, po); 
    case 8
        lam = ReshNeyaUravKoTeplMaxwell(laef, lavo, po); %не подходит 
    case 9
        lam = ReshNeyaUravKoTeplBurgFrik(laef, lavo, po); %не подходит
    case 10
        lam = ReshNeyaUravKoTeplShumanVoss(laef, lavo, po); %не подходит
    case 11
        lam = ReshNeyaUravKoTeplBruggman(laef, lavo, po); %не подходит
    case 12
        lam = ReshNeyaUravKoTeplGoringChirchillSpher(laef, lavo, po); %не подходит
    case 13
        lam = ReshNeyaUravKoTeplGoringChirchillCyl(laef, lavo, po); %не подходит
    case 14
        lam = ReshNeyaUravKoTeplGoringChirchillCylPlotSl(laef, lavo); %плотный слой %не подходит
    case 15
        lam = ReshNeyaUravKoTeplMaxwAcken(laef, lavo, po); %не подходит
    case 16
        lam = ReshNeyaUravKoTeplSwift(laef, lavo); %не подходит
    case 17
        lam = MetodHashinaShtrikmanaMin(laef,lavo,po);
    case 18
        lam = MetodHashinaShtrikmanaMax(laef,lavo,po);
end
srla(k) = urovOtsechen(lam,lavo,po);
end
laefm=laefm';
f = 1;
    for k=1:length(tem)
        if ((srla(k) == 0) && (f > 0))
            f = 0;
        end
    end
    if (f==0)
        for k=1:length(srla)
            srla(k)=0;
        end
    end
srlam=srla;
end

function fl = urovOtsechen(pro,ref,por)
po=30:10:90; po(length(po)+1)=95; po=1e-2*po;
uom=[220,150,110,80,60,45,30,20];
urot=opredToch(uom,por,po,length(uom)); %определить точнее
if (pro<=0)
    f=0;
else
if (abs(pro)>=abs(urot*ref))
    f=0;
else f=1;
end
end
if (f>0)
fl=pro;
else
    fl=0;
end
end

function ot = opredToch(ktpts,te,tref,n)
f=1; p=1;
for k=1:n
    if ((tref(k)>te) && (f>0))
		p=k; f=0; 
    end
end
ktptsp=0; ko=0;
if (f==0)
if ((p>1) && (p<=n)) 
    ko=(ktpts(p)-ktpts(p-1))/(tref(p)-tref(p-1));
    ktptsp=ktpts(p-1)+ko*(te-tref(p-1));
end
if (p==1)
    ko=(ktpts(2)-ktpts(1))/(tref(2)-tref(1));
	ktptsp=ktpts(1)+ko*(te-tref(1)); 
end
end
if ((f>0) && (p==1))
    ko=(ktpts(n)-ktpts(n-1))/(tref(n)-tref(n-1));
	ktptsp=ktpts(n-1)+ko*(te-tref(n-1));
end
ot=ktptsp;
end

function res = ReshNeyaUravMaxKoTepl(lae, lan, po) %lan - КТП непрерывной фазы (воздух), ищем КТП твердого скелета
ladb=1e4; lada=1e-9; ra=1e2; ep=1e-7; h=0; la=0; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=lan*po+(1-po)*lada-lae;
fb=lan*po+(1-po)*ladb-lae;
fc=lan*po+(1-po)*ladc-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = ReshNeyaUravMinKoTepl(lae, lan, po)
ladb=1e5; lada=1e-9; ra=1e5; ep=1e-7; h=0; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=lan*lada/(po*lada+(1-po)*lan)-lae;
fb=lan*ladb/(po*ladb+(1-po)*lan)-lae;
fc=lan*ladc/(po*ladc+(1-po)*lan)-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = ReshNeyaUravKoTeplLiechtenecker(lae, lan, po)
ladb=1e4; lada=1e-9; ra=1e2; ep=1e-7; la=0; h=0; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=(lada^(1-po))*(lan^po)-lae;
fb=(ladb^(1-po))*(lan^po)-lae;
fc=(ladc^(1-po))*(lan^po)-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = ReshNeyaUravKoTeplMaxwell(lae, lan, po)
ladb=1e9; lada=1e-9; ra=1e9; ep=1e-7; h=0; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=lan*(lada+2*lan-2*(1-po)*(lan-lada))/(lada+2*lan+(1-po)*(lan-lada))-lae;
fb=lan*(ladb+2*lan-2*(1-po)*(lan-ladb))/(ladb+2*lan+(1-po)*(lan-ladb))-lae;
fc=lan*(ladc+2*lan-2*(1-po)*(lan-ladc))/(ladc+2*lan+(1-po)*(lan-ladc))-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = ReshNeyaUravKoTeplBurgFrik(lae, lan, po)
ladb=1e5; lada=1e-9; ra=1e2; ep=1e-7;
la=0; h=0; f=0; hko=1e4;
f(1)=1/8; f(2)=f(1); f(3)=1-f(1)-f(2); %де Вриз - про почвы
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;  
fba=0; fbb=0; fbc=0;
for k=1:3
     fba=fba+(1+((lada/lan-1)*f(k)))^(-1);
     fbb=fbb+(1+((ladb/lan-1)*f(k)))^(-1);
     fbc=fbc+(1+((ladc/lan-1)*f(k)))^(-1);
end
fba=fba/3; fbb=fbb/3; fbc=fbc/3;
fa=lan*(1+(1-po)*(fba*lada/lan-1))/(1+(1-po)*(fba-1))-lae;
fb=lan*(1+(1-po)*(fbb*ladb/lan-1))/(1+(1-po)*(fbb-1))-lae;
fc=lan*(1+(1-po)*(fbc*ladc/lan-1))/(1+(1-po)*(fbc-1))-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = ReshNeyaUravKoTeplShumanVoss(lae, lan, po)
ladb=1e5; lada=1e-9; ra=1e5; ep=1e-7; la=0; h=0; lama=0;
p=PoiskPShumanVoss(po); hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2; 
lamaa=lan*lada/(lan+p*(lan-lada));
lamaa=lamaa*(1+p*(1+p)*(lan-lada)*log(lan*(1+p)/p/lada)/(lan+p*(lan-lada)));
fa=lan*(po^3)+(1-po^3)*lamaa-lae;
lamab=lan*ladb/(lan+p*(lan-ladb));
lamab=lamab*(1+p*(1+p)*(lan-ladb)*log(lan*(1+p)/p/ladb)/(lan+p*(lan-ladb)));
fb=lan*(po^3)+(1-(po^3))*lamab-lae;
lamac=lan*ladc/(lan+p*(lan-ladc));
lamac=lamac*(1+p*(1+p)*(lan-ladc)*log(lan*(1+p)/p/ladc)/(lan+p*(lan-ladc)));
fc=lan*(po^3)+(1-po^3)*lamac-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function poisp = PoiskPShumanVoss(po)
pb=1e5; pa=1e-9; ra=1e5; ep=1e-7;
la=0; h=0; hko=1e4;
while ((ra>ep) && (h<hko))
pc=(pa+pb)/2;
fa=(pa^2+pa)*log(1+1/pa)-pa-po;
fb=(pb^2+pb)*log(1+1/pb)-pb-po;
fc=(pc^2+pc)*log(1+1/pc)-pc-po;
la = VybCAB(fa,fb,fc,pa,pb,pc);
pa=la(1); pb=la(2); pc=la(3);
ra=abs(fa-fb);
h=h+1;
end
poisp=pc;
end

function res = ReshNeyaUravKoTeplGoringChirchillSpher(lae, lan, po)
ladb=1e4; lada=1e-9; ra=1e4; ep=1e-7; h=0; vd=1-po; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2; 
fa=lae/lan-(2+lada/lan-2*vd*(1-lada/lan))/(2+lada/lan+vd*(1-lada/lan));
fb=lae/lan-(2+ladb/lan-2*vd*(1-ladb/lan))/(2+ladb/lan+vd*(1-ladb/lan));
fc=lae/lan-(2+ladc/lan-2*vd*(1-ladc/lan))/(2+ladc/lan+vd*(1-ladc/lan));
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = ReshNeyaUravKoTeplGoringChirchillCyl(lae, lan, po)
ladb=1e4; lada=1e-9; ra=1e2; ep=1e-7; la=0; h=0; vd=1-po; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2; 
fa=((1+lada/lan)/(1-lada/lan)-vd)/((1+lada/lan)/(1-lada/lan)+vd)-lae/lan;
fb=((1+ladb/lan)/(1-ladb/lan)-vd)/((1+ladb/lan)/(1-ladb/lan)+vd)-lae/lan;
fc=((1+ladc/lan)/(1-ladc/lan)-vd)/((1+ladc/lan)/(1-ladc/lan)+vd)-lae/lan;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = ReshNeyaUravKoTeplBruggman(lae, lan, po)
%Модель Бруггмана-Ханаи - только для круглых частиц
ladb=1e4; lada=1e-9; ra=1e4; ep=1e-7; la=0; h=0; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2; 
fa=po-(lada-lae)*((lan/lae)^(1/3))/(lada-lan);
fb=po-(ladb-lae)*((lan/lae)^(1/3))/(ladb-lan);
fc=po-(ladc-lae)*((lan/lae)^(1/3))/(ladc-lan);
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = ReshNeyaUravKoTeplMaxwAcken(lae, lan, po)
ladb=1e5; lada=1e-9; ra=1e5; ep=1e-7; h=0; vd=1-po; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;
fa=(1-lan/lada)/(2*lan/lada+1);
fa=(1+2*vd*fa)/(1-vd*fa)-lae/lan;
fb=(1-lan/ladb)/(2*lan/ladb+1);
fb=(1+2*vd*fb)/(1-vd*fb)-lae/lan;
fc=(1-lan/ladc)/(2*lan/ladc+1);
fc=(1+2*vd*fc)/(1-vd*fc)-lae/lan;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function res = KoTeplFrickEllipsoidChasti(lae, lan, po)
ladb=1e4; lada=1e-9; ra=1e4; ep=1e-7; la=0; h=0; s=0;
AsRa=99e-2;
%Для пор в виде эллипсоидов
%b=2; a=4; c=1;
%alp=acos(c/a); 
%p=sqrt((a^2-b^2)/(a^2-c^2)); m=p^2;  %0<=m<=1
%в Матлабе модуль и параметр m связаны: m=k^2
%labc(1)=(ellipticF(alp,m)-ellipticE(alp,m))*2/(a^2-b^2)/sqrt(a^2-c^2);
%alp=asin(c/b);
%[peiK,peiE]=ellipke(m); %полный эллипт. интеграл 2 рода, Фрике, 1953
%labc(3)=2*a/(a^2-c^2)/b/c-2/(b^2-c^2)/sqrt(a^2-c^2)*(peiE-ellipticE(alp,k));
%labc(2)=2/(a*b*c)-(labc(1)+labc(2));
%требуется вычисление неполных эллиптических интегралов
%labc(2)=ellipticE(alp,m)*2*sqrt(a^2-c^2)/(a^2-b^2)/(b^2-c^2)-...
%2/(a^2-b^2)/sqrt(a^2-c^2)*ellipticF(alp,m)-2*c/(b^2-c^2)/a/b;
%labc(3)=ellipticE(alp,m)*2/(c^2-b^2)/sqrt(a^2-c^2)+2/(b^2-c^2)*b/a/c;
%b=2; a=AsRa*b; c=b; M=opredMEllipIntegSpher(a,b); M=M*a*(b^2); la=2/a/b^2-2*M %Фрике, 1924
%Для пор в виде сфероидов - симметрия вращения относительно одной оси
b=2; a=AsRa*b; c=b;
if (a>b)
e=a^2-b^2; p0=sqrt(e)/a; M=-2*p0+log(abs((1+p0)/(1-p0))); M=M/(e^(3/2));
else
e=b^2-a^2; p0=sqrt(e)/a; M=-p0+atan(p0); M=M*2/(e^(3/2));
end
%lb=a/abs(a^2-b^2)/b^2+log(abs(a-sqrt(abs(a^2-b^2)))/abs(a+sqrt(abs(a^2-b^2))))/2/(abs(a^2-b^2))^(3/2)
%la=-2/a/abs(a^2-b^2)-log(abs(a-sqrt(abs(a^2-b^2)))/abs(a+sqrt(abs(a^2-b^2))))/(abs(a^2-b^2))^(3/2)
lb=2/a/b/c-M; lb=lb/2;
labc=[M,lb,lb];
for k=1:3
    xa(k)=2/a/b/c/labc(k)-1;
end
%xa=xa'
nko=1e3;
while ((ra>ep) && (h<nko))
ladc=(lada+ladb)/2;
sa=0; sb=0; sc=0;
for k=1:3
    pea(k)=(1+xa(k))/(xa(k)+lan/lada);
    sa=sa+pea(k);
    peb(k)=(1+xa(k))/(xa(k)+lan/ladb);
    sb=sb+peb(k);
    pec(k)=(1+xa(k))/(xa(k)+lan/ladc);
    sc=sc+pec(k);
end
fa=lada+po*(lan-lada)*sa/3-lae; %объемная доля растворенного вещества (эллипсоидов)
fb=ladb+po*(lan-ladb)*sb/3-lae;
fc=ladc+po*(lan-ladc)*sc/3-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end

function repl = opredMEllipIntegSpher(a,b)
if (a<b)
    phi=acos(a/b);
    M=phi-sin(2*phi)/2;
    M=M*cos(phi)/((sin(phi))^3);
else
    phis=acos(b/a);
    M=(1+sin(phis))/(1-sin(phis)); M=abs(M);
    M=((cos(phis))^2)/((sin(phis))^3)*log(M)/2;
    M=1/((sin(phis))^2)-M;
end
repl=M;
end

function [ vy ] = VybCAB(fa,fb,fc,xa,xb,xc)
if (((fa*fc)<0) && ((fb*fc)>0))
    xb=xc; 
end
if (((fa*fc)>0) && ((fb*fc)<0))
    xa=xc; 
end
vy = [xa,xb,xc];
end

function  fb = opredFiBol(delb,a,lm)
nb=1e2;
r=2*delb/sqrt(pi);
b=2*a/sqrt(pi);
bb=b/lm;
rr=r/lm;
sn=0;
for k=0:nb
    n=2*k+1;
    I1b=besseli(1,n*pi*bb);
    I1r=besseli(1,n*pi*rr);
    K1b=besselk(1,n*pi*bb);
    K1r=besselk(1,n*pi*rr);
   sn=sn+I1b/n^2/I1r*(I1r*K1b-K1r*I1b); 
end
fb=1-16*sn/pi^2;
end

%взаимопроникающие компоненты
function [ rde ] = opredRdeltaSrVzaPro(por,laef,lavo)
L=1e0;
Pa=1e1;
Pdelt=1e2;
delb=PoiskDelta(por,L);
a=delb/Pa;
la1=raschTeplVzaimoPronKomp(delb,L,laef,lavo); %поиск КТП  твердого каркаса - грубая оценка
delm=L/Pdelt; %дельта малая - толщина трещин
lm=2*L-delm; %эль малая
Fb=opredFiBol(delb,a,lm);
la1u=utochLamTverKark(delm,a,Fb,lm,la1,delb,lavo,laef,L); %уточнение КТП твердого каркаса - тонкая оценка
rde=[la1,la1u];
end

%поиск дельты большой
function pd = PoiskDelta(m2,L)
ca=0; cb=1e0; ra=abs(ca-cb);
ep=1e-7; h=0; nit=1e2;
while ((ra>ep) && (h<nit))
cc=(ca+cb)/2;
fa=2*(ca^3)-3*(ca^2)+1-m2; 
fb=2*(cb^3)-3*(cb^2)+1-m2; 
fc=2*(cc^3)-3*(cc^2)+1-m2; 
pode = VybCAB(fa,fb,fc,ca,cb,cc);
ca=pode(1); cb=pode(2); cc=pode(3);
ra=abs(ca-cb);
h=h+1;
end
pd=cc*L; %cc - относительный размер бруса
end

%взаимпопроникающие компоненты
function vpk = raschTeplVzaimoPronKomp(del,L,laef,la_v)
R=1/laef/L; R3=1/la_v/del; R4=L/la_v/(L-del)^2; 
la1b=1e4; la1a=1e-9; ra=1e4; ep=1e-7; h=0; hko=1e4;
while ((ra>ep) && (h<hko))
la1c=(la1a+la1b)/2;    
R1a=L/la1a/del^2;
R1b=L/la1b/del^2;
R1c=L/la1c/del^2;
R2a=1/la1a/(L-del);
R2b=1/la1b/(L-del);
R2c=1/la1c/(L-del);
fa=1/R-1/R1a-2/(R2a+R3)-1/R4;
fb=1/R-1/R1b-2/(R2b+R3)-1/R4;
fc=1/R-1/R1c-2/(R2c+R3)-1/R4;
lam1 = VybCAB(fa,fb,fc,la1a,la1b,la1c);
la1a=lam1(1); la1b=lam1(2); la1c=lam1(3);
ra=abs(fa-fb);
h=h+1;
end
vpk=la1c;
end

%структура с переменным сечением компонент
function lu = utochLamTverKark(delm,a,phib,lm,la1,delb,lavo,laef,L)
R=1/laef/L; b=2*a/sqrt(pi); delsh=delm/L;
r=2*delb/sqrt(pi); la2=lavo; y=b/r;
laz=lavo; R6=delm/2/laz/((delb^2)-(a^2)); %laz - КТП компоненты, заполняющей трещины (воздух)
R4=L/la2/((L-delb)^2); R3=1/delb/la2;
phim=phib*(y^2)/(phib-(y^2));
la1a=1e-9; la1b=la1+1e3;
ra=1e5; ep=1e-7; h=0; nit=1e2;
while ((ra>ep) && (h<nit))
la1c=(la1a+la1b)/2;
lamsha=la1a; %термосопротивление шейки
R1a=delm/2/a^2/lamsha; 
Rp0a=lm/la1a/pi/(b^2); %полное сопротивление квадратного бруса
R7a=2*Rp0a*phim;
lamshb=la1b; %термосопротивление шейки
R1b=delm/2/a^2/lamshb;
Rp0b=lm/la1b/pi/(b^2); %полное сопротивление квадратного бруса
R7b=2*Rp0b*phim;
lamshc=la1c; %термосопротивление шейки
R1c=delm/2/a^2/lamshc;
Rp0c=lm/la1c/pi/(b^2); %полное сопротивление квадратного бруса
R7c=2*Rp0c*phim;
R2a=2*Rp0a*phib; %тепловое сопротивление восьмой части центрального бруса
R2b=2*Rp0b*phib;
R2c=2*Rp0c*phib;
R5a=1/(L-delb)/la1a;
R5b=1/(L-delb)/la1b;
R5c=1/(L-delb)/la1c;
Rdsa=1/(R1a+R2a)+1/(R6+R7a); Rdsa=1/Rdsa;
Rdsb=1/(R1b+R2b)+1/(R6+R7b); Rdsb=1/Rdsb;
Rdsc=1/(R1c+R2c)+1/(R6+R7c); Rdsc=1/Rdsc;
Rdelmina=RdeltaMin(R1a,R6,la1a,delb,L,delsh);
Rdelmaxa=RdeltaMax(la1a,a,delb,L,delsh,R1a,R6);
Rdelminb=RdeltaMin(R1b,R6,la1b,delb,L,delsh);
Rdelmaxb=RdeltaMax(la1b,a,delb,L,delsh,R1b,R6);
Rdelminc=RdeltaMin(R1c,R6,la1c,delb,L,delsh);
Rdelmaxc=RdeltaMax(la1c,a,delb,L,delsh,R1c,R6);
if ((Rdsa>=Rdelmina) && (Rdsa<=Rdelmaxa))
fa=1/Rdsa+1/R4+2/(R3+R5a)-1/R;
end
fa=1/Rdsa+1/R4+2/(R3+R5a)-1/R;
if ((Rdsb>=Rdelminb) && (Rdsb<=Rdelmaxb))
fb=1/Rdsb+1/R4+2/(R3+R5b)-1/R;
end
fb=1/Rdsb+1/R4+2/(R3+R5b)-1/R;
if ((Rdsc>=Rdelminc) && (Rdsc<=Rdelmaxc))
fc=1/Rdsc+1/R4+2/(R3+R5c)-1/R;    
end
fc=1/Rdsc+1/R4+2/(R3+R5c)-1/R; 
lam1 = VybCAB(fa,fb,fc,la1a,la1b,la1c);
la1a=lam1(1); la1b=lam1(2); la1c=lam1(3);
ra=abs(fa-fb);
h=h+1;
end
lu=la1c;
end

function rdm = RdeltaMin(R1,R6,la1,delb,L,delsh)
rdma=(1/R1+1/R6);
rdm=1/rdma+(L-delsh/2)/la1/delb^2;
end

function rdm = RdeltaMax(la1,a,delb,L,dsh,R1,R6)
R2=L*(1-dsh/2)/la1/a^2;
R7=L*(1-dsh/2)/la1/(delb^2-a^2);
rdma=1/(R1+R2)+1/(R6+R7);
rdm=1/rdma;
end

function ps = ReshNeyaUravKoTeplGoringChirchillCylPlotSl(laef, lavo)
x0=(2+0.7)/2/2; x0=x0*1e-3; len=1;
%len=4;
%for k=1:len
%kCar(k)=10^(-4+(k-1)*2);
%end
lan=lavo;
for k=1:len
%kC=kCar(k)
kC=1e0; ladb=1e5; lada=1e-8; ra=1e5;
ep=1e-7; h=1; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
kBa=(lada/kC/(lan-lada))^(1/3);
kBb=(ladb/kC/(lan-ladb))^(1/3);
kBc=(ladc/kC/(lan-ladc))^(1/3);
CBa=log(sqrt(kBa^2-kBa*x0+x0^2)/(kBa+x0));
CBa=CBa+sqrt(3)*atan((2*x0-kBa)/sqrt(3*kBa));
CBa=CBa-sqrt(3)*atan(1/sqrt(3));
CBa=CBa*pi*lan*lada/6/(lan-lada)/kC/kBa;
CBa=CBa+lada*(1-pi*(x0^2)/4);
CBb=log(sqrt(kBb^2-kBb*x0+x0^2)/(kBb+x0));
CBb=CBb+sqrt(3)*atan((2*x0-kBb)/sqrt(3*kBb));
CBb=CBb-sqrt(3)*atan(1/sqrt(3));
CBb=CBb*pi*lan*ladb/6/(lan-ladb)/kC/kBb;
CBb=CBb+ladb*(1-pi*(x0^2)/4);
CBc=log(sqrt(kBc^2-kBc*x0+x0^2)/(kBc+x0));
CBc=CBc+sqrt(3)*atan((2*x0-kBc)/sqrt(3*kBc));
CBc=CBc-sqrt(3)*atan(1/sqrt(3));
CBc=CBc*pi*lan*ladc/6/(lan-ladc)/kC/kBc;
CBc=CBc+ladc*(1-pi*(x0^2)/4);
fa=CBa-laef;
fb=CBb-laef;
fc=CBc-laef;
lam1 = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=lam1(1); ladb=lam1(2); ladc=lam1(3);
ra=abs(fa-fb);
h=h+1;
end
ladc=ladc';
end
ps=ladc;
end
%оценка вилки - для изотропных в среднем образцов
function mini = MetodHashinaShtrikmanaMin(laef,lan,por)
lada=1e-8; ladb=1e6; ra=1e6; ep=1e-7; h=1; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
ta=1/(lada-lan)+(1-por)/3/lada; ta=lan+por/ta;
tb=1/(ladb-lan)+(1-por)/3/ladb; tb=lan+por/tb;
tc=1/(ladc-lan)+(1-por)/3/ladc; tc=lan+por/tc;
fa=ta-laef;
fb=tb-laef;
fc=tc-laef;
lam1 = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=lam1(1); ladb=lam1(2); ladc=lam1(3);
ra=abs(fa-fb);
h=h+1;
end
mini=ladc;
end

function maxi = MetodHashinaShtrikmanaMax(laef,lan,por)
lada=1e-8; ladb=1e6; ra=1e6; ep=1e-7; h=1; hko=1e4;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
ta=1/(lada-lan)+por/3/lada; ta=lada+(1-por)/ta;
tb=1/(ladb-lan)+por/3/ladb; tb=ladb+(1-por)/tb;
tc=1/(ladc-lan)+por/3/ladc; tc=ladc+(1-por)/tc;
fa=ta-laef;
fb=tb-laef;
fc=tc-laef;
lam1 = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=lam1(1); ladb=lam1(2); ladc=lam1(3);
ra=abs(fa-fb);
h=h+1;
end
maxi=ladc;
end

function res = ReshNeyaUravKoTeplSwift(lae, lan)
ladb=1e5; lada=1e-9; ra=1e5; ep=1e-7; h=0; la=0; hko=1e5;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=(log(abs(lan/lada)))*((lan/lada-1)^(-2));
fa=((lan/lada-1)^(-1))-fa;
fa=0.577*pi*fa+0.093-lae/lan;
fb=(log(abs(lan/ladb)))*((lan/ladb-1)^(-2));
fb=((lan/ladb-1)^(-1))-fb;
fb=0.577*pi*fb+0.093-lae/lan;
fc=(log(abs(lan/ladc)))*((lan/ladc-1)^(-2));
fc=((lan/ladc-1)^(-1))-fc;
fc=0.577*pi*fc+0.093-lae/lan;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
res=ladc;
end