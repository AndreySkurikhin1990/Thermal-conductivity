%определение КТП твердой части шамота иными способами
function [ srelam ] = opredTvChaKoeTepSrInSpSha(vyb, po, tem, laefm, x0, lavoz) %пористость шамота
%format long g; srla=0; ep=1e-8; pm=[5,12.5,17.5,22.5]*1e-2; am=[1.5,2,2.4,2.6]; pap=polyfit(pm,am,2); apo=polyval(pap,po)*po; urso=1/(1-apo); 
up=1; ep=1e-8; uo=urovPod(po);
lamax=1e5;
lamin=1e-9;
hko=1e5;
lear=length(tem);
for k=1:lear
laef=laefm(k);
lavo=lavoz(k);
switch (vyb)
    case (1)
        lam = ReshNeyaUravMaxKoTepl(laef,lavo,po,ep,lamax,lamin,hko); %максимальный КТП - слои идут параллельно тепловому потоку
    case (2)
        lam = ReshNeyaUravMinKoTepl(laef,lavo,po,ep,lamax,lamin,hko); %минимальный КТП - слои идут перпендикулярно тепловому потоку
    case (3)
        lam = ReshNeyaUravKoTeplLiechtenecker(laef, lavo, po, ep,lamax,lamin,hko);
    case (4)
        %disp('4');
        lam = ReshNeyaUravKoTeplMaxwell(laef, lavo, po, ep,lamax,lamin,hko); %не подходит 
    case (5)
        lam = ReshNeyaUravKoTeplBurgFrik(laef, lavo, po, ep,lamax,lamin,hko); %не подходит
    case (6)
        %disp('6');
        lam = ReshNeyaUravKoTeplShumanVoss(laef, lavo, po, ep,lamax,lamin,hko); %не подходит
    case (7)
        lam = ReshNeyaUravKoTeplGoringChirchillSpher(laef, lavo, po, ep,lamax,lamin,hko); %не подходит
    case (8)
        lam = ReshNeyaUravKoTeplGoringChirchillCyl(laef, lavo, po, ep,lamax,lamin,hko); %не подходит
    case (9)
        lam = ReshNeyaUravKoTeplBruggman(laef, lavo, po, ep,lamax,lamin,hko); %не подходит
    case (10)
        lam = ReshNeyaUravKoTeplMaxwAcken(laef, lavo, po, ep,lamax,lamin,hko); %не подходит
    case (11)
        lam = KoTeplFrickEllipsoidChasti(laef, lavo, po, ep,lamax,lamin,hko); 
    case (12)
        lam = opredRdeltaSrVzaPro(po,laef,lavo,ep,1); 
    case (13)
        %disp('13');
        lam = opredRdeltaSrVzaPro(po,laef,lavo,ep,2);
    case (14)
        %disp('14');
        lam = ReshNeyaUravKoTeplGoringChirchillCylPlotSl(laef, lavo, ep, x0,lamax,lamin,hko); %плотный слой %не подходит
    case (15)
        lam=MetodHashinaShtrikmanaMin(laef,lavo,po, ep,lamax,lamin,hko);
    case (16)
        lam=MetodHashinaShtrikmanaMax(laef,lavo,po, ep,lamax,lamin,hko);
    case (17)
        %disp('17');
        lam = ReshNeyaUravKoTeplSwift(laef, lavo, ep,lamax,lamin,hko); %не подходит
    case (18)
        lam=SysPlastShakhPor(laef, lavo, po, ep,lamax,lamin,hko);
    case (19)
        lam=MetodStarostina(laef, lavo, po, ep,lamax,lamin,hko);
    case (20) 
        lam=MetodRusselNeprTverTelo(laef, lavo, po, ep,lamax,lamin,hko);
    case (21) 
        lam=MetodRusselNeprVozd(laef, lavo, po, ep,lamax,lamin,hko);
    otherwise
        disp('Incorrect number');
end
srla(k)=lam;
end
stsrla=srla;
srla=proverkakvi(srla, laefm, lavoz, up, uo, lear, ep);
switch vyb
    case 4 
        %disp('4');
        srla=srla';
    case 6
        %disp('6');
        stsrla=stsrla';
        srla=srla';
    case 13
        %disp('13');
        stsrla=stsrla';
        srla=srla';
    case 14
        %disp('14');
        srla=srla';
    case 17
        %disp('17');
        srla=srla';
end
srelam=srla;
end

function fli = urovPodder(pro,ref,urpo)
f=1; 
if (pro<0) 
    f=-1; fl=0;
else
    if (pro<(urpo*ref))
    f=-1; 
    else
        f=1;
    end
end
if (f>0) 
    fl=pro;
else
    fl=0;
end
fli=fl;
end

function fli = urovOtsech(pro,ref,urot)
fl=0; f=1; 
if (pro<=0) 
    pro=0; fl=0;
else
    if (abs(pro)>=abs(urot*ref))
    f=0; 
    else
        f=1;
    end
end
if (f>0) 
    fl=pro;
end
fli=fl;
end

function res = ReshNeyaUravMaxKoTepl(lae, lan, po, ep,ladb,lada,hko) %lan - КТП непрерывной фазы (воздух), ищем КТП твердого скелета
ra=1e2; h=0;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=lan*po+(1-po)*lada-lae;
fb=lan*po+(1-po)*ladb-lae;
fc=lan*po+(1-po)*ladc-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(lada-ladb);
h=h+1;
end
lae=lae'; lan=lan'; ladc=ladc';
res=ladc;
end

function res = ReshNeyaUravMinKoTepl(lae, lan, po, ep,ladb,lada,hko)
ra=1e5; h=0;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=lan*lada/(po*lada+(1-po)*lan)-lae;
fb=lan*ladb/(po*ladb+(1-po)*lan)-lae;
fc=lan*ladc/(po*ladc+(1-po)*lan)-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(lada-ladb);
h=h+1;
end
%disp('2');
ladc=ladc';
res=ladc;
end

function res = ReshNeyaUravKoTeplLiechtenecker(lae, lan, po, ep,ladb,lada,hko)
ra=1e2; h=0;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=(lada^(1-po))*(lan^po)-lae;
fb=(ladb^(1-po))*(lan^po)-lae;
fc=(ladc^(1-po))*(lan^po)-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(lada-ladb);
h=h+1;
end
%disp('3');
ladc=ladc';
res=ladc;
end

function res = ReshNeyaUravKoTeplMaxwell(lae, lan, po, ep,ladb,lada,hko)
ra=1e9; h=0;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2;    
fa=lan*(lada+2*lan-2*(1-po)*(lan-lada))/(lada+2*lan+(1-po)*(lan-lada))-lae;
fb=lan*(ladb+2*lan-2*(1-po)*(lan-ladb))/(ladb+2*lan+(1-po)*(lan-ladb))-lae;
fc=lan*(ladc+2*lan-2*(1-po)*(lan-ladc))/(ladc+2*lan+(1-po)*(lan-ladc))-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(lada-ladb);
h=h+1;
end
%disp('4');
ladc=ladc';
res=ladc;
end

function res = ReshNeyaUravKoTeplBurgFrik(lae, lan, po, ep,ladb,lada,hko)
ra=1e2; h=0; f=0;
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
ra=abs(lada-ladb);
h=h+1;
end
%disp('5');
ladc=ladc';
res=ladc;
end

function res = ReshNeyaUravKoTeplShumanVoss(lae, lan, po, ep,ladb,lada,hko)
ra=1e5; h=0;
p=PoiskPShumanVoss(po, ep,ladb,lada,hko);
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
ra=abs(lada-ladb);
h=h+1;
end
%disp('6');
ladc=ladc';
res=ladc;
end

function poisp = PoiskPShumanVoss(po, ep,pb,pa,hko)
ra=1e5; h=0;
while ((ra>ep) && (h<hko))
pc=(pa+pb)/2;
fa=(pa^2+pa)*log(1+1/pa)-pa-po;
fb=(pb^2+pb)*log(1+1/pb)-pb-po;
fc=(pc^2+pc)*log(1+1/pc)-pc-po;
la = VybCAB(fa,fb,fc,pa,pb,pc);
pa=la(1); pb=la(2); pc=la(3);
ra=abs(pa-pb);
h=h+1;
end
poisp=pc;
end

function res = ReshNeyaUravKoTeplGoringChirchillSpher(lae, lan, po, ep,ladb,lada,hko)
ra=1e4; h=0; vd=1-po;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2; 
fa=lae/lan-(2+lada/lan-2*vd*(1-lada/lan))/(2+lada/lan+vd*(1-lada/lan));
fb=lae/lan-(2+ladb/lan-2*vd*(1-ladb/lan))/(2+ladb/lan+vd*(1-ladb/lan));
fc=lae/lan-(2+ladc/lan-2*vd*(1-ladc/lan))/(2+ladc/lan+vd*(1-ladc/lan));
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(lada-ladb);
h=h+1;
end
%disp('7');
ladc=ladc';
res=ladc;
end

function res = ReshNeyaUravKoTeplGoringChirchillCyl(lae, lan, po, ep,ladb,lada,hko)
ra=1e2; h=0; vd=1-po;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2; 
fa=((1+lada/lan)/(1-lada/lan)-vd)/((1+lada/lan)/(1-lada/lan)+vd)-lae/lan;
fb=((1+ladb/lan)/(1-ladb/lan)-vd)/((1+ladb/lan)/(1-ladb/lan)+vd)-lae/lan;
fc=((1+ladc/lan)/(1-ladc/lan)-vd)/((1+ladc/lan)/(1-ladc/lan)+vd)-lae/lan;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(lada-ladb);
h=h+1;
end
%disp('8');
ladc=ladc';
res=ladc;
end

function res = ReshNeyaUravKoTeplBruggman(lae, lan, po, ep,ladb,lada,hko) %Модель Бруггмана-Ханаи - только для круглых частиц
ra=1e4; h=0;
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2; 
fa=po-(lada-lae)*((lan/lae)^(1/3))/(lada-lan);
fb=po-(ladb-lae)*((lan/lae)^(1/3))/(ladb-lan);
fc=po-(ladc-lae)*((lan/lae)^(1/3))/(ladc-lan);
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(lada-ladb);
h=h+1;
end
%disp('9');
ladc=ladc';
res=ladc;
end

function res = ReshNeyaUravKoTeplMaxwAcken(lae, lan, po, ep,ladb,lada,hko)
ra=1e5; h=0; vd=1-po;
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
ra=abs(lada-ladb);
h=h+1;
end
%disp('10');
ladc=ladc';
res=ladc;
end

function res = KoTeplFrickEllipsoidChasti(lae, lan, po, ep,ladb,lada,hko)
ra=1e4; h=0; s=0;
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
while ((ra>ep) && (h<hko))
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
ra=abs(lada-ladb);
h=h+1;
end
%disp('11');
ladc=ladc';
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

function fb= opredFiBolN(x,y)
n=11;
y1=[0,.047,.113,.195,.3,.413,.563,.725,.832,.932,1]; 
y01=[0,.261,.411,.526,.618,.7,.774,.837,.9,.958,1];
fib=0:1e-1:1; f1=1; f2=1; q1=1; q2=1;
for k=1:n
    if ((f1>0) && (y<y01(k))) 
        q1=k; f1=0; break;
    end
	if ((f2>0) && (y<y1(k))) 
        q2=k; f2=0; break;
    end
end
if (q1==1) 
    q1=2; 
end
if (q2==1) 
    q2=2;
end
fb01=fib(q1-1)+(fib(q1)-fib(q1-1))*(y-y01(q1-1))/(y01(q1)-y01(q1-1));
fb1=fib(q2-1)+(fib(q2)-fib(q2-1))*(y-y1(q2-1))/(y1(q2)-y1(q2-1));
fb=(fb1+(fb01-fb1)*(x-1e0)/(1e-1-1e0));
end

function fi = opredFiBolNN(nu,m2)
cb=1e2; ca=-1e1; ra=1e2; ep=1e-9; h=0; kit=100;
while ((ra>ep) && (h<kit))
cc=(ca+cb)/2;
fa=2*(ca^3)-3*(ca^2)+1-m2;
fb=2*(cb^3)-3*(cb^2)+1-m2;
fc=2*(cc^3)-3*(cc^2)+1-m2;
if ((fc*fb>0) && (fa*fc<0)) 
    cb=cc; 
end
    if ((fc*fa>0) && (fb*fc<0)) 
        ca=cc; 
    end
ra=abs(ca-cb); h=h+1;
end
fi1=abs(cc/(2-cc)); fi1=sqrt(fi1);
fi=(nu^2)*(7.5-11*nu+4.5*(nu^2))*(1-fi1)+fi1;
end
%взаимопроникающие компоненты
function rde = opredRdeltaSrVzaPro(por,laef,lavo,ep,vyvm)
L=1e0;
Pa=1e1;
Pdelt=1e2;
delb=PoiskDelta(por,L,ep);
a=delb/Pa;
la1=raschTeplVzaimoPronKomp(delb,L,laef,lavo,ep); %поиск КТП твердого каркаса - грубая оценка
delm=L/Pdelt; %дельта малая - толщина трещин
lm=2*L-delm; %эль малая
%disp('Func Bess');
Fb=opredFiBol(delb,a,lm);
%disp('Prib po graf');
x=2*delb/sqrt(pi)/lm;
y=a/delb;
fbn=opredFiBolN(x,y);
if (la1>ep) 
    fbnn=lavo/la1; 
    fbnn=opredFiBolNN(fbnn, por); 
else fbnn=0;
end
epf=1e-2;
if (abs(fbn-fbnn)<epf) 
    Fb=(fbn+fbnn)/2; 
else Fb=fbn; %lam2 - пора, lam1 - твердый каркас
end
%disp('Prib po poly');
fbnn=fbnn';
Fb=Fb';
la1u=utochLamTverKark(delm,a,Fb,lm,la1,delb,lavo,laef,L,ep); %уточнение КТП твердого каркаса - тонкая оценка
if (vyvm==1)
%disp('12');
rde=la1;
elseif (vyvm==2)
%disp('13');
rde=la1u;
else
    rde=0;
end
end
%поиск дельты большой
function pd = PoiskDelta(m2,L,ep)
ca=0; cb=1e0; ra=abs(ca-cb);
h=0; nit=100;
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
function vpk = raschTeplVzaimoPronKomp(del,L,laef,la_v,ep)
R=1/laef/L; R3=1/la_v/del; R4=L/la_v/(L-del)^2; 
la1b=1e4; la1a=1e-9; ra=1e4; h=0; hko=100;
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
ra=abs(la1a-la1b);
h=h+1;
end
vpk=la1c;
end
%структура с переменным сечением компонент
function lu = utochLamTverKark(delm,a,phib,lm,la1,delb,lavo,laef,L,ep)
R=1/laef/L; b=2*a/sqrt(pi); delsh=delm/L;
r=2*delb/sqrt(pi); la2=lavo; y=b/r;
laz=lavo; R6=delm/2/laz/((delb^2)-(a^2)); %laz - КТП компоненты, заполняющей трещины (воздух)
R4=L/la2/((L-delb)^2); R3=1/delb/la2;
phim=phib*(y^2)/(phib-(y^2));
la1a=1e-9; la1b=la1+1e5;
ra=1e5; h=0; nit=100;
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
ra=abs(la1a-la1b);
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

function srp = SreRazPorSha()
rpr = 0;
raspr=[21.8,21.75,21.625,21.25,20.75,20.25,19.75,19,18.25,17.5,12,9.875,8.75,8.25,7.25,6.25,5,4.25,3.25,0.875,0.325,0];
for k=2:length(raspr)
    rpr(k)=(raspr(k-1)-raspr(k))*1e2/raspr(1);
end
prgr=[0.05,0.1,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50];
legr(1)=0;
for k=2:length(prgr)
    legr(k)=prgr(k-1);
end
sr = (legr + prgr) / 2; %средний размер пор в микронах
le = length(sr);
srpo = 0;
for k=1:le
    if (rpr(k) > 0)
    srpo = srpo + rpr(k) * sr(k) / 1e2; %объемная доля поры заданного размера в полном (во всем) объеме, в процентах
    end
end
srp = srpo*1e-6;
end

function ps = ReshNeyaUravKoTeplGoringChirchillCylPlotSl(laef, lavo, ep, x0,ladb,lada,hko)
%len=1;
%len=4;
%for k=1:len
%kCar(k)=10^(-4+(k-1)*2);
%end
laef=laef';
lan=lavo';
%for k=1:len
%kC=kCar(k)
kC=1e0; %kC=1e2; 
ra=1e5; h=0; 
while ((ra>ep) && (h<hko))
ladc=(lada+ladb)/2; 
kBa=lada/kC/(lan-lada);
kBa=sign(kBa)*(abs(kBa))^(1/3); 
%kBa=real(kBa); 
kBa=abs(kBa);
kBb=ladb/kC/(lan-ladb);
kBb=sign(kBb)*(abs(kBb))^(1/3);
if (h<10)
    kBb=kBb';
end
%kBb=real(kBb);
kBb=abs(kBb);
kBc=ladc/kC/(lan-ladc);
kBc=sign(kBc)*(abs(kBc))^(1/3); 
%kBc=real(kBc); 
kBc=abs(kBc);
CBa=kBa^2-kBa*x0+x0^2; 
CBa=abs(CBa);
CBa=log(sqrt(CBa)/abs(kBa+x0));
CBa=CBa+sqrt(3)*atan((2*x0-kBa)/sqrt(3*kBa));
CBa=CBa-sqrt(3)*atan(1/sqrt(3));
CBa=CBa*pi*lan*lada/6/(lan-lada)/kC/kBa;
CBa=CBa+lada*(1-pi*(x0^2)/4);
CBb=kBb^2-kBb*x0+x0^2;
CBb=abs(CBb);
CBb=log(sqrt(CBb)/abs(kBb+x0));
CBb=CBb+sqrt(3)*atan((2*x0-kBb)/sqrt(3*kBb));
CBb=CBb-sqrt(3)*atan(1/sqrt(3));
CBb=CBb*pi*lan*ladb/6/(lan-ladb)/kC/kBb;
CBb=CBb+ladb*(1-pi*(x0^2)/4);
CBc=kBc^2-kBc*x0+x0^2;
CBc=abs(CBc);
CBc=log(sqrt(CBc)/abs(kBc+x0));
CBc=CBc+sqrt(3)*atan((2*x0-kBc)/sqrt(3*kBc));
CBc=CBc-sqrt(3)*atan(1/sqrt(3));
CBc=CBc*pi*lan*ladc/6/(lan-ladc)/kC/kBc;
CBc=CBc+ladc*(1-pi*(x0^2)/4);
fa=CBa-laef;
fb=CBb-laef;
fc=CBc-laef;
%lam1 = VybCAB(fa,fb,fc,lada,ladb,ladc);
%lada=lam1(1); ladb=lam1(2); ladc=lam1(3);
if ((fc*fb>0) && (fa*fc<0)) 
    ladb=ladc; 
end
    if ((fc*fa>0) && (fb*fc<0)) 
        lada=ladc; 
    end
ra=abs(lada-ladb);
h=h+1;
end
%disp('14');
ladc=ladc';
%end
ps=ladc;
end
%оценка вилки - для изотропных в среднем образцов
function mini = MetodHashinaShtrikmanaMin(laef,lan,por,ep,ladb,lada,hko)
ra=1e6; h=0;
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
ra=abs(lada-ladb);
h=h+1;
end
%disp('15');
ladc=ladc';
mini=ladc;
end

function maxi = MetodHashinaShtrikmanaMax(laef,lan,por,ep,ladb,lada,hko)
ra=1e6; h=0;
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
ra=abs(lada-ladb);
h=h+1;
end
%disp('16');
ladc=ladc';
maxi=ladc;
end

function res = ReshNeyaUravKoTeplSwift(lae, lan, ep,ladb,lada,hko)
ra=1e5; h=0; la=0;
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
ra=abs(lada-ladb);
h=h+1;
end
%disp('17');
ladc=ladc';
res=ladc;
end

function fli = urovOtsechen(pro,ref,urot)
%po=30:10:90; po(length(po)+1)=95; po=1e-2*po; %uom=[220,150,110,80,60,45,30,20]; %urot=opredKTPTKTochSha(uom,por,po,length(uom)); %определить точнее
if (pro<=0)
    f=-1; pro=0; fl=0;
else
    if (pro>(urot*ref))
    f=-1;
else
    f=1;
    end
end
if (f>0)
fl=pro;
else
    fl=0;
end
fli=fl;
end

function t = SysPlastShakhPor(laef,lam2,po,ep,lamb,lama,nit)
ra=1e0; k=0;
while ((ra>ep) && (k<nit))
lamc=(lama+lamb)/2e0;
if (po<=5e-1)
fa=lam2*(4*po/(1+lam2/lama)+lama*(1-2*po)/lam2)-laef;
fb=lam2*(4*po/(1+lam2/lamb)+lamb*(1-2*po)/lam2)-laef;
fc=lam2*(4*po/(1+lam2/lamc)+lamc*(1-2*po)/lam2)-laef;
elseif (po<1e0)
fa=lam2*(4*(1-po)/(1+lam2/lama)+2*(po-1))-laef;
fb=lam2*(4*(1-po)/(1+lam2/lamb)+2*(po-1))-laef;
fc=lam2*(4*(1-po)/(1+lam2/lamc)+2*(po-1))-laef;
end
if ((fc*fb>0) && (fa*fc<0)) 
    lamb=lamc; 
elseif ((fc*fa>0) && (fb*fc<0)) 
    lama=lamc; 
end
k=k+1; ra=abs(lama-lamb); 
end
%disp('18');
labc=lamc';
t=lamc;
end

function t = MetodStarostina(laef,lam2,po,ep,lamb,lama,nit)
ra=1e0; k=0;
while ((ra>ep) && (k<nit))
lamc=(lama+lamb)/2e0;
fa=(lama^2)*(po^(2/3))+(lam2-lama)*lama;
fa=fa/(lama+((po^(2/3))-po)*(lam2-lama))-laef;
fb=(lamb^2)*(po^(2/3))+(lam2-lamb)*lamb;
fb=fb/(lamb+((po^(2/3))-po)*(lam2-lamb))-laef;
fc=(lamc^2)*(po^(2/3))+(lam2-lamc)*lamc;
fc=fc/(lamc+((po^(2/3))-po)*(lam2-lamc))-laef;
if ((fc*fb>0) && (fa*fc<0)) 
    lamb=lamc; 
elseif ((fc*fa>0) && (fb*fc<0)) 
    lama=lamc; 
end
k=k+1; ra=abs(lama-lamb);
end
%disp('19');
lamc=lamc';
t=lamc;
end

function t = MetodRusselNeprTverTelo(laef,lam2,po,ep,lamb,lama,nit)
ra=1e0; k=0;
while ((ra>ep) && (k<nit))
lamc=(lama+lamb)/2e0;
fa=lama*po+(lama/lam2)*(1-(po^(2/3)));
fa=fa/(po-(po^(2/3))+(lam2/lama)*(1-(po^(2/3))+po))-laef;
fb=lamb*po+(lamb/lam2)*(1-(po^(2/3)));
fb=fb/(po-(po^(2/3))+(lam2/lamb)*(1-(po^(2/3))+po))-laef;
fc=lamc*po+(lamc/lam2)*(1-(po^(2/3)));
fc=fc/(po-(po^(2/3))+(lam2/lamc)*(1-(po^2/3)+po))-laef;
if ((fc*fb>0) && (fa*fc<0)) 
    lamb=lamc; 
elseif ((fc*fa>0) && (fb*fc<0)) 
    lama=lamc; 
end
k=k+1; 
ra=abs(lama-lamb);
end
%disp('20');
lamc=lamc';
t=lamc;
end

function t = MetodRusselNeprVozd(laef,lam2,po,ep,lamb,lama,nit)
ra=1e0; k=0;
while ((ra>ep) && (k<nit))
lamc=(lama+lamb)/2e0;
fa=lama*((1-po)^2/3)+1-((1-po)^(2/3));
fa=fa/((lama/lam2)*(((1-po)^(2/3))-1+po)+(2-((1-po)^(2/3))-po))-laef;
fb=lamb*((1-po)^(2/3))+1-((1-po)^(2/3));
fb=fb/((lamb/lam2)*(((1-po)^(2/3))-1+po)+(2-((1-po)^(2/3))-po))-laef;
fc=lamc*((1-po)^(2/3))+1-((1-po)^(2/3));
fc=fc/((lamc/lam2)*(((1-po)^(2/3))-1+po)+(2-((1-po)^(2/3))-po))-laef;
if ((fc*fb>0) && (fa*fc<0)) 
    lamb=lamc; 
elseif ((fc*fa>0) && (fb*fc<0)) 
    lama=lamc; 
end
k=k+1; 
ra=abs(lama-lamb); 
end
%disp('21');
lamc=lamc';
t=lamc; 
end

function t = opredUrovPodderM03(por)
n=3; 
pam(3)=-8.93227577791347; 
pam(2)=8.89444783404518; 
pam(1)=1.06833435021354;
p=0; s=0; 
for k=1:n 
    s=s+pam(k)*pow(por,p); 
    p=p+1;
end
t=1/(1-s*por);
end

function opr = opredTeploprovVozd(T)
te0=273.15;
te = arrTempAir();
te=te+te0;
lamb = koefTeploprovAir();
n = length(te);
lam=opredKTPTKTochSha(lamb, te, T, n);
opr = lam;
end

function ktpo = opredKTPTKTochSha(ktptks, te, temp, n)
f = 1; p = 0; ep=1e-4; ktp = 0;
if ((temp>te(1)) && (temp<te(n)))
for k = 1:n
if ((te(k) >= temp) && (f>0))
        p = k; f = 0; break;
end
end
elseif (temp>=te(n))
    p=n; f=0;
elseif (temp<=te(1))
    p=2; f=0;
end
if ((f==0) && (p>1))
    x2=te(p);
    x1=te(p-1);
	dt = x2 - x1;
if (abs(dt) > ep)
    y2=ktptks(p);
    y1=ktptks(p - 1);
    b=y1;
			ko = (y2 - y1) / dt;
            if (p==n)
                b=y2;
            end
			ktp = b + ko*(temp - x1);
else 
    ktp=0;
end
end
ktpo=ktp;
end

function u = urovPod(po)
porex=95e-2; pormakvi=53e-2; 
porist=28e-2; porg=5e-1; 
if (po>6e-1) 
     m2=porex;
elseif (po>36e-2)
    m2=pormakvi;
else
    m2=po;
end
 if (po<porist)
     uo=opredUrovPodderM03(po); 
elseif (po<porg) 
    uo=1e0/(1e0-po); 
else
     uo=1/(1-m2);
 end
u=uo;
end

function [ vm ] = proverkakvi(ta, laefm, lavo, up, uo, n, ep)
for k=1:n
	ta(k)=urovOtsechen(ta(k), laefm(k), uo); 
	ta(k)=urovPodder(ta(k), lavo(k), up);
end
f=1; 
for k=1:n
    if ((ta(k)<ep) && (f>0)) 
        f=0; break;
    end
end
if (f==0) 
    for k=1:n 
        ta(k)=0; 
    end
end
vm=ta;
end