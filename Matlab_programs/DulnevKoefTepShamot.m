%программа определяет КТП твердого каркаса шамота методом Дульнева
function [ dul ] = DulnevKoefTepShamot(vyb, por, tem, te0, vybsha, srp,gra)
if (vybsha == 0) 
    koef=[-0.435e-9,0.685e-6,0.134e-3,0.725]; 
end
if (vybsha == 2)
    koef=[-0.377e-9,0.918e-6,-0.338e-3,0.77];
end
n=length(tem); ep=1e-8;
for k=1:n
ts=tem(k)-te0;
lavo(k)=opredTeploprovVozd(ts);
laef(k)=polyval(koef,ts);
end
if (vyb<gra)
du=BezUchetaRaspedPor(por, tem-te0, laef, lavo,vyb,ep,srp); %1 - Кубы, адиаб., 2 - Кубы, изот., 3 - По Одолевскому 1, 4 - По Одолевскому 2, 5 - Цилиндры, адиаб., 6 - Цилиндры, изот.
else
du=UchetRapsredPorPoRazm(por, tem-te0, laef, lavo,vyb,ep); %7 - Кубы, адиаб., 8 - Кубы, изотерм., 9 - Цилиндры, адиаб., 10 - Цилиндры, изот., 11 - По Одолевскому 1, 12 - По Одолевскому 2
end
dul=du;
end

function fli = urovPodder(pro,ref,urpo)
fl=0; f=1; 
if (pro<=0) 
    f=0; 
else
    if (abs(pro)<=abs(urpo*ref))
    f=0; 
    else f=1;
    end
end
if (f>0) 
    fl=pro;
end
fli=fl;
end

function fl = urovOtsechen(pro,ref,urot)
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
%без учета распределения пор
function [ BURP ] = BezUchetaRaspedPor(por, tem, lameff, lamvoz,vy,ep,l)
lamvy=0;
up=1; 
if (po>5e-1) 
     m2=95e-2; 
else m2=po;
 end
 if (po<28e-2) 
     uo=opredUrovPodderM03(por); 
 else
     uo=1/(1-m2);
 end
%l = srRaPor();
n=length(tem);
lb = l / (por^(1 / 3));
lamadi = 0; lamizo = 0; lamodom=0; lamodos=0; lamcyladi=0; lamcylizo = 0;
for k=1:n
laef=lameff(k);
lavo=lamvoz(k);
lamadi(k) = DulKoefTep1Adi(lavo, por, laef, ep, k, vy); %адиабатное разбиение, полость - прямоугольный параллелепипед
lamizo(k) = DulKoefTep1Izoterm(lavo, por, laef, ep); %изотермическое разбиение, полость - прямоугольный параллелепипед
lamodom(k) = DulKoefTep1OdolevMatr(lavo, por, laef, ep); %по Одолевскому (матрицы)
lamodos(k) = DulKoefTep1OdolevMatr(lavo, por, laef, ep); %по Одолевскому (стат. смесь)
lamcyladi(k) = DulKoefTep1AdiCyl(lavo, por, laef, lb, ep); %адиабатное разбиение, полость - цилиндр
lamcylizo(k) = DulKoefTep1IzotermCyl(lavo, por, laef, lb, ep); %изотермическое разбиение, полость - цилиндр, ось цилиндра идет параллельно тепловому потоку
lamcylizoper(k) = DulKoefTep1CylPer(lavo, por, laef, lb, ep); %изотермическое разбиение, полость - цилиндр, ось цилиндра направлена перпендикулярно тепловому потоку
end
for k=1:n
    if (vy==1) %19
        lamvy(k)=urovOtsechen(lamadi(k),lameff(k),uo); %куб, адиаб.
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==2) %20
        lamvy(k)=urovOtsechen(lamizo(k),lameff(k),uo); %куб, изотерм.
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==3) %21
        lamvy(k)=urovOtsechen(lamodom(k),lameff(k),uo); %по Одолевскому (матрицы)
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==4) %22
        lamvy(k)=urovOtsechen(lamodos(k),lameff(k),uo); %по Одолевскому (стат. смесь)
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==5) %23
        lamvy(k)=urovOtsechen(lamcyladi(k),lameff(k),uo); %цилиндр., адиаб.
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==6) %24
        lamvy(k)=urovOtsechen(lamcylizo(k),lameff(k),uo); %цилиндр., изотерм.
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==7) %25
        lamvy(k)=urovOtsechen(lamcylizoper(k),lameff(k),uo); %цилиндр., изотерм.
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
end
BURP=lamvy';
end
%учет распределения пор
function [ URP ] = UchetRapsredPorPoRazm(po, T, lameff, lamvoz, vy, ep)
rpr = 0;
if (po>5e-1) 
     m2=95e-2; 
else m2=po;
 end
 if (po<28e-2) 
     uo=opredUrovPodderM03(por); 
 else
     uo=1/(1-m2);
 end
lamvy=0;
raspr=[21.8,21.75,21.625,21.25,20.75,20.25,19.75,19,18.25,17.5,12,9.875,8.75,8.25,7.25,6.25,5,4.25,3.25,0.875,0.325,0];
for k=2:length(raspr)
    rpr(k)=(raspr(k-1)-raspr(k))*1e2/raspr(1);
end
prgr=[0.05,0.1,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50];
legr(1)=0;
for k=2:length(prgr)
    legr(k)=prgr(k-1);
end
sr = (legr + prgr)*1e-6 / 2; %средний размер пор в метрах
le = length(sr);
j = 1; vo = 0; pw = 0; mx = 0;
for k=1:le
    if (rpr(k) > 0)
    pw(j) = rpr(k) * po / 1e2; %объемная доля поры заданного размера в полном (во всем) объеме
    vo(j) = sr(k)^3 / pw(j); %все поры - кубы
    mx = mx + vo(j);
    j = j + 1;
    end
end
mx=mx';
lcub = mx^(1/3); %оценка максимального значения l - размера образца
lb = lcub / (po^(1/3)); %оценка размера поры
nT=length(T);
nvo=length(pw);
lamcubadia = 0; lamcubizo = 0; lamcyladi=0; lamcylizo = 0; 
lamcub = 0; lamodo1=0; lamodo2=0; lamcyl = 0;
for k=1:nT
    lavo=lamvoz(k); lame=lameff(k); lamadit = lame; lamizot = lame; lamodot1 = lame; 
    lamodot2 = lame; lamcyladit = lame; lamcylizot = lame; lamcylpert=lame;
    for j=1:nvo
lamadit = DulKoefTep1Adi(lavo, pw(j), lamadit, ep, k, vy);
lamizot = DulKoefTep1Izoterm(lavo, pw(j), lamizot, ep);
lamodot1 = DulKoefTep1OdolevMatr(lavo, pw(j), lamodot1, ep);
lamodot2 = DulKoefTep1OdolevStatSm(lavo, pw(j), lamodot2, ep);
lamcyladit = DulKoefTep1AdiCyl(lavo, pw(j), lamcyladit, lb, ep);
lamcylizot = DulKoefTep1IzotermCyl(lavo, pw(j), lamcylizot, lb, ep);
lamcylpert=DulKoefTep1CylPer(lavo, pw(j), lamcylpert, lb, ep); 
    end
    lamcubadia(k)=lamadit; lamcubizo(k)=lamizot;
        lamcub(k) = (lamadit+lamizot)/2;
        lamodo1(k) = lamodot1; lamodo2(k) = lamodot2;
    lamcyladi(k)=lamcyladit; lamcylizo(k)=lamcylizot;
        lamcyl(k) = (lamcylizot+lamcyladit)/2;
        lamcylper(k)=lamcylpert;
end
for k=1:nT
    lavo=lamvoz(k);
    if (vy==8) %26 проверить
        lamvy(k)=urovOtsechen(lamcubadia(k),lameff(k),uo); %куб, адиаб.
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==9) %27
        lamvy(k)=urovOtsechen(lamcubizo(k),lameff(k),uo); %куб., изотерм.
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==10) %28 проверить
        lamvy(k)=urovOtsechen(lamcyladi(k),lameff(k),uo); %цил., адиаб. 
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==11) %29 проверить
        lamvy(k)=urovOtsechen(lamcylizo(k),lameff(k),uo); %цил., изотерм.
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==12) %30
        lamvy(k)=urovOtsechen(lamodo1(k),lameff(k),uo); %по Одолевскому (матрицы)
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==13) %31
        lamvy(k)=urovOtsechen(lamodo2(k),lameff(k),uo); %по Одолевскому (стат. смесь)
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
    if (vy==14) %32
        lamvy(k)=urovOtsechen(lamcylper(k),lameff(k),uo); %цилиндры - ось перпендикулярно тепловому потоку
        lamvy(k)=urovPodder(lamvy(k),lavo,up);
    end
end
URP = lamvy';
end
function SRP = srRaPor()
rp = 0;
raspr=[21.8,21.75,21.625,21.25,20.75,20.25,19.75,19,18.25,17.5,12,9.875,8.75,8.25,7.25,6.25,5,4.25,3.25,0.875,0.325,0];
for k=2:length(raspr)
    drasp(k)=(raspr(k-1)-raspr(k))*1e2/raspr(1);
end
prgr=[0.05,0.1,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50];
legr(1)=0;
for k=2:length(prgr)
    legr(k)=prgr(k-1);
end
for k=1:length(prgr)
    srra(k)=(legr(k)+prgr(k))/2;
end
s=0;
for k=1:length(raspr)
    s=s+srra(k)*drasp(k)*1e-2;
end
SRP = s;
end
function du1 = DulKoefTep1Adi(lam2, m2, lame, ep, u, y)
lamb = 1e4;
lama = 1e-6;
k=0;
nit=1e2;
ra=1;
while ((ra > ep) && (k<nit))
    lamc = (lama + lamb) / 2;
    nua = lam2 / lama;
    nub = lam2 / lamb;
    nuc = lam2 / lamc;
    fa = (nua - (nua - 1) * (m2^(1 / 3)) * (1 - m2^(2 / 3))) / (nua - (m2^(1/3)) * (nua - 1)) - lame / lama;
    fb = (nub - (nub - 1) * (m2^(1 / 3)) * (1 - m2^(2 / 3))) / (nub - (m2^(1/3)) * (nub - 1)) - lame / lamb;
    fc = (nuc - (nuc - 1) * (m2^(1 / 3)) * (1 - m2^(2 / 3))) / (nuc - (m2^(1/3)) * (nuc - 1)) - lame / lamc;
    if ((fc*fb > 0)  && (fa*fc < 0))
            lamb=lamc; 
    end
    if ((fc*fa > 0)  && (fb*fc < 0)) 
        lama=lamc; 
    end
    k=k+1;
    ra=abs(lama - lamb);
end
if ((u==1) && (y==7))
lamc=lamc';
end
du1 = lamc;
end
function du2 = DulKoefTep1Izoterm(lam2, m2, lame, ep)
lamb = 1e4;
lama = 1e-6;
k=0;
ra=1;
nit=1e2;
while ((ra > ep) && (k<nit))
    lamc = (lama + lamb) / 2;
    nua = lam2 / lama;
    nub = lam2 / lamb;
    nuc = lam2 / lamc;
    fa = (1 + (nua - 1) * (m2^(2 / 3)) ) / (1 + (m2^(2 / 3)) * (nua - 1)) * (1 - m2^(1 / 3)) - lame / lama;
    fb = (1 + (nub - 1) * (m2^(2 / 3)) ) / (1 + (m2^(2 / 3)) * (nub - 1)) * (1 - m2^(1 / 3)) - lame / lamb;
    fc = (1 + (nuc - 1) * (m2^(2 / 3)) ) / (1 + (m2^(2 / 3)) * (nuc - 1)) * (1 - m2^(1 / 3)) - lame / lamc;
    if ((fc*fb > 0)  && (fa*fc < 0))
        lamb=lamc; 
    end
    if ((fc*fa > 0) && (fb*fc < 0)) 
        lama=lamc; 
    end
    k=k+1;
    ra=abs(lama - lamb);
end
du2 = lamc;
end
%матричная гетерогенная система - одна фаза образует связную матрицу при любой объемной концентрации этой фазы, 
%система имеет включения в виде кубов, центры которых образуют простую кубическую решетку, ребра параллельны
function du31 = DulKoefTep1OdolevMatr(lam2, m2, lame, ep)
lamb = 1e4; lama = 1e-9; k=0; nit=1e2; ra=1e4;
while ((ra > ep) && (k < nit))
    lamc = (lama + lamb) / 2;
    nua = 1 - lama / lam2; %lam2 - КТП воздуха
    nub = 1 - lamb / lam2;
    nuc = 1 - lamc / lam2;
    fa = 1 - (1 - m2) / (1 / nua - m2 / 3) - lame / lam2;
    fb = 1 - (1 - m2) / (1 / nub - m2 / 3) - lame / lam2;
    fc = 1 - (1 - m2) / (1 / nuc - m2 / 3)  - lame / lam2;
    if ((fc*fb > 0) && (fa*fc < 0)) 
        lamb=lamc; 
    end
    if ((fc*fa > 0) && (fb*fc < 0)) 
        lama=lamc; 
    end
    k=k+1;
    ra=abs(lama - lamb);
end
du31 = lamc;
end
%нет регулярной структуры, статистическая смесь, частицы распределены хаотически
function du32 = DulKoefTep1OdolevStatSm(lam2, m2, lame, ep)
lamb = 1e4; lama = 1e-9; k=0; nit=1e3; ra=1e4; %lam2 - КТП воздуха, ищем КТП твердого скелета
while ((ra > ep) && (k < nit))
    lamc = (lama + lamb) / 2;
    v1 = m2; v2 = 1 - m2;
    nua = ((3*v1 - 1) * lama + (3*v2 - 1) * lam2) / 4;
    nub = ((3*v1 - 1) * lamb+ (3*v2 - 1) * lam2) / 4;
    nuc = ((3*v1 - 1) * lamc + (3*v2 - 1) * lam2) / 4;
    fa = nua + sqrt (nua^2 + lama * lam2 / 2) - lame;
    fb = nub + sqrt (nub^2 + lamb * lam2 / 2) - lame;
    fc = nuc + sqrt (nuc^2 + lamc * lam2 / 2) - lame;
    if ((fc*fb > 0) && (fa*fc < 0)) 
        lamb=lamc; 
    end
    if ((fc*fa > 0) && (fb*fc < 0)) 
        lama=lamc; 
    end
    k=k+1;
    ra=abs(lama - lamb);
end
du32 = lamc;
end
function du4 = DulKoefTep1AdiCyl(lam2, m2, lame, d, ep)
fra = 0.9;
lamb = 1e4;
lama = 1e-6;
k=0;
h = d; 
ra=1;
nit=1e2;
while ((ra > ep) && (k<nit))
    lamc = (lama + lamb) / 2;
    h1 = fra * d;
    d1 = d * sqrt(4 * m2 / pi / fra);
    F1 = ( pi / 4 ) * d1^2; 
    R1a = (h - h1) / 2 / F1 / lama;
    R1b = (h - h1) / 2 / F1 / lamb;
    R1c = (h - h1) / 2 / F1 / lamc;
    R2 = (h - h1) / 2 / F1 / lam2;
    F12 = ( pi / 4 ) * (d^2 - d1^2); 
    R3a = h / F12 / lama;
    R3b = h / F12 / lamb;
    R3c = h / F12 / lamc;
    F = (pi / 4) * d^2; 
    R = h / lame / F;
    fa = R3a * (2 * R1a + R2) / (2 * R1a + R2 + R3a) - R;
    fb = R3b * (2 * R1b + R2) / (2 * R1b + R2 + R3b) - R;
    fc = R3c * (2 * R1c + R2) / (2 * R1c + R2 + R3c) - R;
    if ((fc * fb > 0) && (fa * fc < 0) )
            lamb=lamc; 
    end
    if ((fc * fa > 0) && (fb * fc < 0))
            lama=lamc; 
    end
    k=k+1;
ra=abs(lama - lamb);
end
du4 = lamc;
end
function du5 = DulKoefTep1IzotermCyl(lam2, m2, lame, d, ep)
fra = 0.9;
lamb = 1e5;
lama = 1e-6;
k=0; 
h = d;
ra=1;
nit=1e2;
while ((ra > ep) && (k<nit))
    lamc = (lama + lamb) / 2;
    h1 = fra * d;
    d1 = d * sqrt(4 * m2 / pi / fra);
    F = ( pi / 4 ) * d^2; 
    R1a = (h - h1) / 2 / F / lama;
    R1b = (h - h1) / 2 / F / lamb;
    R1c = (h - h1) / 2 / F / lamc;
    F1 = ( pi / 4 ) * d1^2; 
    R2 = h1 / F1 / lam2;
    F12 = ( pi / 4 ) * (d^2 - d1^2); 
    R3a = h1 / F12 / lama;
    R3b = h1 / F12 / lamb;
    R3c = h1 / F12 / lamc;
    F = ( pi / 4 ) * d^2; 
    R = h / lame / F;
    fa = 2 * R1a + R3a * R2 / (R2 + R3a) - R;
    fb = 2 * R1b + R3b * R2 / (R2 + R3b) - R;
    fc = 2 * R1c + R3c * R2 / (R2 + R3c) - R;
    if ((fc * fb > 0) && (fa * fc < 0))
            lamb=lamc; 
    end
    if ((fc * fa > 0) && (fb * fc < 0) )
            lama=lamc; 
    end
    k=k+1;
    ra=abs(lama - lamb);
end
du5 = lamc;
end
function du6 = DulKoefTep1CylPer(lavo, por, laef, lb, ep) %lam1 - КТП воздуха, lam2 - КТП ТТ
lamb=1e5; lam1=lavo;
lama=1e-6;
k=0; ra=1;
nit=1e2; dt=1; h=1;
Qo=laef*h*dt;
rb=sqrt(por/pi)*lb;
while ((ra>ep) && (k<nit))
    lamc=(lama + lamb)/2;
    bb=dt*lam1*lama*h*rb;
    bm=(2*(lama-lam1)*rb);
    am=(lam1*lb);
    if (abs(am)>abs(bm))
        ab=sqrt(am^2-bm^2);
        int1=2*(atan((am+bm)/ab)-atan(bm/ab))/ab;
    else
        ab=sqrt(bm^2-am^2);
        int1=(log((am+bm-ab)/(am+bm+ab))-log((bm-ab)/(bm+ab)))/ab;
    end
    Q1=bb*(pi/2)/bm-am*bb/bm*int1; Q1=2*Q1;
    Q2=lama*h*(lb-2*rb)/lb;
    fa=Qo-Q1-Q2;
    bb=dt*lam1*lamb*h*rb;
    bm=(2*(lamb-lam1)*rb);
    if (abs(am)>abs(bm))
        ab=sqrt(am^2-bm^2);
        int1=2*(atan((am+bm)/ab)-atan(bm/ab))/ab;
    else
        ab=sqrt(bm^2-am^2);
        int1=(log((am+bm-ab)/(am+bm+ab))-log((bm-ab)/(bm+ab)))/ab;
    end
    Q1=bb*(pi/2)/bm-am*bb/bm*int1; Q1=2*Q1;
    Q2=lamb*h*(lb-2*rb)/lb;
    fb=Qo-Q1-Q2;
    bb=dt*lam1*lamc*h*rb;
    bm=(2*(lamc-lam1)*rb);
    if (abs(am)>abs(bm))
        ab=sqrt(am^2-bm^2);
        int1=2*(atan((am+bm)/ab)-atan(bm/ab))/ab;
    else
        ab=sqrt(bm^2-am^2);
        int1=(log((am+bm-ab)/(am+bm+ab))-log((bm-ab)/(bm+ab)))/ab;
    end
    Q1=bb*(pi/2)/bm-am*bb/bm*int1; Q1=2*Q1;
    Q2=lamc*h*(lb-2*rb)/lb;
    fc=Qo-Q1-Q2;
    if ((fc*fb>0) && (fa*fc<0))
            lamb=lamc; 
    end
    if ((fc*fa>0) && (fb*fc<0) )
            lama=lamc; 
    end
    k=k+1;
    ra=abs(lama - lamb);
end
du6 = lamc;
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