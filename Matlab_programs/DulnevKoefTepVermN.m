%��������� ���������� ��� �������� ������� ����������� ������� ��������
function [ dul ] = DulnevKoefTepVermN(vyb, tem, laef,  srp, gra, rp, legr, prgr, sr, por, lavo)
n=length(tem); gra=gra' ; vyb=vyb'; rp=rp'; te0=273.15;
lamax=1e6; lamin=1e-9; hko=10000; tocras=1e-8;
if (vyb<gra)
duln=BezUchetaRaspedPor(por, laef, lavo, vyb, srp, lamax, lamin, hko, tocras, n); %1 - ����, �����., 2 - ����, ����., 3 - �� ����������� 1, 4 - �� ����������� 2, 5 - ��������, �����., 6 - ��������, ����.
else
duln=UchetRapsredPorPoRazm(por, laef, lavo, vyb, rp, legr, prgr, sr, lamax, lamin, hko, tocras, n, tem); %7 - ����, �����., 8 - ����, �������., 9 - ��������, �����., 10 - ��������, ����., 11 - �� ����������� 1, 12 - �� ����������� 2
end
dul=duln;
end

function fl = urovOtsechen(pro,ref,urot)
if (pro<0)
    f=0;
elseif (pro>(urot*ref))
    f=0;
else f=1;
end
if (f>0)
fl=pro;
else
    fl=0;
end
end

%��� ����� ������������� ���
function [ BURP ] = BezUchetaRaspedPor(por, lameff, lamvoz, vy, l, lamax, lamin, hko, tocras, n)
up=1; por=por'; uo=urovPod(por);
lb = l / (por^(1 / 3));
lamadi = 0; lamizo = 0; lamodom=0; lamodos=0; lamcyladi=0; lamcylizo = 0; lamcylper=0;
for k=1:n
lamadi(k) = DulKoefTep1Adi(lamvoz(k), por, lameff(k), lamax, lamin, hko, tocras); %���������� ���������, ������� - ������������� ��������������
lamizo(k) = DulKoefTep1Izoterm(lamvoz(k), por, lameff(k), lamax, lamin, hko, tocras); %�������������� ���������, ������� - ������������� ��������������
lamodom(k) = DulKoefTep1OdolevMatr(lamvoz(k), por, lameff(k), lamax, lamin, hko, tocras); %�� ����������� (�������)
lamodos(k) = DulKoefTep1OdolevStatSm(lamvoz(k), por, lameff(k), lamax, lamin, hko, tocras); %�� ����������� (����. �����)
lamcyladi(k) = DulKoefTep1AdiCyl(lamvoz(k), por, lameff(k), lb, lamax, lamin, hko, tocras); %���������� ���������, ������� - �������
lamcylizo(k) = DulKoefTep1IzotermCyl(lamvoz(k), por, lameff(k), lb, lamax, lamin, hko, tocras); %�������������� ���������, ������� - �������
lamcylper(k) = DulKoefTep1CylPer(lamvoz(k), por, lameff(k),lb, lamax, lamin, hko, tocras); %��� �������� ��������������� ��������� ������
end
for k=1:n
switch (vy)
    case 0
        lamvy(k)=lamadi(k); %���, �����. - 23 %1
    case 1 
        lamvy(k)=lamizo(k); %���, �������. - 24 %2
    case 2 
        lamvy(k)=lamodom(k); %�� ����������� (�������) - 25 %3
    case 3
        lamvy(k)=lamodos(k); %�� ����������� (����. �����) - 26 %4
    case 4 
        lamvy(k)=lamcyladi(k); %�������., �����. - 27 %5
    case 5 
        lamvy(k)=lamcylizo(k); %�������., �������. - 28 %6
    case 6
        lamvy(k)=lamcylper(k); %��� �������� ��������������� ��������� ������ - 29 %7
    otherwise
        lamvy(k)=0;
end
%laef=lameff(k);%lavo=lamvoz(k); %lamvy(k)=urovOtsechen(lamvy(k),laef,uo); lamvy(k)=urovPodder(lamvy(k),lavo,up);
end
lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, n, tocras);
switch (vy)
    case 0
        %disp('23'); lamvy=lamvy'; %���, �����.
    case 1
        %disp('24'); lamvy=lamvy'; %���., �������.
    case 2
        %disp('25'); lamvy=lamvy'; %�� ����������� (�������)
    case 3
        %disp('26'); lamvy=lamvy'; %�� ����������� (����. �����)
    case 4
        %disp('27'); lamvy=lamvy'; %���., �����. 
    case 5
        %disp('28'); lamvy=lamvy'; %���., �������.
    case 6
        %disp('29'); lamvy=lamvy'; %��� �������� ��������������� ��������� ������
end
BURP=lamvy';
end

%���� ������������� ���
function [ URP ] = UchetRapsredPorPoRazm(po, lameff, lamvoz, vy, rp, legr, prgr, srk, lamax, lamin, hko, tocras, nT, tem)
le = length(legr); po=po'; %k=1; legr(k)=0; sr=0; 
for k=1:length(srk)
sr(k)=((legr(k)+prgr(k))/2e0+srk(k))/2e0;
end
e=1e-6;
j = 1; ep=tocras; hf=1e0;
if (po<hf)
    w=hf;
else w=1e2;
end
pw=zeros(1,le);
vo=zeros(1,le);
for k=1:le    
    pw(k) = rp(k) * po / w; %�������� ���� ���� ��������� ������� � ������ (�� ����) ������
    if (rp(k) > e)
        vo(k) = (sr(k)^3) / pw(k); %��� ���� - ����, ����� ��� ��������� ���������    
    else
        vo(k)=0;
    end
end
mx=sum(vo); up=1; uo=urovPod(po); 
if (vy==13)
pw=pw;
vo=vo;
sr=sr;
rp=rp';
lamvoz=lamvoz;
tem=tem;
uo=uo;
up=up;
po=po;
end
lcub = mx^(1/3); lcub=lcub/length(vo); %������ ������������� �������� l - ������� �������
lb = lcub / (po^(1/3)); %������ ������� ����
nvo=length(vo);
for k=1:nT
    lavo=lamvoz(k); lame=lameff(k); lamadit=lame; lamizot=lame;
    lamodo1(k)=lame; lamodo2(k)=lame; lamcyladi(k)=lame; lamcylizo(k)=lame; lamcylper(k)=lame;
    for j=1:nvo
        if (pw(j)>0)
lamadit = DulKoefTep1Adi(lavo, pw(j), lamadit, lamax, lamin, hko, tocras); %1
lamizot = DulKoefTep1Izoterm(lavo, pw(j), lamizot, lamax, lamin, hko, tocras); %2
lamodo1(k) = DulKoefTep1OdolevMatr(lavo, pw(j), lamodo1(k), lamax, lamin, hko, tocras); %3
lamodo2(k) = DulKoefTep1OdolevStatSm(lavo, pw(j), lamodo2(k), lamax, lamin, hko, tocras); %4
lamcyladi(k) = DulKoefTep1AdiCyl(lavo, pw(j), lamcyladi(k), lb, lamax, lamin, hko, tocras); %5
lamcylizo(k) = DulKoefTep1IzotermCyl(lavo, pw(j), lamcylizo(k), lb, lamax, lamin, hko, tocras); %6
lamcylper(k) = DulKoefTep1CylPer(lavo,pw(j),lamcylper(k),lb, lamax, lamin, hko, tocras); %7
        end
    end
        lamcubadia(k)=lamadit; lamcubizo(k)=lamizot;
        lamcub(k) = (lamcubadia(k)+lamcubizo(k))/2; %lamodo1(k) = lamodot1; %lamodo2(k) = lamodot2; %lamcyladi(k)=lamcyladit; %lamcylizo(k)=lamcylizot;
        lamcyl(k) = (lamcylizo(k)+lamcyladi(k))/2; %lamcylper(k)=lamcylpert;
end
for k=1:nT
switch (vy)
    case 7
        lamvy(k)=lamcubadia(k); %���, �����. - 30
    case 8
        lamvy(k)=lamcubizo(k); %���., �������. - 31
    case 9
        lamvy(k)=lamodo1(k); %�� ����������� (�������) - 32
    case 10
        lamvy(k)=lamodo2(k); %�� ����������� (����. �����) - 33
    case 11
        lamvy(k)=lamcyladi(k); %���., �����. - 34
    case 12
        lamvy(k)=lamcylizo(k); %���., �������. - 35
    case 13
        lamvy(k)=lamcylper(k); %��� �������� ��������������� ��������� ������ - 36
    otherwise 
        lamvy(k)=0;
end
end
lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, nT, tocras);
switch (vy)
    case 7
        %disp('30'); %lamvoz=lamvoz'
        %lamcubadia=lamcubadia'
        %lamvy=lamvy' %���, �����.
    case 8
        %disp('31'); 
        %lamcubizo=lamcubizo'
        %lamvy=lamvy' %���., �������.
    case 9
        %disp('32'); 
        %lamodo1=lamodo1' %�� ����������� (�������)
        %lamvy=lamvy'
    case 10
        %disp('33'); lamvy=lamvy' %�� ����������� (����. �����)
    case 11
        %disp('34'); lamvy=lamvy' %���., �����. 
    case 12
        %disp('35'); 
        %lamcylizo=lamcylizo' 
        %lamvy=lamvy' %���., �������.
    case 13
        %disp('36'); 
        %lamcylper=lamcylper'
        %lamvy=lamvy' %��� �������� ��������������� ��������� ������
end
URP = lamvy';
end

function u = urovPod(po)
porex=95e-2; pormakvi=53e-2; pomi=36e-2;
porist=28e-2; porg=5e-1; povi=6e-1;
if (po>povi) 
     m2=porex;
elseif (po>pomi)
    m2=pormakvi;
else
    m2=po;
end
 if (po<porist)
     uo=opredUrovPodderM03(po); 
elseif (po<porg) 
    uo=1e0/(1e0-po); 
else
     uo=1e0/(1e0-m2);
 end
u=uo;
end

function du1 = DulKoefTep1Adi(lam2, m2, lame, lamb, lama, nit, ep)
k=0; ra=1e4;
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
%k=k'
%ra=ra' 
%fc=fc' 
%lamc=lamc' 
du1 = lamc;
end

function du2 = DulKoefTep1Izoterm(lam2, m2, lame, lamb, lama, nit, ep)
k=0; ra=1e4;
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
%disp('23'); lamc=lamc';
du2 = lamc;
end
%��������� ������������ ������� - ���� ���� �������� ������� ������� ��� ����� �������� ������������ ���� ����, 
%������� ����� ��������� � ���� �����, ������ ������� �������� ������� ���������� �������, ����� �����������
function du31 = DulKoefTep1OdolevMatr(lam2, m2, lame, lamb, lama, nit, ep)
k=0; ra=1e4;
while ((ra > ep) && (k < nit))
    lamc = (lama + lamb) / 2;
    nua = 1 - lama / lam2; %lam2 - ��� �������
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
%lamc=lamc';
du31 = lamc;
end
%��� ���������� ���������, �������������� �����, ������� ������������ ����������
function du32 = DulKoefTep1OdolevStatSm(lam2, m2, lame, lamb, lama, nit, ep)
k=0; ra=1e4; %lam2 - ��� �������, ���� ��� �������� �������
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
%lamc=lamc';
du32 = lamc;
end

function du4 = DulKoefTep1AdiCyl(lam2, m2, lame, d, lamb, lama, nit, ep)
fra = 0.9; k=0; h = d; ra=1e4;
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
%lamc=lamc';
du4 = lamc;
end

function du5 = DulKoefTep1IzotermCyl(lam2, m2, lame, d, lamb, lama, nit, ep)
fra = 0.9; k=0; h = d; ra=1e4; 
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
%lamc=lamc';
du5 = lamc;
end

function du6 = DulKoefTep1CylPer(lavo,por,laef,lb, lamb, lama, nit, ep) %�������� ����� ��������������� ��� ��������
k=0; ra=1e4; lam1=lavo; dt=1e0; h=1; Qo=laef*h*dt; rb=sqrt(por/pi)*lb; %lam1 - ��� �������, lam2 - ��� �� 
%por=por'
%lb=lb'
while ((ra>ep) && (k<nit))
    lamc=(lama+lamb)/2;
	bb=dt*lam1*lama*h*rb; 
    bm=2*(lama-lam1)*rb; 
    am=lam1*lb;
    if (abs(am)>=abs(bm)) 
        ab=sqrt((am^2)-(bm^2)); 
        int1=2*(atan((am+bm)/ab)-atan(bm/ab))/ab;
    else
        ab=sqrt(bm^2-am^2); 
        int1=(log(abs(am+bm-ab)/abs(am+bm+ab))-log(abs(bm-ab)/abs(bm+ab)))/ab;
    end
    Q1=bb*(pi/2)/bm-am*bb/bm*int1; 
    Q1=2*Q1; Q2=lama*h*(lb-2*rb)/lb; 
    fa=Qo-Q1-Q2;
    bb=dt*lam1*lamb*h*rb;
    bm=2*(lamb-lam1)*rb;
if (abs(am)>=abs(bm)) 
        ab=sqrt(am^2-bm^2); 
        int1=2*(atan((am+bm)/ab)-atan(bm/ab))/ab;
    else
        ab=sqrt(bm^2-am^2); 
        int1=(log(abs(am+bm-ab)/abs(am+bm+ab))-log(abs(bm-ab)/abs(bm+ab)))/ab;
end
    Q1=bb*(pi/2)/bm-am*bb/bm*int1; 
    Q1=2*Q1; Q2=lamb*h*(lb-2*rb)/lb; 
    fb=Qo-Q1-Q2;
    bb=dt*lam1*lamc*h*rb; 
    bm=2*(lamc-lam1)*rb;
if (abs(am)>=abs(bm)) 
        ab=sqrt(am^2-bm^2); 
        int1=2*(atan((am+bm)/ab)-atan(bm/ab))/ab; 
    else
        ab=sqrt(bm^2-am^2); 
        int1=(log(abs(am+bm-ab)/abs(am+bm+ab))-log(abs(bm-ab)/abs(bm+ab)))/ab;
end
    Q1=bb*(pi/2)/bm-am*bb/bm*int1; 
    Q1=2*Q1; 
    Q2=lamc*h*(lb-2*rb)/lb; 
    fc=Qo-Q1-Q2;
    if ((fc*fb>0) && (fa*fc<0)) 
        lamb=lamc; 
    end
    if ((fc*fa>0) && (fb*fc<0)) 
        lama=lamc;
    end
	k=k+1; 
    ra=abs(lama-lamb);
end
%lamc=lamc'
du6=lamc;
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

function [ vm ] = proverkakvi(ta, laefm, lavo, up, uo, n, ep)
lvmi=min(lavo);
for k=1:n
	ta(k)=urovOtsechen(ta(k), laefm(k), uo); 
	ta(k)=urovPodder(ta(k), lavo(k), up);
end
f=1;
for k=1:n
    if ((ta(k)<lvmi) && (f>0)) 
        f=-1; break;
    end
end
if (f<0) 
    for k=1:n 
        ta(k)=0; 
    end
end
vm=ta;
end

function urpo = urovPodder(pro,ref,urpo)
f=1; 
if (pro<0) 
    f=0; 
elseif (pro<(urpo*ref)) 
    f=0;
else f=1;
end
if (f>0) 
    fl=pro;
else
    fl=0;
end
urpo=fl;
end