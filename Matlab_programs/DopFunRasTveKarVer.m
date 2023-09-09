function [ t ] = DopFunRasTveKarVer(tem, no, m22, r, laefm, lentem, lavo, stchm)
ta=0; 
switch (no)
    case (0)
    ta=MetodDulnevaSigalovoyPolyDispVer(tem,m22,r,laefm,stchm);
    case (1)
    ta=MetodDulnevaSigalovoyDopVer(tem,m22,r,laefm,stchm);
    case (2)
    ta=MetodDulnevaSigalovoyBezIspVer(tem,m22,r,laefm,stchm);
    case (3)
    ta=VasFrayObsh(laefm, lavo, m22);
end
up=1;
uo=urovPod(m22);
for k=1:length(ta)
	ta(k)=urovOtsechen(ta(k), laefm(k), uo); 
	ta(k)=urovPodder(ta(k), lavo(k), up);
end
t=ta;
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

function [ lam ] = MetodDulnevaSigalovoyPolyDispVer(tem,m22,r,laefm,stchm) %m2 - пористость
m20 = 95e-2; porex=m20; por=m22; %эксперимент
m2 = m22 - m20; mg=26e-2; %межзерновая пористость (внешняя)
if (m2<mg)
    m2=0.365; %межзерновая пористость по Нижегородову
end
te0=273.15; %dan=arrTem_VVF2(); %dann=arrKTP_VVF2(); %фракция 2-0,7 мм, измерения на немецкой установке %kektp=polyfit(dan,dann,2);
kektp=polyfit(tem,laefm,2);
A=(0.74/(1-m2))^(1/3); %alpha - размер воздушной оболочки (ореола)
gu=9.8; hv=30e-3; rov=125;
sigb=2.156e5; %предел прочности при сжатии материала частиц, Па
ka=7/5; H0=101325; H=1e5; kB=1.380658e-23;
d=(0.6*0.21+0.65*0.79)*1e-10; sigse=pi*(d^2)/4;
E=8.428e6; %модуль Юнга вермикулита
akk=9e-1; %коэффициент аккомодации
ko1=0.2595/(1-0.2595); ko2=0.2595-ko1; fp=0.2595/m2-(ko1*m2+ko2);
%fileID = fopen('Stepen_chernoty_ver.txt','r'); formatSpec='%f'; stchm=fscanf(fileID,formatSpec); fclose(fileID);
for k=1:length(tem)
ladb=1e5; lada=1e-5; ra=1e2; ep=1e-7; h=0; kit=1e3;
while ((ra>ep) && (h<kit))
ladc=(lada+ladb)/2;
T=tem(k);
laef=polyval(kektp,T);
stch=opredKTPTKTochSha(stchm,tem,T,length(tem));
izco=2*(5.67e-2)*(stch^2)*3/4;
    lam0=opredKTPVozd(T); %температура в К
    n=H/kB/T;
    dliSvoPro=1/sqrt(2)/n/sigse; %длина свободного пробега молекулы воздуха
    cPr=opredPrVozd(T-te0); %число Прандтля
    N=sqrt(m2^2-10*m2+9); N=(m2+3+N)/2/m2;
    debo=rov*gu*hv; %удельная нагрузка
    %dem1=2*r*(A-pi/4);
    dem1=r*(A-pi/4);
    dem2=7*r*(debo^(7/18))*sqrt(sigb)/sqrt(N)/(E^(8/9));
    dem4=2*r*(sqrt(2)*A-1);
    dem5=2*sqrt(2*r)*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4*ka/(ka+1)/cPr*((2-akk)/akk)*H0*dliSvoPro;
    lamg1=lam0/(1+B/H/dem1);
    lamg2=lam0/(1+B/H/dem2);
    lamg4=lam0/(1+B/H/dem4);
    lamg5=lam0/(1+B/H/dem5);
    sigm1=lamg1*pi*r*(A*log(A/(A-1))-1);
    sigm2=lamg2*sqrt(N)*(E^(4/9))*(debo^(1/18))*r/3.2/sqrt(sigb);
    sigm2=sigm2*fp;
    sigm3=2e-2*sqrt(N)*(E^(8/9))*(debo^(11/18))*r/(sigb^(3/2));
    sigm3a=lada*sigm3;
    sigm3b=ladb*sigm3;
    sigm3c=ladc*sigm3;
    sigm4=lamg4*pi*r*(sqrt(2)*A-1)/2; %большая пора - цилиндр
    if (m2>=0.74)
        sigm5=lamg5*2*r*(2*(A^2)-pi)/sqrt(2)/A;
        rr=1+(A^2)/pi;
    else
        sigm5=0;
        rr=(1+((sqrt(2)*A-1)^2)/2);
    end
    gi=izco*((T/1e2)^3); 
    sigmr=pi*(r^2)*gi*rr;
    fa=(sqrt(2)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3a)/A/r-laef;
    fb=(sqrt(2)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3b)/A/r-laef;
    fc=(sqrt(2)*(sigm1+sigm2+sigm4/2+sigm5+sigmr)+sigm3c)/A/r-laef;
    la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
up=1;
m2=porex; 
if (por>5e-1)
    m2=porex;
else m2=po;
end
if (por<28e-2) 
uo=opredUrovPodderM03(por); 
else uo=1e0/(1e0-m2);
end
T=tem(k);
laef=polyval(kektp,T);
lam0=opredKTPVozd(T); %температура в К
ladc=urovOtsechen(ladc,laef,1/m22);
ladc=urovPodder(ladc,lam0,1+m22);
ladi(k)=ladc;
end
lam=ladi;
end

function [ lam ] = MetodDulnevaSigalovoyDopVer(tem,m22,r,laefm,stchm) %m2 - пористость %Из статьи "Теплопроводность моно- и полидисперсных материалов"
m20 = 0.95; porex=m20; por=m22; %эксперимент
m2 = m20 - m22;  te0=273.15; %dan=arrTem_VVF2()+te0; %dann=arrKTP_VVF2(); %фракция 2-0,7 мм, измерения на немецкой установке %kektp=polyfit(dan,dann,2);
mg=26e-2; 
if (m2<mg)
    m2=0.365; %межзерновая пористость по Нижегородову
end
kektp=polyfit(tem,laefm,2);
A=(0.74/(1-m2))^(1/3); %alpha - размер воздушной оболочки (ореола)
gu=9.8; hv=30e-3; rov=125;
sigb=2.156e5; %предел прочности при сжатии материала частиц, Па
ka=7/5; H0=101325; H=1e5; kB=1.380658e-23;
d=(0.6*0.21+0.65*0.79)*1e-10; sigse=pi*(d^2)/4;
E=8.428e6; %модуль Юнга вермикулита
akk=0.9; %коэффициент аккомодации
ko1=0.2595/(1-0.2595); ko2=0.2595-ko1; fp=0.2595/m2-(ko1*m2+ko2);
%fileID = fopen('Stepen_chernoty_ver.txt','r'); formatSpec='%f'; stchm=fscanf(fileID,formatSpec); fclose(fileID);
for k=1:length(tem)
ladb=1e5; lada=1e-5; ra=1e2; ep=1e-7; h=0; kit=1e3;
while ((ra>ep) && (h<kit))
ladc=(lada+ladb)/2;
T=tem(k);
laef=polyval(kektp,T);
stch=opredKTPTKTochSha(stchm,tem,T,length(tem));
izco=2*(5.67e-2)*(stch^2)*3/4;
    lam0=opredKTPVozd(T); %температура в К
    n=H/kB/T;
    dliSvoPro=1/sqrt(2)/n/sigse; %длина свободного пробега молекулы воздуха
    cPr=opredPrVozd(T-te0); %число Прандтля
    %N=sqrt(m2^2-10*m2+9); N=(m2+3+N)/2/m2;
    dem1=2*r*(A-pi/4);
    dem3=2*r*(1.41*A-1);
    dem5=2*sqrt(2)*r*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4*ka/(ka+1)*((2-akk)/akk)*H0*dliSvoPro/cPr;
    X=4.45*(A*log(A/(A-1))-1)/(1+B/H/dem1);
    Zs=2.23*(1.41*A-1)/(1+B/H/dem3);
    Zss=2.23/(1+B/H/dem3)/(1.41*A-1);
    V=(A^2-1.57)/(1+B/H/dem5)/A;
    G=izco*((T/1e2)^3); 
    lamrs=4.45*G*r*(1+0.5*(1.41*A-1)^2)/A;
    lamrss=4.45*G*r*(1+A^2/pi);
    if (m2>=0.74)
    fa=lam0*(X+Zss+V)/A+lamrss+lada-laef;    
    fb=lam0*(X+Zss+V)/A+lamrss+ladb-laef;
    fc=lam0*(X+Zss+V)/A+lamrss+ladc-laef;    
    elseif (m2>=0.26)
    fa=lam0*(X+Zs)/A+lamrs+lada-laef;    
    fb=lam0*(X+Zs)/A+lamrs+ladb-laef;
    fc=lam0*(X+Zs)/A+lamrs+ladc-laef;
    end
    la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
up=1;
m2=porex; 
if (por<5e-1)
    m2=po;
end
if (por<28e-2) 
uo=opredUrovPodderM03(por); 
else uo=1e0/(1e0-m2);
end
T=tem(k);
laef=polyval(kektp,T);
lam0=opredKTPVozd(T); %температура в К
ladc=urovOtsechen(ladc,laef,uo);
ladc=urovPodder(ladc,lam0,up);
ladi(k)=ladc;
end
lam=ladi;
end

function [ lam ] = MetodDulnevaSigalovoyBezIspVer(tem,m22,r,laefm,stchm) %m2 - пористость %Из статьи "Теплопроводность зернистых систем"
m20 = 0.95; porex=m20; por=m22; %эксперимент
m2 = m20 - m22; te0=273.15; %dan=arrTem_VVF2()+te0; %dann=arrKTP_VVF2(); %фракция 2-0,7 мм, измерения на немецкой установке %kektp=polyfit(dan,dann,2);
mg=26e-2; %межзерновая пористость (внешняя)
if (m2<mg)
    m2=0.365; %межзерновая пористость по Нижегородову
end
kektp=polyfit(tem,laefm,2);
A=((1-0.2595)/(1-m2))^(1/3); %alpha - размер воздушной оболочки (ореола)
gu=9.8; hv=30e-3; rov=125; 
sigb=2.156e5; %предел прочности при сжатии материала частиц, Па
ka=7/5; H0=101325; H=1e5; kB=1.380658e-23;
d=(0.6*0.21+0.65*0.79)*1e-10; sigse=pi*(d^2)/4; %считаем, что воздух состоит из азота и кислорода
E=8.428e6; %модуль Юнга вермикулита
akk=0.9; %коэффициент аккомодации
ko1=0.2595/(1-0.2595); ko2=0.2595-ko1; fp=0.2595/m2-(ko1*m2+ko2);
%fileID = fopen('Stepen_chernoty_ver.txt','r'); formatSpec='%f'; stchm=fscanf(fileID,formatSpec); fclose(fileID);
for k=1:length(tem)
ladb=1e5; lada=1e-5; ra=1e2; ep=1e-7; h=0; kit=1e3;
while ((ra>ep) && (h<kit))
ladc=(lada+ladb)/2;
T=tem(k);
laef=polyval(kektp,T-te0);
stch=opredKTPTKTochSha(stchm,tem,T,length(tem));
izco=2*(5.67e-2)*(stch^2)*3/4;
    lam0=opredKTPVozd(T); %температура в К
    n=H/kB/T;
    dliSvoPro=1/sqrt(2)/n/sigse; %длина свободного пробега молекулы воздуха
    cPr=opredPrVozd(T-te0); %число Прандтля
    N=sqrt(m2^2-10*m2+9); N=(m2+3+N)/2/m2;
    debo=rov*gu*hv; %удельная нагрузка
    %dem1=2*r*(A-pi/4);
    dem1=r*(A-pi/4);
    dem2=7*r*(debo^(7/18))*sqrt(sigb)/sqrt(N)/(E^(8/9));
    dem4=2*r*(1.41*A-1);
    dem5=2.82*r*A;
    akk=opredKoefAkkomodLandau(T); 
    B=4*ka/(ka+1)/cPr*((2-akk)/akk)*H0*dliSvoPro;
    ch1=4.45*(A*log(A/(A-1))-1)/(1+B/H/dem1);
    ch2=fp*0.44*sqrt(N)*(E^(4/9))*(debo^(1/18))/sqrt(sigb)/(1+B/H/dem2);
    ch3=2.23*(1.41*A-1)/(1+B/H/dem4);
    ch123=lam0*(ch1+ch2+ch3)/A;
    ch4a=(2e-2)*fp*lada*sqrt(N)*(E^(8/9))*(debo^(11/18))/(sigb^(3/2))/A;
    ch4b=(2e-2)*fp*ladb*sqrt(N)*(E^(8/9))*(debo^(11/18))/(sigb^(3/2))/A;
    ch4c=(2e-2)*fp*ladc*sqrt(N)*(E^(8/9))*(debo^(11/18))/(sigb^(3/2))/A;
    gi=izco*((T/1e2)^3);
    if ((m2>=0.2595) && (m2<=0.74))
        ch5=4.45*gi*r*(1+((1.41*A-1)^2)/2)/A;
        fa=ch123+ch4a+ch5-laef;
        fb=ch123+ch4b+ch5-laef;
        fc=ch123+ch4c+ch5-laef;        
    end
    if (m2>0.74)
        ch3=2.23/(1+B/H/dem4)/(1.41*A-1);
        ch6=((A^2)-1.58)/A/(1+B/H/dem5);
        ch5=4.45*gi*r*(1+(A^2)/pi)/A;
        fa=(ch1+ch2+ch3+ch6)*lam0/A+ch4a+ch5-laef;
        fb=(ch1+ch2+ch3+ch6)*lam0/A+ch4b+ch5-laef;
        fc=(ch1+ch2+ch3+ch6)*lam0/A+ch4c+ch5-laef;
    end
    la = VybCAB(fa,fb,fc,lada,ladb,ladc);
    lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
up=1;
m2=porex;
if (por<5e-1)
    m2=po;
end
if (por<28e-2) 
uo=opredUrovPodderM03(por); 
else uo=1e0/(1e0-m2);
end
T=tem(k);
laef=polyval(kektp,T);
lam0=opredKTPVozd(T); %температура в К
ladc=urovOtsechen(ladc,laef,1/m22);
ladc=urovPodder(ladc,lam0,1+m22);
ladi(k)=ladc;
end
lam=ladi;
end

function ak = opredKoefAkkomodLandau(T)
PP=6.6260755e-34/2/pi;
kB=1.380658e-23;
NA=6.0221409e23;
R=kB*NA;
mu=29e-3;
a=3e-9;
gamv=7/5;
m=mu/NA;
ro=opredPlotnVozd(T);
c=sqrt(gamv*R*T/mu);
alpha=(a*c)^2;
alpha=(kB*T)/alpha;
alpha=alpha^(3/2);
alpha=alpha/6/ro;
alpha=alpha/sqrt(2*pi*m);
%alpha=(kB*T/PP/c)^3;
%alpha=1.7*m*alpha/ro;
ak=alpha;
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

function opr = opredKTPVozd(T)
te0=273.15;
te = arrTempAir()+te0;
lamb = koefTeploprovAir();
n = length(te);
lam=opredKTPTKTochSha(lamb,te,T,n);
opr = lam;
end

function [ vy ] = VybCAB(fa,fb,fc,xa,xb,xc)
if (((fa*fc)<0) && ((fb*fc)>0))
    xb=xc; 
end
if (((fa*fc)>0) && ((fb*fc)<0) )
    xa=xc;
end
vy = [xa,xb,xc];
end

function fl = urovOtsechen(pro,ref,urot)
if (pro<0)
    f=0;
elseif (pro>(urot*ref))
    f=0;
else
    f=1;
end
if (f>0)
    fl=pro;
else
    fl=0;
end
end

function fl = urovPodder(pro,ref,urpo)
f=1; 
if (pro<0) 
    f=0; 
elseif (pro<(urpo*ref)) 
    f=0;
else f=1;
end
if (f>0) 
    fl=pro;
else fl=0;
end
urpo=fl;
end

function t = opredUrovPodderM03(por)
n=3; 
pam(3)=-8.93227577791347; 
pam(2)=8.89444783404518; 
pam(1)=1.06833435021354;
p=0; s=0; 
for k=1:n 
    s=s+pam(k)*(por^p); 
    p=p+1;
end
t=1/(1-s*por);
end