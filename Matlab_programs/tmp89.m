%расчет доли пор  ¬» разных марок по заданным диапазонам
function t = tmp89()
format shortg;
%t = raschkvi();
%t = epssrver();
t=opredcp();
end

function otv = opredcp()
kk=1e-2; k3=1e-3; te0=273.15;
mualo=101.96*k3; mucao=56.08*k3; mucro=152*k3; mufeo3=159.7*k3;
mufeo=69.06*k3; mumgo=40.31*k3; musio=60.09*k3; muzro=123.22*k3; muo=15.999*k3;
muma=[musio, mualo, mufeo3, mucao, mumgo, muo];
kcpalo=[154.96, -16.168, 7.120, -20.817];
kcpcao=[57.68, -1.324, 1.560, -4.418];
kcpfeo3=[176.60, -24.0, 7.2, -19.843];
kcpfeo=[45.31, 14.9, -2.6, -0.374];
kcpmgo=[63.24, -7.632, 2.88, -7.263];
kcpsio=[77.09, 3.384, -0.16, -10.558];
%sio=53.94*kk; alo=42.68*kk; feo3=2.29*kk; cao=0.65*kk; mgo=0.44*kk; ro=400; lam=0.2; lam=0.24; %shamot
%sio=50.51*kk; alo=23.97*kk; feo3=2.88*kk; cao=3.6*kk; mgo=19.03*kk; ro=440; lam=0.16; %lam=0.12; %itom-440
%sio=0.03*kk; alo=99.54*kk; feo3=0.03*kk; cao=0.03*kk; mgo=1-sio-alo-feo3-cao; ro=1.3e3; lam=0.8; %kt13
%sio=(74.8+88.15)/2*kk; alo=(3.34+9.75)/2*kk; feo3=(2.37+5.26)/2*kk; cao=(0.47+0.85)/2*kk; mgo=(0.61+1.71)/2*kk; ro=4e2; lam=0.15; %lam=0.14; %diatomit
sio=36.8*kk; alo=13.2*kk; feo3=8.3*kk; cao=2*kk; mgo=21.7*kk; %ro=25e1; 
ro=1e2; %lam=0.23; 
lam=0.16; %vermiculit
%sio=43*kk; alo=54*kk; feo3=1.8*kk; cao=1*kk; mgo=1e-8*kk; ro=340; lam=0.16; %lam=0.23; %mkrv-340
%sio=92*kk; alo=2.5*kk; feo3=3*kk; cao=1.5*kk; mgo=0.5*kk; ro=7e2; lam=0.34; lam=0.39; %dinas
obso=sio+alo+feo3+cao+mgo; 
sio=sio/obso; alo=alo/obso; feo3=feo3/obso; cao=cao/obso; mgo=mgo/obso;
%ro_sio=2.29; ro_alo=3.50; ro_mgo=3.49; ro_feo=2.24; ro_cao=3.35;
sodox=[sio, alo, feo3, cao, mgo]/kk
%sodox=[sio/ro_sio, alo/ro_alo, feo3/ro_feo, cao/ro_cao, mgo/ro_mgo];
%ssodox=sum(sodox);
%sodox=sodox/ssodox
T=35+te0;
%T=350+te0
%T=600+te0
%k=1;
%tem=[301, 367, 466, 615, 765, 914, 1064, 1264];
%T=tem(k)
%lamb=[0.070, 0.080, 0.103, 0.150, 0.207, 0.283, 0.373, 0.477]; lam=lamb(k);
%t1=620;	l1=0.180; t2=891; l2=0.245; lam=l1+(l2-l1)/(t2-t1)*(T-t1)
ep=1e-3;
T=ep*T;
A=[1, 1, 1, 1, 1, 1; 
    musio*(1-1/sio), mualo, mufeo3, mucao, mumgo, muo; 
    musio, mualo*(1-1/alo), mufeo3, mucao, mumgo, muo; 
    musio, mualo, mufeo3*(1-1/feo3), mucao, mumgo, muo; 
    musio, mualo, mufeo3, mucao*(1-1/cao), mumgo, muo;
    musio, mualo, mufeo3, mucao, mumgo*(1-1/mgo), muo];
b=[1, 0, 0, 0, 0, 0];
x=inv(A)*b';
eps=1e-6;
muob=0;
su=0;
for k=1:length(x)
    if (x(k)<eps)
        x(k)=0;
    end
    muob=muob+x(k)*muma(k);
    su=su+x(k);
end
x=x'/su
muob=muob';
vsio=kcpsio(1)+kcpsio(2)*T+kcpsio(3)*(T^2)+kcpsio(4)/T;
valo=kcpalo(1)+kcpalo(2)*T+kcpalo(3)*(T^2)+kcpalo(4)/T;
vfeo3=kcpfeo3(1)+kcpfeo3(2)*T+kcpfeo3(3)*(T^2)+kcpfeo3(4)/T;
vcao=kcpcao(1)+kcpcao(2)*T+kcpcao(3)*(T^2)+kcpcao(4)/T;
vmgo=kcpmgo(1)+kcpmgo(2)*T+kcpmgo(3)*(T^2)+kcpmgo(4)/T;
a=x(1); b=x(2); c=x(3); d=x(4); e=x(5); f=x(6); 
vsio=vsio*a;
valo=valo*b;
vfeo3=vfeo3*c;
vcao=vcao*d;
vmgo=vmgo*e;
cp=(vsio+valo+vfeo3+vcao+vmgo)/muob
tempprov=lam/cp/ro
%cp=1143; %shamot
%cp=907; %dinas 350
%cp=968; %dinas 600
%bsr=sqrt(lam*ro*cp)
otv=0;
end

function t = raschkvi()
vybkvi=10;
dkvi400=[0, 0.66, 2.08, 3.37, 7.24, 20.22, 22.70, 29.46, 38.16, 46.05 ...
53.62, 61.99, 69.08, 75.84, 84.54, 89.00, 91.38, 93.46, 95.24, 96.36 ...
96.97, 97.57, 97.88];
rkvi400=[0.053, 0.631, 0.774, 1.000, 1.442, 5.374, 6.610, 10.000 ...
13.092, 16.408, 20.347, 26.640, 32.407, 43.481, 59.785, 71.225 ...
84.243, 100.000, 128.792, 149.908, 193.070, 213.634, 261.567];
dkvi500=[0, 0.72, 1.99, 3.97, 5.68, 7.21, 9.11, 11.54, 13.01 ...
14.83, 16.80, 18.78, 19.97, 21.35, 23.18, 26.10, 28.21, 31.13 ...
34.13, 37.13, 41.38, 42.38, 52.67, 57.09, 61.59, 64.32, 67.29 ...
72.00, 74.08, 76.88, 79.18, 80.79, 82.58, 84.48, 86.43, 88.46 ...
91.02, 92.75, 93.60, 94.62, 95.54, 96.54, 97.14, 97.26, 97.31, 97.38];
rkvi500=[0.055, 0.445, 0.550, 0.640, 0.710, 0.802, 1.000, 1.482 ...
1.944, 2.348, 2.861, 3.521, 3.863, 4.665, 4.945, 6.015, 6.700 ...
7.786, 8.873, 10.000, 11.080, 11.256, 11.981, 13.127, 14.429 ...
15.924, 18.065, 22.529, 25.618, 30.105, 33.878, 37.381, 41.736 ...
46.477, 53.062, 61.330, 71.807, 84.073, 91.316, 100.000, 117.438 ...
137.916, 162.240, 190.638, 223.378, 252.815];
rkvi600=[0.061, 0.606, 0.784, 0.921, 1.000, 1.238, 1.648, 2.296 ...
3.209, 4.477, 5.288, 6.247, 7.352, 8.682, 9.633, 10.321, 11.711 ...
12.674, 14.260, 16.050, 17.383, 18.013, 18.738, 19.882, 21.687 ...
24.608, 25.803, 27.377, 31.436, 35.390, 41.447, 48.541, 56.848 ...
66.578, 77.972, 84.381, 100.000, 123.240, 149.563, 190.423, 252.405];
dkvi600=[0, 0.32, 2.16, 3.95, 4.80, 6.82, 8.58, 10.87, 13.69, 17.13 ...
19.19, 21.71, 24.56, 27.53, 30.11, 31.45, 35.37, 37.83, 41.60 ...
46.51, 50.11, 53.41, 56.73, 62.98, 67.13, 72.48, 73.99, 76.70, 79.65, 82.17 ...
84.83, 87.27, 89.72, 91.39, 92.81, 93.85, 94.69, 96.36, 96.43, 96.80, 96.87];
dkvi700=[0, 0.92, 1.74, 5.16, 9.44, 14.32, 17.07, 18.72, 19.67 ...
20.87, 21.70, 23.19, 25.29, 28.25, 29.55, 31.61, 39.42, 44.54 ...
49.56, 55.09, 60.21, 64.78, 66.47, 69.34, 71.52, 73.83, 75.68 ...
77.25, 78.78, 79.87, 81.63, 83.15, 84.57, 85.79, 87.23, 88.32 ...
89.54, 90.42, 91.74, 92.85, 93.24, 93.89, 94.98, 96.95, 97.38];
rkvi700=[0.060, 1.000, 1.105, 1.452, 1.865, 2.394, 2.947, 3.479 ...
4.109, 4.852, 5.730, 6.766, 7.990, 9.435, 10.000, 11.081, 13.716 ...
15.441, 17.367, 19.548, 22.002, 24.340, 25.761, 27.923, 30.219 ...
32.702, 35.390, 38.299, 41.447, 44.854, 48.541, 52.531, 56.848 ...
61.521, 66.578, 72.050, 77.972, 84.381, 91.316, 100.000, 107.542 ...
117.537, 138.150, 176.041, 233.572];
dkvi800=[0, 0.55, 1.21, 3.35, 6.15, 9.89, 14.85, 21.66, 23.05 ...
25.36, 27.18, 27.81, 28.26, 28.89, 29.34, 29.96, 30.57, 32.10 ...
34.25, 36.86, 40.17, 41.70, 44.14, 45.57, 47.73, 49.76, 51.66 ...
53.33, 55.37, 57.78, 59.19, 60.83, 63.85, 68.45, 70.99, 74.36 ...
77.03, 78.97, 81.75, 85.02, 87.66, 88.88, 91.05, 92.99, 93.95 ...
94.67, 95.59, 95.89, 96.47, 96.63, 96.77, 96.93, 97.06, 97.22];
rkvi800=[0.061, 0.348, 0.417, 0.471, 0.554, 0.648, 0.762, 0.952 ...
1.000, 1.161, 1.348, 1.489, 2.216, 2.704, 2.986, 3.299, 3.644 ...
4.446, 5.424, 6.618, 8.075, 8.919, 10.000, 10.820, 11.708, 12.669 ...
13.708, 14.833, 16.050, 17.367, 18.792, 19.242, 22.002, 25.761 ...
27.874, 30.161, 32.636, 35.314, 41.346, 48.410, 56.679 ...
61.330, 71.807, 84.073, 90.971, 100.000, 110.181, 121.398 ...
133.757, 147.374, 162.378, 178.909, 197.123, 217.191];
dkvi900=[0, 2.42, 3.16, 5.35, 8.55, 11.57, 16.32, 19.11 ...
23.43, 26.02, 29.26, 32.92, 35.05, 38.04, 42.48, 47.37 ...
50.82, 54.28, 59.35, 63.63, 67.90, 72.15, 75.40, 78.69 ...
82.51, 85.60, 87.01, 89.55, 91.82, 92.94, 93.78, 94.86 ...
95.20, 95.88, 96.07, 96.23, 96.41, 96.57, 96.75];
rkvi900=[0.063, 0.823, 1.000, 1.127, 1.254, 1.430 ...
1.711, 2.046, 2.926, 3.715, 4.716, 5.987, 6.746, 8.067 ...
10.000, 12.546, 14.053, 15.759, 18.690, 22.166, 26.288 ...
31.176, 34.931, 41.346, 48.410, 56.679, 61.330, 71.807 ...
84.073, 90.971, 100.000, 110.181, 121.398, 133.757 ...
147.374, 162.378, 178.909, 197.123, 217.191];
dkvi1000=[0, 1.61, 4.46, 15.00, 27.47, 40.37, 43.25, 53.90 ...
63.18, 69.73, 75.24, 80.23, 85.72, 90.21, 93.07, 96.12, 98.36, 100.00];
rkvi1000=[0.063, 0.820, 1.000, 2.046, 4.186, 8.469, 10.000, 14.609 ...
19.412, 23.462, 28.358, 34.275, 41.427, 50.071, 60.519, 73.147, 88.410, 235.429];
x=[150, 95, 85, 75, 65, 55, 45, 35, 25, 15, 7.5, 4, 2, 0.75, 0.3, 0.055];
m=length(x);
switch (vybkvi)
    case (4)
        xr=rkvi400;
        yr=dkvi400;
    case (5)
        xr=rkvi500;
        yr=dkvi500;
    case (6)
        xr=rkvi600;
        yr=dkvi600;
    case (7)
        xr=rkvi700;
        yr=dkvi700;
    case (8)
        xr=rkvi800;
        yr=dkvi800;
    case (9)
        xr=rkvi900;
        yr=dkvi900;
    case (10)
        xr=rkvi1000;
        yr=dkvi1000;
end
n=(length(xr)+length(yr))/2;
for k=1:m
    y(k)=opredKTPTKTochSha(yr, xr, x(k), n);
end
y=y';
y=round(y)
t=0;
end

function ktpo = opredKTPTKTochSha(ktptks, te, temp, n)
f = 1; p = 0; ep=1e-10; ktp = 0;
if ((temp>te(1)) && (temp<te(n)))
for k = 1:n
if ((te(k) >= temp) && (f>0))
        p = k; f = 0; 
        break;
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
            if ((p==n) && (temp>te(p)))
                b=y2;
            elseif ((p==n) && (temp<=te(p)))
                b=y1;
            end
			ktp = b + ko*(temp - x1);
else 
    ktp=0;
end
end
if (ktp<0)
    ktp=0;
end
ktpo=ktp;
end

function e = epssrver()
kk=1e-3; wsi=368*kk; wal=132*kk; wmg=217*kk; vyfv=0
if (vyfv==0) %por = 0.6635; %пористость вермикулита фракции 2-0,7 мм
    wal=140*kk; wsi=458*kk; wmg=290*kk;
elseif (vyfv==1) %por=0.5575; %пористость вермикулита фракции 8-4 мм
    wal=129*kk; wsi=467*kk; wmg=282*kk;
end
wo=wsi+wal+wmg;
dkoscvm=oprdopkoefoslabstchver(1,wal/wo,wmg/wo,wsi/wo);
tdkoscvm=oprdopkoefoslabstchver(2,wal/wo,wmg/wo,wsi/wo);
kuscvm=oprdopkoefoslabstchver(3,wal/wo,wmg/wo,wsi/wo);
tkuscvm=oprdopkoefoslabstchver(4,wal/wo,wmg/wo,wsi/wo);
te0=273.15; te = 2e2:1e2:12e2; te=te+te0;
lete=length(te);
stchm=zeros(1, lete);
for k=1:lete
    ts=te(k);
dkusct=opredKTPTKTochSha(dkoscvm, tdkoscvm, ts, length(tdkoscvm));
dkusct=ProvAdek(dkusct);
dkosce=opredKTPTKTochSha(kuscvm, tkuscvm, ts, length(tkuscvm));
dkosce=ProvAdek(dkosce);
ume=dkusct*dkosce;
eps=epsisred(ts)*ume;
stchm(k)=eps;
end
stchm=stchm'
te=te'
e=0;
end

function r = ProvAdek(m)
if (m>1e0) 
    m=1e0; 
end
if (m<0) 
    m=0;
end
r=m;
end

function es = epsisred(T)
npp=Kramers_n_uk();
dv=dliny_voln();
lenp=length(npp);
eps=zeros(1, lenp);
iz0=zeros(1, lenp);
iz1=zeros(1, lenp);
dl=zeros(1, lenp);
for k=1:lenp
    eps(k)=epsilnu(npp(k));
end
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*(c0^2);
c2=PP*c0/PB;
c1t=c1; c2t=c2;
for j=1:lenp
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
iz0(j)=c1t/((la^5)*(exp(c2t/(la*T))-1));
iz1(j)=eps(j)*iz0(j);
c1t=c1;
c2t=c2;
end
chi=trapz(dl, iz1);
zna=trapz(dl, iz0);
es=chi/zna;
end

function [ dlv ] = dliny_voln()
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k); 
end
dlv = dl;
end

function [ stch ] = oprdopkoefoslabstchver(vyma,wal,wmg,wsi)
tnd = 6e2; dtd = 2e2; dkosckl=6; dkosckt(1) = tnd;
for k=2:dkosckl
        dkosckt(k) = dkosckt(k - 1) + dtd;
end
dkosckm=opdopkoefoslstch(dkosckl);
tnosck=3e2; dtosck=1e2; dmkoosck=14;
tkusck(1) = tnosck;
for k = 2:dmkoosck
    tkusck(k) = tkusck(k - 1) + dtosck;
end
tnack = 2e2; tek0 = 273.15; tnak = tnack + tek0; 
etek(1) = tnak; detek = 1e2; cemk = 11;
for k = 2:cemk
    etek(k) = etek(k - 1) + detek;
end
kusck = koefoslab(wmg, wsi, wal, tkusck);
if (vyma==1)
    stch=dkosckm';
elseif (vyma==2)
    stch=dkosckt';
elseif (vyma==3)
    stch=kusck';
elseif (vyma==4)
    stch=tkusck';
end
end

function [ om ] = opdopkoefoslstch(dkosckl)
k = 1; 
dkoscvm(k) = 4.68; k=k+1;
dkoscvm(k) = 4.69; k=k+1; 
dkoscvm(k) = 5.65; k=k+1; 
dkoscvm(k) = 13.17; k=k+1; 
dkoscvm(k) = 20.2; k=k+1; 
dkoscvm(k) = 27.81;
ko=1e-2;
dkosckl=(length(dkoscvm)+dkosckl)/2;
for k=1:dkosckl
		tm = dkoscvm(k)*ko; 
        dkoscvm(k) = 1e0 - tm;
end
om=dkoscvm;
end

function koef = koefoslab(wmg, wsi, wal, tere)
mgo=oslabintegrstepchernalsimg(0,tere);
alo=oslabintegrstepchernalsimg(2,tere);
sio=oslabintegrstepchernalsimg(1,tere);
wo=wmg+wsi+wal;
for k=1:length(tere)
kuo(k)=(mgo(k)*wmg+alo(k)*wal+sio(k)*wsi);
end
koef=kuo/wo;
end

function [ t ] = oslabintegrstepchernalsimg(vy,te)
format long g;
n=length(te);
scma=[1,0.940,0.870,0.801,0.736,0.676,0.635,0.590,0.567,0.543,0.530,0.525,0.515,0.507]; %степень черноты магнезита
mdma=[0.98,0.02]; %MgO - магнезит
%-------
scsh=[1,0.976,0.949,0.905,0.859,0.812,0.774,0.737,0.709,0.681,0.661,0.639,0.626,0.620]; %степень черноты шамота
scksh=[1,0.980,0.951,0.920,0.883,0.853,0.821,0.790,0.767,0.746,0.730,0.715,0.705,0.692]; %степень черноты корундошамота
scktik=[1,0.983,0.936,0.867,0.819,0.721,0.659,0.593,0.541,0.490,0.453,0.429,0.403,0.384]; %степень черноты каолинового теплоизол€ционного кирпича
scmu=[1,0.984,0.941,0.882,0.813,0.751,0.695,0.641,0.594,0.558,0.530,0.499,0.479,0.462]; %степень черноты муллита
sckr=[1,0.984,0.953,0.917,0.854,0.808,0.756,0.711,0.578,0.523,0.495,0.468,0.448,0.429]; %степень черноты кремнезема
%-------
mdsh=[0.56,0.396]; mdsh=[mdsh(1)/(mdsh(1)+mdsh(2)),mdsh(2)/(mdsh(1)+mdsh(2))]; mdsh=mdsh'; %SiO2, Al2O3 - шамот
mdksh=[0.28,0.70]; mdksh=[mdksh(1)/(mdksh(1)+mdksh(2)),mdksh(2)/(mdksh(1)+mdksh(2))]; mdksh=mdksh'; %SiO2, Al2O3 - корундошамот
mdktik=[0.57,0.4]; mdktik=[mdktik(1)/(mdktik(1)+mdktik(2)),mdktik(2)/(mdktik(1)+mdktik(2))]; mdktik=mdktik'; %SiO2, Al2O3 - каол. т/и кирпич
mdmu=[0.28,0.72]; %SiO2, Al2O3 - муллит
mdkr=[0.985,0.01]; mdkr=[mdkr(1)/(mdkr(1)+mdkr(2)),mdkr(2)/(mdkr(1)+mdkr(2))]; mdkr=mdkr'; %SiO2, Al2O3 - кремнезем
for j=1:n
kor=0; k=1;
kor=RasKor(mdsh(1),mdsh(2),mdksh(1),mdksh(2),scsh(j),scksh(j),kor,k); k=k+1;
kor=RasKor(mdsh(1),mdsh(2),mdktik(1),mdktik(2),scsh(j),scktik(j),kor,k); k=k+1;
kor=RasKor(mdsh(1),mdsh(2),mdmu(1),mdmu(2),scsh(j),scmu(j),kor,k); k=k+1;
kor=RasKor(mdksh(1),mdksh(2),mdktik(1),mdktik(2),scksh(j),scktik(j),kor,k); k=k+1;
kor=RasKor(mdksh(1),mdksh(2),mdmu(1),mdmu(2),scksh(j),scmu(j),kor,k); k=k+1;
kor=RasKor(mdktik(1),mdktik(2),mdmu(1),mdmu(2),scktik(j),scmu(j),kor,k); k=k+1;
kor=RasKor(mdsh(1),mdsh(2),mdkr(1),mdkr(2),scsh(j),sckr(j),kor,k); k=k+1;
kor=RasKor(mdksh(1),mdksh(2),mdkr(1),mdkr(2),scksh(j),sckr(j),kor,k); k=k+1;
kor=RasKor(mdktik(1),mdktik(2),mdkr(1),mdkr(2),scktik(j),sckr(j),kor,k); k=k+1;
kor=RasKor(mdmu(1),mdmu(2),mdkr(1),mdkr(2),scmu(j),sckr(j),kor,k); k=k+1;
usr=UsredMas(kor,k);
msio(j)=usr(1);
malo(j)=usr(2);
end
msio=msio';
malo=malo';
for k=1:length(scma)
    scma(k)=(scma(k)-msio(k)*mdma(2))/mdma(1);
end
scma=scma';
if (vy==0)
    t=scma;
elseif (vy==1)
    t=msio;
elseif (vy==2)
    t=malo;
end
end

function [ ko ] = RasKor(a11,a12,a21,a22, b1,b2,kor,k)
mas(1,1)=a11; mas(1,2)=a12; 
mas(2,1)=a21; mas(2,2)=a22; 
b(1)=b1; b(2)=b2; 
x=inv(mas)*b';
if ((x(1)>=0) && (x(2)>=0) && (x(1)<=1) && (x(2)<=1))
    kor(k,1)=x(1); kor(k,2)=x(2);
else
    kor(k,1)=0; kor(k,2)=0;
end
ko=kor;
end

function [ usr ] = UsredMas(kor,k)
s1=0; s2=0; q=0;
for j=1:k-1
    if ((kor(j,1)>0) && (kor(j,2)>0))
        s1=s1+kor(j,1); 
        s2=s2+kor(j,2); 
        q=q+1;
    end
end
s1=s1/q;
s2=s2/q;
usr=[s1,s2];
end