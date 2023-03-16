function [ vm ] = rasPorpoRazkvi(vypl,ide)
cvmi=6; f=cvmi; 
rpn=1e0; dp=1e0; k0=1e-2;
s=0; l=0; p=1e-6; h=p; r=0;
rp=poisrasprporkvi(vypl,1); 
raspr=poisrasprporkvi(vypl,2); 
prgr=poisrasprporkvi(vypl,3);
legr=poisrasprporkvi(vypl,4);
qg=length(legr);
r=0; srra=0; sp=0;
for k=1:qg
    p=prgr(k); 
    l=legr(k); 
    raspr(k)=raspr(k)*k0; 
    srra(k)=(p+l)/2e0; 
    s=s+srra(k)*raspr(k); 
    sp=sp+raspr(k);
    r=r+h; 
end
s=(s/sp+srRaPorkvi(vypl))/2e0;
if (ide==1)
    mu(1)=s; %средний размер пор
elseif (ide==2)
    mu(1)=r; %максимальный размер  пор
elseif (ide==3)
    mu=raspr; %распределение пор, %
elseif (ide==4)
    mu=rp; %распределение пор, см3
elseif (ide==5)
    mu=srra; %средний размер (массив)
elseif (ide==6)
    mu=legr; %левая граница
elseif (ide==7)
    mu=prgr; %правая граница
else
    mu=[0];
end
vm=mu;
end

function [ vm ] = poisrasprporkvi(vybmar,vybmas)
if (vybmar==4) 
    raspr0=rasprPorpoRazmAbskvi400(); 
    prgr0=LegrRasPorpoRazmkvi400(); 
elseif (vybmar==5) 
    raspr0=rasprPorpoRazmAbskvi500(); 
    prgr0=LegrRasPorpoRazmkvi500();
elseif (vybmar==6)
    raspr0=rasprPorpoRazmAbskvi600(); 
    prgr0=LegrRasPorpoRazmkvi600(); 
elseif (vybmar==7)
    raspr0=rasprPorpoRazmAbskvi700(); 
    prgr0=LegrRasPorpoRazmkvi700(); 
elseif (vybmar==8)
    raspr0=rasprPorpoRazmAbskvi800();
    prgr0=LegrRasPorpoRazmkvi800();
elseif (vybmar==9)
    raspr0=rasprPorpoRazmAbskvi900(); 
    prgr0=LegrRasPorpoRazmkvi900();
elseif (vybmar==10)
    raspr0=rasprPorpoRazmAbskvi1000();
    prgr0=LegrRasPorpoRazmkvi1000();
end
vm=obrabMaskvi(raspr0,prgr0,vybmas);
end

function [ vm ] = obrabMaskvi(raspr0,prgr0,vyb)
prk=300; nh=prk; n=nh; qg=n;
prmi=1; k0=1e2; ht=1e-6; 
n1=length(prgr0);
mpg0=max(prgr0); mrp0=min(raspr0); 
for k=1:n 
    ras01(k)=0; 
    ras02(k)=0;
end
k=1; legr0(k)=0; 
for k=2:n1
    legr0(k)=prgr0(k-1);
end
for k=1:qg
    prgr01(k)=k;
end
mpg01=max(prgr01); p=2; ri=1;
for k=2:qg
    if ((prgr0(k)>=prmi) && (ri>0))
        p=k; ri=-1; 
        break;
    end
end
pr=raspr0(p); w=raspr0(p-1); pr=pr-w;
s=prgr0(p); ko=prgr0(p-1); s=s-ko;
s=pr/s; r=w+s*(prmi-ko);
v=1; w0=raspr0(v); ras01(v)=w0; v=v+1; 
ras01(v)=r; pr=r; v=v+1; v0=v;
u=1; r=(w0-r)/w0; ras02(u)=r*k0; u=u+1;
for k=v0-1:qg %идем по prgr01, legr01
    f=1; x=-1;
        for q=p:n1 %идем по prgr0
        if ((legr0(q)<=prgr01(k)) && (f>0) && (prgr0(q)>=prgr01(k)))
			x=q; f=-1; break; 
        end
        end
if ((x>0) && (f<0))
s=raspr0(x); r=raspr0(x-1); s=s-r; 
r=prgr0(x); w=legr0(x); r=r-w; ko=s/r;
r=prgr01(k); w=legr0(x); r=r-w;
w=raspr0(x-1); s=w+ko*r;
else
r=mpg01-mpg0; w=-mrp0; r=w/r;
w=prgr01(k)-mpg0;
s=mrp0+r*w;
end
w=pr; w=(w-s)*k0/w0;
if (v<qg) 
    ras01(v)=s; v=v+1;
end
pr=s;
if (u<qg) 
    ras02(u)=w; u=u+1; 
end
end  
s=0; w=1e2; 
s=sum(ras02);
prgr01=ht*prgr01;
ras02=ras02*w/s;
k=1; legr01(k)=0; mu(k)=0;
for k=2:qg
    legr01(k)=prgr01(k-1); mu(k)=0;
end
eps=1e-7;
for k=2:qg
    if (ras02(k)<eps)
        ras02(k)=ras02(k-1);
    end
end
s=sum(ras02);
ras02=ras02*w/s;
if (vyb==1)
    mu=ras01; 
elseif (vyb==2)
    mu=ras02; 
elseif (vyb==3) 
    mu=prgr01;
elseif (vyb==4)
    mu=legr01;
else mu=[0];
end
vm=mu;
end

function rpr4 = rasprPorpoRazmAbskvi400()
rpr4=[12.936,12.851,12.667,12.5,12.0,10.320,10.0,9.125,8.0,6.980,6.0,4.917,4.0,3.125,2.0,1.423,1.115,0.846,0.615,0.471,0.392,0.314,0.275];
end

function rp4 = LegrRasPorpoRazmkvi400()
rp4=[0.053,0.631,0.774,1.0,1.442,5.374,6.61,10.0,13.092,16.408,20.347,26.64,32.407,43.481,59.785,71.225,84.243,100.0,128.792,149.908,193.07,213.634,261.567];
end

function opr5 = rasprPorpoRazmAbskvi500()
opr5=[10.414,10.339,10.207,10.0,9.822,9.663,9.465,9.212,9.059,8.870,8.664,8.458,8.335,8.19,8.0 ... 
7.695,7.477,7.172,6.859,6.547,6.105,6.0,4.929,4.469,4.0,3.716,3.406,2.916,2.699,2.408,2.168 ...
2.0,1.814,1.616,1.414,1.201,0.935,0.755,0.667,0.561,0.464,0.36,0.298,0.286,0.28,0.273];
end

function rpr5 = LegrRasPorpoRazmkvi500()
rpr5=[0.055,0.445,0.55,0.64,0.71,0.802,1.0,1.482,1.944,2.348,2.861,3.521,3.863,4.665,4.945,6.015,6.7,7.786 ...
8.873,10.0,11.08,11.256,11.981,13.127,14.429,15.924,18.065,22.529,25.618,30.105,33.878,37.381,41.736 ...
46.477,53.062,61.33,71.807,84.073,91.316,100.0,117.438,137.916,162.24,190.638,223.378,252.815];
end

function rpr6 = rasprPorpoRazmAbskvi600()
rpr6=[8.585,8.558,8.4,8.246,8.173,8.0,7.849,7.652,7.41,7.115,6.938,6.721,6.477,6.222,6.0,5.886,5.549,5.338 ...
5.014,4.593,4.283,4.0,3.715,3.178,2.822,2.362,2.233,2.0,1.747,1.531,1.302,1.093,0.883,0.739,0.617 ...
0.528,0.456,0.313,0.306,0.275,0.269];
end

function rp6 = LegrRasPorpoRazmkvi600()
rp6=[0.061,0.606,0.784,0.921,1.0,1.238,1.648,2.296,3.209,4.477,5.288,6.247,7.352,8.682 ...
9.633,10.321,11.711,12.674,14.260,16.050,17.383,18.013,18.738,19.882,21.687,24.608 ...
25.803,27.377,31.436,35.390,41.447,48.541,56.848,66.578,77.972,84.381,100.0,123.24,149.563,190.423,252.405];
end

function rpr7 = rasprPorpoRazmAbskvi700()
rpr7=[5.678,5.626,5.579,5.385,5.142,4.865,4.709,4.615,4.561,4.493,4.446,4.361,4.242,4.074,4.0 ... 
3.883,3.44,3.149,2.864,2.55,2.259,2.0,1.904,1.741,1.617,1.486,1.381,1.292,1.205,1.143,1.043,0.957 ...
0.876,0.807,0.725,0.663,0.594,0.544,0.469,0.406,0.384,0.347,0.285,0.173,0.149];
end

function rp7 = LegrRasPorpoRazmkvi700()
rp7=[0.06,1.0,1.105,1.452,1.865,2.394,2.947,3.479,4.109,4.852,5.73,6.766,7.99,9.435,10.0 ...
11.081,13.716,15.441,17.367,19.548,22.002,24.34,25.761,27.923,30.219,32.702,35.39,38.299 ...
41.447,44.854,48.541,52.531,56.848,61.521,66.578,72.050,77.972,84.381,91.316,100.0,107.542 ...
117.537,138.15,176.041,233.572];
end

function rpr8 = rasprPorpoRazmAbskvi800()
rpr8=[5.106,5.078,5.044,4.935,4.792,4.601,4.348,4.0,3.929,3.811,3.718,3.686,3.663,3.631,3.608,3.576 ...
3.545,3.467,3.357,3.224,3.055,2.977,2.852,2.779,2.669,2.565,2.468,2.383,2.279,2.156,2.084,2.0,1.846 ...
1.611,1.481,1.309,1.173,1.074,0.932,0.765,0.63,0.568,0.457,0.358,0.309,0.272,0.225,0.21,0.18,0.172 ...
0.165,0.157,0.15,0.142];
end

function rp8 = LegrRasPorpoRazmkvi800()
rp8=[0.061,0.348,0.417,0.471,0.554,0.648,0.762,0.952,1.0,1.161,1.348,1.489,2.216,2.704,2.986 ...
3.299,3.644,4.446,5.424,6.618,8.075,8.919,10.0,10.82,11.708,12.669,13.708,14.833,16.050,17.367 ...
18.792,19.242,22.002,25.761,27.874,30.161,32.636,35.314,41.346,48.41,56.679,61.33,71.807 ...
84.073,90.971,100.0,110.181,121.398,133.757,147.374,162.378,178.909,197.123,217.191];
end

function rpr9 = rasprPorpoRazmAbskvi900()
rpr9=[4.374,4.268,4.236,4.14,4.0,3.868,3.66,3.538,3.349,3.236,3.094,2.934,2.841,2.71,2.516 ...
2.302,2.151,2.0,1.778,1.591,1.404,1.218,1.076,0.932,0.765,0.63,0.568,0.457,0.358,0.309 ...
0.272,0.225,0.21,0.18,0.172,0.165,0.157,0.15,0.142];
end

function rp9 = LegrRasPorpoRazmkvi900()
rp9=[0.063,0.823,1.0,1.127,1.254,1.43,1.711,2.046,2.926,3.715,4.716,5.987,6.746,8.067,10.0 ...
12.546,14.053,15.759,18.69,22.166,26.288,31.176,34.931,41.346,48.41,56.679,61.33,71.807,84.073 ...
90.971,100.0,110.181,121.398,133.757,147.374,162.378,178.909,197.123,217.191];
end

function rpr10 = rasprPorpoRazmAbskvi1000()
rpr10=[3.607,3.549,3.446,3.066,2.616,2.151,2.047,1.663,1.328 ...
1.092,0.893,0.713,0.515,0.353,0.25,0.14,0.059];
end

function rp10 = LegrRasPorpoRazmkvi1000()
rp10=[0.063,0.820,1.0,2.046,4.186,8.469,10.0,14.609,19.412,23.462,28.358,34.275,41.427 ...
50.071,60.519,73.147,88.410];
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

function srp = srRaPorkvi(vybmar)
prgr=poisrasprporkvi(vybmar, 3); 
legr=poisrasprporkvi(vybmar, 4);
raspr=poisrasprporkvi(vybmar, 2);
sp=sum(raspr);
srra=(legr+prgr)/2e0; %средний размер пор
s=0;
n=(length(srra)+length(raspr))/2;
for k=1:n
    s=s+srra(k)*raspr(k);
end
srp=s/sp;
end