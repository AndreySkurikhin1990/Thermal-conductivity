function [ stch ] = oprStCherkvi(vmkvi,etek,vyma,fl)
ko=1e-2; tnd = 6e2; dtd = 2e2; dkosckl=6; dkosckt(1) = tnd;
for k=2:dkosckl
        dkosckt(k) = dkosckt(k - 1) + dtd;
end
dkosckm=opdopkoefoslstch(vmkvi,dkosckl);
tnosck=3e2; dtosck=1e2; dmkoosck=14;
kusck=oprmassialmgkvi(vmkvi);
wal=kusck(1); wmg=kusck(2); wsi=kusck(3);
tkusck(1) = tnosck; kusck=0;
for k = 2:dmkoosck
    tkusck(k) = tkusck(k - 1) + dtosck;
end
kusck = koefoslab(wmg, wsi, wal, tkusck);
tnack = 2e2; tek0 = 273.15; tnak = tnack + tek0; etek(1) = tnak; detek = 1e2; cemk = 11;
for j = 2:cemk
    etek(j) = etek(j - 1) + detek;
end
if (fl==1)
for k=1:cemk
        s = epsisredkvi(etek(k), tkusck, kusck, dmkoosck, dkosckt, dkosckm, dkosckl, vmkvi);
        stchsrkvi(k) = s; 
end
end
switch (vyma)
    case (1)
    stch=stchsrkvi';
    case (2)
    stch=[wsi, wal, wmg];
    case (3)
    stch=kusck';
    case (4)
    stch=tkusck';
    case (5)
    stch=dkosckm';
    case (6)
    stch=dkosckt';
    otherwise 
        stch=[0];
end
end

function [ om ] = oprmassialmgkvi(vmkvi)
if (vmkvi == 4)
		salok = 33e0; smgok = 15e0; ssiok = 52e0; 
		wal = 25e0; wsi = 11e0; wmg = 4e1;
elseif (vmkvi == 5) 
		salok = 34e0; smgok = 11e0; ssiok = 54e0;
		wal = 28e0; wmg = 8e0; wsi = 44e0;
elseif (vmkvi == 6) 
		salok = 36e0; smgok = 9e0; ssiok = 55e0;
		wal = 3e1; wsi = 7e0; wmg = 45e0; 
elseif (vmkvi == 7)
		salok = 37e0; smgok = 8e0; ssiok = 55e0; 
		wal = 31e0; wmg = 6e0; wsi = 45e0; 
elseif (vmkvi == 8)
		salok = 38e0; smgok = 7e0; ssiok = 55e0; 
		wal = 3e1; wmg = 5e0; wsi = 45e0;
elseif (vmkvi == 9)
		salok = 39e0; smgok = 6e0; ssiok = 55e0; 
		wal = 32e0; wmg = 5e0; wsi = 45e0;
elseif (vmkvi == 10) 
		salok = 39e0; smgok = 6e0; ssiok = 55e0; 
		wal = 32e0; wmg = 4e0; wsi = 45e0;
else
        disp('Net takoy marki KVI!');
end
ko=1e-2;
om=[wal,wmg,wsi,salok,smgok,ssiok]*ko;
end

function [ om ] = opdopkoefoslstch(vmkvi,dkosckl)
dkosckm=0;
k = 1; 
        dkosckm(k) = 3e0; k=k+1; 
        dkosckm(k) = 6e0; k=k+1; 
        dkosckm(k) = 7e0; k=k+1; 
		dkosckm(k) = 15e0; k=k+1; 
        dkosckm(k) = 2e1; k=k+1; 
        dkosckm(k) = 26e0;
if (vmkvi>6)
            k=2; dkosckm(k)=7e0;
end
if (vmkvi>9)
            k=3; dkosckm(k)=8e0;
end
ko=1e-2;
for k=1:dkosckl
		tm = dkosckm(k)*ko; 
        dkosckm(k) = 1e0 - tm;
end
om=dkosckm;
end

function [ t ] = oslabintegrstepchernalsimg(vy,te)
format long g;
n=length(te);
scma=[1,0.940,0.870,0.801,0.736,0.676,0.635,0.590,0.567,0.543,0.530,0.525,0.515,0.507]; %степень черноты магнезита
mdma=[0.98,0.02]; %MgO - магнезит
%-------
scsh=[1,0.976,0.949,0.905,0.859,0.812,0.774,0.737,0.709,0.681,0.661,0.639,0.626,0.620]; %степень черноты шамота
scksh=[1,0.980,0.951,0.920,0.883,0.853,0.821,0.790,0.767,0.746,0.730,0.715,0.705,0.692]; %степень черноты корундошамота
scktik=[1,0.983,0.936,0.867,0.819,0.721,0.659,0.593,0.541,0.490,0.453,0.429,0.403,0.384]; %степень черноты каолинового теплоизоляционного кирпича
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

function es = epsisredkvi(T, tdkusctm, dkusctm, dmkoosck, dkoscet, dkoscem, dkosckl, vk)
fileID=fopen('DlinaVolny_kvi.txt','r'); 
formatSpec='%f';
dv=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=Kramers_Kronig_kvi(vk); 
dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, length(tdkusctm));
dkusct=ProvAdek(dkusct);
dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, length(dkoscet));
dkosce=ProvAdek(dkosce);
ume=dkusct*dkosce;
epsil=epsilnu(npp);
epssr=usrednen(T, epsil, dv, npp);
es=ume*epssr;
end

function [ np ] = Kramers_Kronig_kvi(vk)
if (vk==4)
    np=Kramers_n_kvi400(); 
elseif (vk==5)
    np=Kramers_n_kvi500(); 
elseif (vk==6)
    np=Kramers_n_kvi600(); 
elseif (vk==7)
    np=Kramers_n_kvi700(); 
elseif (vk==8)
    np=Kramers_n_kvi800(); 
elseif (vk==9)
    np=Kramers_n_kvi900(); 
elseif (vk==10)
    np=Kramers_n_kvi1000(); 
else
    np=[0];
end
end

function [ e ] = epsilnu(npp)
lear=length(npp);
for k=1:lear
	n=abs(npp(k));
	eps=(4e0*n+2e0)/3e0/((n+1e0)^2e0);
	eps=eps+2e0*(n^3e0)*((n^2e0)+2e0*n-1e0)/((n^2e0)+1e0)/((n^4e0)-1e0);
	eps=eps-8e0*(n^4e0)*((n^4e0)+1e0)*log(n)/((n^2e0)+1e0)/(((n^4e0)-1e0)^2e0);
	eps=eps-(n^2e0)*log((n-1e0)/(n+1e0))*(((n^2e0)-1e0)^2e0)/(((n^2.0)+1.0)^3.0); 
    ep(k)=eps;
end
e=ep;
end

function knus = usrednen(tem,knuSha,dv,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*(c0^2);
c2=PP*c0/PB;
ct1=0; ct2=0; dl=0;
for k=1:length(npp)
    ct1=ct1/(npp(k)^2);
    ct2=ct2/npp(k);
    lambda=dv(k)/npp(k);
    dl(k)=lambda;
    Ib(k)=ct1/((lambda^5)*(exp(ct2/(lambda*tem))-1));
    Ib(k)=ProvAd(Ib(k));
    Ibc(k)=knuSha(k)*Ib(k);
ct1=c1;
ct2=c2;
end
nc=trapz(dl,Ibc);
nz=trapz(dl,Ib);
knus=nc/nz;
end

function mk = ProvAd(m)
ep=1e-30;
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
if (abs(m)<ep)  
    m=0;  
end; 
mk=real(abs(m));
end