function [ oktptks ] = opredKTPTverKarkitom(vmitom, tem, y0)
hf=1e0;
k=0; hps=rasPorpoRazITOM(vmitom, k);
k=3; por=hps(k);
k=3; rp=rasPorpoRazITOM(vmitom, k);
k=5; sr=rasPorpoRazITOM(vmitom, k);
k=6; legr=rasPorpoRazITOM(vmitom, k);
k=7; prgr=rasPorpoRazITOM(vmitom, k);
k=1; hps=rasPorpoRazITOM(vmitom, k);
k=1; srp=hps(k); %y0=30e-3; %qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); 
lentem=length(tem);
x0=srp/(por^(1/3))/2;
srk=zeros(1,lentem);
    k=0; tholv = napMasEKTPitomNac(vmitom, k, y0, tem);
	k=k+1; tgorv = napMasEKTPitomNac(vmitom, k, y0, tem);
	k=k+1; qobv = napMasEKTPitomNac(vmitom, k, y0, tem);
	k=k+1; laefm = napMasEKTPitomNac(vmitom, k, y0, tem)
	k=k+1; temvs = napMasEKTPitomNac(vmitom, k, y0, tem)
dkosc=oprmassialmgitom(vmitom,0);
dkosct=oprmassialmgitom(vmitom,1);
stchi=oprmassialmgitom(vmitom,2);
wsio=stchi(1); walo=stchi(2); wmgo=stchi(3);
kusci=oprStCheritom(vmitom, temvs, 3, 0, dkosc, dkosct, walo, wmgo, wsio);
tkusci=oprStCheritom(vmitom, temvs, 4, 0, dkosc, dkosct, walo, wmgo, wsio);
stchi=oprStCheritom(vmitom, temvs, 1, 1, dkosc, dkosct, walo, wmgo, wsio)
%stchi=[0.774201168152912, 0.726556391091177,  0.679236566513966,  0.629831873704107, 0.586608357104253 ...
%0.54493572669103, 0.498913352639034, 0.454633399441586, 0.420445232169707, 0.393504043083303, 0.331449013299856];
for k=1:lentem
    lavo(k)=opredTeploprovVozd(tem(k));
end
na=1; cvp=na+21; cve=cvp+7*2; cved=cve+1; vyb=0; e=1e-6; q=1;
for k=na:cved
    if (k==na)
        srk=DulnevZernitom(por, temvs, laefm, x0, lavo, vyb, stchi);
    elseif (k<=cvp)
        srk = opredTvChaKoeTepSrInSpSha(k-na, por, temvs, laefm, x0, lavo);
    elseif (k<=cve) 
        srk=DulnevKoefTepVermN(k-cvp-1, temvs, laefm, srp, abs(cve-cvp)/2, rp, legr, prgr, sr, por, lavo);
    elseif (k<=cved) 
        srk=VasFrayObsh(laefm, lavo, por);
    end
    srk=proverka(srk, laefm);
    if (sum(srk)>e)
        k=k'
        srk=srk
        no(q)=k;
        for j=1:length(srk)
        srktp(q,j)=srk(j);
        end
        q=q+1;
    end
end
srktp=srktp
n=length(temvs);
srk=zeros(1,n);
no=no
for j=1:n
    s=0; t=0;
for k=1:length(no)
    if (srktp(k,j)>e)
        s=s+srktp(k,j);
        t=t+hf;
    end
end
if (t>e)
srk(j)=s/t;
end
end
srk=srk
oktptks=srk;
end

function [ om ] = oprmassialmgitom(vybitom, vvm)
ko=1e-2; tnd = 6e2; dtd = 2e2; dkoscil = 6;
pori440 = (8e1+82e0)*ko/2e0;
pori620 = (75e0+78e0)*ko/2e0;
pori860 = (65e0+68e0)*ko/2e0;
pori1000 = (62e0+65e0)*ko/2e0;
dkoscit(1) = tnd; 
for k = 2:dkoscil
    dkoscit(k) = dkoscit(k - 1) + dtd;
end
switch (vybitom)
    case (0)
		saloi = 26e0*ko; smgoi = 22e0*ko; ssioi = 52e0*ko; poritom = pori440; %для трехкомпонентной смеси
		wal = 23e0*ko; wmg = 19e0*ko; wsi = 49e0*ko; %для многокомпонентной смеси
		dkoscim = [6.07, 5.36, 6.19, 13.48, 19.93, 27.69]; %ИТОМ-440
    case (1)
		saloi = 29e0*ko; smgoi = 16e0*ko; ssioi = 55e0*ko;
		wal = 26e0*ko; wmg = 15e0*ko; wsi = 5e1*ko; poritom = pori620;
		dkoscim = [6.41, 5.53, 6.32, 13.56, 19.86, 27.66]; %ИТОМ-620
    case (2)
		saloi = 33e0*ko; smgoi = 11e0*ko; ssioi = 56e0*ko;
		wal = 3e1*ko; wmg = 1e1*ko; wsi = 52e0*ko; poritom = pori860;
		dkoscim = [7.28, 6.1, 6.83, 14.02, 19.86, 27.22]; %ИТОМ-860
    case (3)
		saloi = 35e0*ko; smgoi = 9e0*ko; ssioi = 56e0*ko; 
		wal = 33e0*ko; wmg = 7e0*ko; wsi = 53e0*ko; poritom = pori1000;
		dkoscim = [7.45, 6.24, 6.95, 14.15, 19.89, 27.09]; %ИТОМ-1000
    otherwise
        disp('Net takoy marki ITOM!');
end
        for k = 1:dkoscil
		tm = dkoscim(k)*ko; 
        dkoscim(k) = 1e0 - tm;
        end
        switch (vvm)
            case (0)
                om=dkoscim;
            case (1)
                om=dkoscit;
            case (2)
                om=[wal,wmg,wsi,saloi,smgoi,ssioi,poritom];
        end
end

function [ stch ] = oprStCheritom(vmitom,etei,vybvykh,fl,dkoscim,dkoscit,wal,wmg,wsi) 
dmkoosci=14; dtosci=1e2; tnosci=3e2;
tkusci = tnosci:dtosci:((dmkoosci-1)*dtosci+tnosci); 
kusci = koefoslab(wmg, wsi, wal, tkusci); %detei = 1e2; cemi = 11; tei0 = 273.15; tnaci = 2e2; etei = tnai:detei:((cemi-1)*detei+tnai); etei = etei + tei0; 
if (fl==1)
for k=1:length(etei)
        s = epsisreditom(etei(k), tkusci, kusci, dkoscit, dkoscim, vmitom);
        stchsritom(k) = s; 
end
end
switch (vybvykh)
    case (1)
    stch=stchsritom';
    case (2)
    stch=[wsi, wal, wmg];
    case (3)
    stch=kusci';
    case (4)
    stch=tkusci';
    case (5)
    stch=dkoscim';
    case (6)
    stch=dkoscit';
    otherwise stch=[0];
end
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

function es = epsisreditom(T, tdkusctm, dkusctm, dkoscet, dkoscem, vi)
if (vi==0)
fileID=fopen('Dliny_voln_itom.txt','r'); 
else fileID=fopen('DlinaVolny_itomN.txt','r'); 
end
formatSpec='%f';
dv=0; dv=fscanf(fileID,formatSpec); 
fclose(fileID);
npp=0; npp=Kramers_Kronig_itom(vi); 
dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, length(tdkusctm));
dkusct=ProvAdek(dkusct);
dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, length(dkoscet));
dkosce=ProvAdek(dkosce);
ume=dkusct*dkosce;
epsil=0; epsil=epsilnu(npp);
%lendv=length(dv); lenep=length(epsil); lennpp=length(npp);
epssr=usrednen(T, epsil, dv, npp);
es=ume*epssr;
end

function [ np ] = Kramers_Kronig_itom(vi)
switch (vi)
    case (0)
    np=Kramers_n_Itom(); 
    case (1)
    np=Kramers_n_Itom620(); 
    case (2)
    np=Kramers_n_Itom860(); 
    case (3)
    np=Kramers_n_Itom1000(); 
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

function [ vm ] = proverka(srk, laefm)
n=length(srk);
f=1;
e=1e-8;
for k=1:n
    r=srk(k);
    t=laefm(k);
    if (r<t)
        f=-1;
        break;
    end
    if (k>1)
        p=srk(k-1);
    if ((abs(p-r)<e) || (p>r))
        f=-1;
        break;
    end
    end
end
if (f<0)
    for k=1:n
        srk(k)=0;
    end
end
vm=srk;
end