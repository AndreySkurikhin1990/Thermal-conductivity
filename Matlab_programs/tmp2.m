function [ w ] = tmp2()
te0=273.15; ko=1e-2; wsio=368e-3; walo=132e-3; wmgo=217e-3; kk=1e-3; %qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); 
vybvesch=1; vybmar=0; vpmf=0; vpkf=0; 
vmivmf=1; vyfv=1; vyuv=1; vysv=2; y0=3e1*1e-3;
por=novNapMas(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf)
    k=0; tholvs = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf);
	k=1; tgorvs = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf);
	k=2; qobvs = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf);
	k=3; ektpvs = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf);
	k=4; tsredvers = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf);
    cemv=length(tsredvers);
	k=0; tholv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv)
	k=k+1; tgorv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv)
	k=k+1; qobv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv)
	k=k+1; temvs = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv)
	k=k+1; laefm = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv)
    lavo=0; dlma=length(temvs); lavo=zeros(1,dlma);
for k=1:dlma
    lavo(k)=opredTeploprovVozd(temvs(k));
end
k=0; la_tk=DulnevZernVer(por, temvs, laefm, wmgo, wsio, walo, vyfv, lavo, k)
%k=k+1; s_ch=DulnevZernVer(por, temvs, laefm, wmgo, wsio, walo, vyfv, lavo, k)
%k=k+1; ume=DulnevZernVer(por, temvs, laefm, wmgo, wsio, walo, vyfv, lavo, k)
%for k=7:15
%r=provFunc2(k);
%end
%w=provFunc();
%w=poiskkb();
%vyb=1;
%w=ktpkvi(vyb);
w=[0];
end

function [ nf ] = DulnevZernVer(por, te, laefm, wmg, wsi, wal, vfv, la_voz, vyb)
%por = 0.6635; %пористость вермикулита фракции 2-0,7 мм %disp('По методу Дульнева для зернистых систем для вермикулита');
k=1; dkoscvm=oprdopkoefoslabstchver(k, wal, wmg, wsi);
k=k+1; tdkoscvm=oprdopkoefoslabstchver(k, wal, wmg, wsi);
k=k+1; kuscvm=oprdopkoefoslabstchver(k, wal, wmg, wsi);
k=k+1; tkuscvm=oprdopkoefoslabstchver(k, wal, wmg, wsi);
up=1; uo=urovPod(por);
lete=length(te);
la_e=zeros(1, lete);
stchm=zeros(1, lete);
lam=zeros(1, lete);
for k=1:lete
    ts=te(k);
    la_e(k)=opredKTPTKTochSha(laefm, te, ts, length(te));
dkusct=opredKTPTKTochSha(dkoscvm, tdkoscvm, ts, length(tdkoscvm));
dkusct=ProvAdek(dkusct);
dkosce=opredKTPTKTochSha(kuscvm, tkuscvm, ts, length(tkuscvm));
dkosce=ProvAdek(dkosce);
ume(k)=dkusct*dkosce;
eps=epsisred(ts)*ume(k);
stchm(k)=eps;
lam(k)=opredDulnLam1Ver(por, ts, eps, la_voz(k), la_e(k), vfv);
end
lam = proverkaVer(lam, la_e, la_voz, up, uo, length(lam));
%stchm=stchm'; f=ZapisFileOptio(stchm); disp('Степень черноты'); lam=lam'; laefm=laefm'; srk=lam';
lam=provnanuli(lam);
switch (vyb)
    case (0)
        nf=lam;
    case (1)
        nf=stchm;
    case (2)
        nf=ume;
    otherwise
        nf=[0];
end
end

function ol1 = opredDulnLam1Ver(po,T,eps,lavo,lae,vfv)
ladb=1e5; lada=1e-5; ra=1e2; ep=1e-7; h=0; la=0; kit=1e3;
while ((ra>ep) && (h<kit))
ladc=(lada+ladb)/2;    
fa=DulnevSloForVer(po,T,eps,lavo,lada,vfv);
fa=fa-lae;
fb=DulnevSloForVer(po,T,eps,lavo,ladb,vfv);
fb=fb-lae;
fc=DulnevSloForVer(po,T,eps,lavo,ladc,vfv);
fc=fc-lae;
la = VybCAB(fa,fb,fc,lada,ladb,ladc);
lada=la(1); ladb=la(2); ladc=la(3);
ra=abs(fa-fb);
h=h+1;
end
ol1=ladc;
end

function lam4 = DulnevSloForVer(m2,T,eps,lavo,lam1,vfv)
switch (vfv)
    case (0)
        r=(2e0+7e-1)*1e-3/2e0/2e0; 
    case (1)
        r=(8e0+4e0)*1e-3/2e0/2e0; 
end
d=2e0*r;
Nk=sqrt(m2^2-10e0*m2+9e0);
Nk=(m2+3e0+Nk)/2/m2;
%y2=3.3e-3*(1-m2)^(-2/9); %pud=0; hsl=30e-3; rona=0.2; ro1=rona/(1-m2); %y2=y2*(pud+9.8*ro1*(1-m2)*hsl)^(1/3); %y1=1e-20; y2=y1;
eta=1e-4; 
y1=10e-4; y2=y1/sqrt(eta); %y1=(10+50)*1e-4/2;
y3=2*sqrt(Nk-1)/Nk;
y4=y3/((1-m2)^(1/3));
hsh=2e-3/2;
%hsh=1e-10; %hshsr=hsh/r;
epf=1e-2;
fbn=opredFiBolN1(y2,y2/y3);
fbnn=opredFiBolNN1(lavo/lam1,m2);
if (fbn>1) 
    fb=fbnn; 
if (fbnn>1)
    fb=0;
end
elseif (fbnn>1) 
    fb=fbn; 
elseif (abs(fbn-fbnn)<epf) 
    fb=(fbn+fbnn)/2; 
elseif (fbn>fbnn) 
    fb=fbnn; 
else fb=fbn;
end
%fbs=opredFiBol(y2/y3,y2); fbn=opredFiBolN(y2,y2/y3); %if (abs(fbn-fbs)>1e-1) fb=(fbs+fbn)/2; else fb=fbs; end
A=y2^2-y1^2;
F=sqrt(1-y2^2);
D=sqrt(1-y3^2);
E=y4^2-y3^2;
delsrsz=d*(hsh+1/Nk);
lamszm=MoleSostTeplVozd(T,delsrsz,lavo);
epspr=eps/(2-eps);
sig=5.668e-8;
lamszl=4*epspr*sig*delsrsz*(T^3);
lamsz=lamszm+lamszl;
nusz=lamsz/lam1;
w=(lavo/lamszm-nusz*D)/(1-nusz);
y12=(y1^2)/(hsh/2+(1-hsh/2)*fb);
DF=abs(w-D)/abs(w-F);
nug=nusz;
numz=MoleSostTeplVozd(T,hsh*r,lavo)/lam1;
nu2sp=nug;
DF=(D-F+w*log(DF))*2*nug/(1-nug);
AF=A/(1-hsh/2-F+hsh/2/numz);
ADF=1/(AF+DF);
ADF=1/(D/(y3^2)+ADF);
ADF=ADF+nu2sp*E+y12;
lam4=ADF*lam1/(y4^2);
end

function la = MoleSostTeplVozd(T,de,lamg)
gam=7/5;
H=101325;
H0=1e5;
kB=1.380658e-23;
n=H/kB/T;
d=(0.6*0.21+0.65*0.79)*1e-10;
sig=pi*(d^2)/4;
dli=1/sqrt(2)/n/sig;
Pr=opredPrVozd(T);
a=opredKoefAkkomodLandau(T);
%a=0.9;
cz=8.42e-3;
Ty=113;
la0=cz/H0/(1+Ty/T);
Kn=la0*H0/H/de;
B=4*gam/(gam+1)*(2-a)/a*Kn/Pr;
la=lamg/(1+B);
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

function t = opredFiBolN1(x, y)
y1=[0,47e-3,113e-3,195e-3,3e-1,413e-3,563e-3,725e-3,832e-3,932e-3,1e0];
y01=[0,261e-3,411e-3,526e-3,618e-3,7e-1,774e-3,837e-3,9e-1,958e-3,1e0];
n=length(y1);
fib(1)=0; 
for k=2:n 
    fib(k)=fib(k-1)+1e-1; 
end
f1=1; f2=1; q1=0; q2=0;
for k=1:n
if ((f1>0) && (y<y01(k))) 
        q1=k; f1=0; break;
end
if ((f2>0) && (y<y1(k))) 
        q2=k; f2=0; break;
end
end
if (q1==0)
    q1=2; 
end
if (q2==0)
q2=2;
end
fb01=fib(q1-1)+(fib(q1)-fib(q1-1))*(y-y01(q1-1))/(y01(q1)-y01(q1-1));
fb1=fib(q2-1)+(fib(q2)-fib(q2-1))*(y-y1(q2-1))/(y1(q2)-y1(q2-1));
fibo=fb1+(fb01-fb1)*(x-1e0)/(1e-1-1e0);
t=fibo;
end

function t = opredFiBolNN1(nu, m2)
cb=1e2; ca=-1e1; ra=1e2; ep=1e-6;
h=0; kit=100;
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
ra=abs(fa-fb); h=h+1;
end
fi1=abs(cc/(2-cc)); 
fi1=sqrt(fi1);
t=(nu^2)*(7.5-11*nu+4.5*(nu^2))*(1-fi1)+fi1;
end

function  fb = opredFiBol(bb,rr)
nb=1e2;
sn=0;
for k=0:nb
    n=2*k+1;
    I1xy=besseli(1,n*pi*bb*rr);
    I1x=besseli(1,n*pi*rr);
    K1xy=besselk(1,n*pi*bb*rr);
    K1x=besselk(1,n*pi*rr);
   sn=sn+I1xy*(I1x*K1xy-K1x*I1xy)/(n^2)/I1x; 
end
fb=1-16*sn/(pi^2);
end

function fibo = opredFiBolN(x,y)
y1=[0,0.047,0.113,0.195,0.3,0.413,0.563,0.725,0.832,0.932,1];
y01=[0,0.261,0.411,0.526,0.618,0.7,0.774,0.837,0.9,0.958,1];
fib=0:1e-1:1e0; f1=1; f2=1; q1=1; q2=1;
for k=2:length(y1)
    if ((f1>0) && (y<y01(k)))
        q1=k; f1=0;
    end
    if ((f2>0) && (y<y1(k)))
        q2=k; f2=0;
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
fibo=fb1+(fb01-fb1)*(x-1e0)/(1e-1-1e0);
end

function fibo = opredFiBolNN(nu,m2)
cb=1e2; ca=-1e1; ra=1e2; ep=1e-9; h=0; la=0; kit=1e3; m2=m2';
while ((ra>ep) && (h<kit))
cc=(ca+cb)/2;
fa=2*(ca^3)-3*(ca^2)+1-m2;
fb=2*(cb^3)-3*(cb^2)+1-m2;
fc=2*(cc^3)-3*(cc^2)+1-m2;
la = VybCAB(fa,fb,fc,ca,cb,cc);
ca=la(1); cb=la(2); cc=la(3);
ra=abs(fa-fb); h=h+1;
end
cc=cc';
fi1=sqrt(cc/(2-cc)); fi1=abs(fi1);
fibo=(nu^2)*(7.5-11*nu+4.5*(nu^2))*(1-fi1)+fi1;
end

function [ pr ] = proverkaVer(ta, laefm, lavo, up, uo, n)
lvmi=min(lavo);
for k=1:n
	ta(k)=urovOtsechen(ta(k), laefm(k), uo); 
	ta(k)=urovPodder(ta(k), lavo(k), up);
end
f=1;
for k=1:n
    if ((ta(k)<lvmi) && (f>0)) 
        f=0; break;
    end
end
if (f==0) 
    for k=1:n 
        ta(k)=0; 
    end
end
pr=ta;
end

function t = ZapisFileOptio(massi)
fid = fopen('Stepen_chernoty_ver.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
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

function opr = opredKTPVozd(T)
te0=273.15;
te = arrTempAir()+te0;
lamb = koefTeploprovAir();
n = length(te);
lam=opredKTPTKTochSha(lamb,te,T,n);
opr = lam;
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

function epsnu = epsilnu(n)
eps=(4*n+2)/3/((n+1)^2);
eps=eps+2*(n^3)*(n^2+2*n-1)/(n^2+1)/(n^4-1);
eps=eps-8*(n^4)*(n^4+1)*log(n)/(n^2+1)/((n^4-1)^2);
eps=eps-(n^2)*log((n-1)/(n+1))*((n^2-1)^2)/((n^2+1)^3);
epsnu=eps;
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

function [ vm ] = provnanuli(srk)
f = 1; ep=1e-6;
for k=1:length(srk)
    if ((srk(k)<ep) && (f > 0))
        f = -1;
        break;
    end
end
if (f < 0)
for j=1:length(te)
    srk(j)=0;
end
end
vm=srk;
end

function [ duln ] = provFunc2(vyb)
format longg;
vmivmf=0; vyfv=0; vyuv=1; vysv=0; 
vybvesch=1; vybmar=0; vpmf=0; vpkf=0; 
por=novNapMas(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf);
k=1; tm=rasPorpoRazVer(vyfv,k);
srp=tm(1); %k=2; tm=0; tm=rasPorpoRazVer(vyfv,k); marp=tm(1)
k=3; rp=rasPorpoRazVer(vyfv, k);
s=sum(rp);
k=5; sr=rasPorpoRazVer(vyfv, k);
k=6; legr=rasPorpoRazVer(vyfv, k);
k=7; prgr=rasPorpoRazVer(vyfv, k);
laefm =[0.102082636939044, 0.136546678058965, 0.174924132522818, 0.217215000330604 ...
    0.263419281482324, 0.313536975977976, 0.367568083817561, 0.425512605001079, 0.487370539528529];
n=length(laefm); te0=273.15; t0=2e2; dt=1e2; tk=(n-1)*dt+t0; tem = t0:dt:tk; tem=tem+te0;
lamax=1e6; lamin=1e-9; hko=1e2; tocras=1e-8;
for k=1:n
ts=tem(k);
lavo(k)=opredTeploprovVozd(ts-te0);
end
dul=UchetRapsredPorPoRazm(por, laefm, lavo, vyb, rp, legr, prgr, sr, lamax, lamin, hko, tocras, n); %7 - Кубы, адиаб., 8 - Кубы, изотерм., 9 - Цилиндры, адиаб., 10 - Цилиндры, изот., 11 - По Одолевскому 1, 12 - По Одолевскому 2
duln=[0];
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

%учет распределения пор
function [ URP ] = UchetRapsredPorPoRazm(po, lameff, lamvoz, vy, rp, legr, prgr, sr, lamax, lamin, hko, tocras, nT)
le = length(legr); po=po'; %sr=(legr+prgr)/2e0;
j = 1; vo = 0; pw = 0; ep=tocras;
if (po<1e0)
    w=1e0;
else w=1e2;
end
for k=1:le
    if (rp(k) > ep)
    pw(j) = rp(k) * po / w; %объемная доля поры заданного размера в полном (во всем) объеме
    vo(j) = (sr(k)^3) / pw(j); %все поры - кубы, объем пор заданного диапазона
    j = j + 1;
    end
end
if (vy==8)
pw=pw
vo=vo
sr=sr
rp=rp
end
mx=sum(vo); up=1; uo=urovPod(po); 
lcub = mx^(1/3); lcub=lcub/length(vo); %оценка максимального значения l - размера образца
lb = lcub / (po^(1/3)); %оценка размера поры
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
        lamvy(k)=lamcubadia(k); %куб, адиаб. - 30
    case 8
        lamvy(k)=lamcubizo(k); %куб., изотерм. - 31
    case 9
        lamvy(k)=lamodo1(k); %по Одолевскому (матрицы) - 32
    case 10
        lamvy(k)=lamodo2(k); %по Одолевскому (стат. смесь) - 33
    case 11
        lamvy(k)=lamcyladi(k); %цил., адиаб. - 34
    case 12
        lamvy(k)=lamcylizo(k); %цил., изотерм. - 35
    case 13
        lamvy(k)=lamcylper(k); %ось цилиндра перпендикулярна тепловому потоку - 36
    otherwise 
        lamvy(k)=0;
end
end
lamvy=proverkakvi(lamvy, lameff, lamvoz, up, uo, nT, tocras);
switch (vy)
    case 7
        %disp('30'); %lamvoz=lamvoz'
        %lamcubadia=lamcubadia'
        %lamvy=lamvy' %куб, адиаб.
    case 8
        %disp('31'); 
        %lamcubizo=lamcubizo'
        %lamvy=lamvy' %куб., изотерм.
    case 9
        %disp('32'); 
        %lamodo1=lamodo1' %по Одолевскому (матрицы)
        %lamvy=lamvy'
    case 10
        %disp('33'); lamvy=lamvy' %по Одолевскому (стат. смесь)
    case 11
        %disp('34'); lamvy=lamvy' %цил., адиаб. 
    case 12
        %disp('35'); 
        %lamcylizo=lamcylizo' 
        %lamvy=lamvy' %цил., изотерм.
    case 13
        %disp('36'); lamvy=lamvy' %ось цилиндра перпендикулярна тепловому потоку
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
%матричная гетерогенная система - одна фаза образует связную матрицу при любой объемной концентрации этой фазы, 
%система имеет включения в виде кубов, центры которых образуют простую кубическую решетку, ребра параллельны
function du31 = DulKoefTep1OdolevMatr(lam2, m2, lame, lamb, lama, nit, ep)
k=0; ra=1e4;
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
%lamc=lamc';
du31 = lamc;
end
%нет регулярной структуры, статистическая смесь, частицы распределены хаотически
function du32 = DulKoefTep1OdolevStatSm(lam2, m2, lame, lamb, lama, nit, ep)
k=0; ra=1e4; %lam2 - КТП воздуха, ищем КТП твердого скелета
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

function du6 = DulKoefTep1CylPer(lavo,por,laef,lb, lamb, lama, nit, ep) %тепловой поток перпендикулярен оси цилиндра
k=0; ra=1e4; lam1=lavo; dt=1e0; h=1; Qo=laef*h*dt; rb=sqrt(por/pi)*lb; %lam1 - КТП воздуха, lam2 - КТП ТТ 
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

function [ vm ] = proverkakvi(ta, laefm, lavo, up, uo, n, ep)
lvmi=min(lavo);
for k=1:n
	ta(k)=urovOtsechen(ta(k), laefm(k), uo); 
	ta(k)=urovPodder(ta(k), lavo(k), up);
end
f=1;
for k=1:n
    if ((ta(k)<lvmi) && (f>0)) 
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

function [ t ] = provFunc()
te0=273.15; y0=30e-3; ko=1e-2; wsio=368e-3; walo=132e-3; wmgo=217e-3; e=1e-8; kk=1e-3;
vmivmf=1; vfv=0; vyuv=1; vysv=0;
tem = 2e2:1e2:12e2; tem=tem+te0;
%qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); 
if (vfv==0) 
    por = 0.6635; %пористость вермикулита фракции 2-0,7 мм
    r = (2+0.7)*kk/2; walo=14*ko; wsio=45.8*ko; wmgo=29*ko;
elseif (vfv==1) 
    por=0.5575; %пористость вермикулита фракции 8-4 мм
    r = (8+4)*kk/2; walo=12.9*ko; wsio=46.7*ko; wmgo=28.2*ko;
end
    tholvs = napMasEKTPVer(vfv, vysv, vmivmf, tem, 0, y0, vyuv);
	tgorvs = napMasEKTPVer(vfv, vysv, vmivmf, tem, 1, y0, vyuv);
	qobvs = napMasEKTPVer(vfv, vysv, vmivmf, tem, 2, y0, vyuv);
	ektpvs = napMasEKTPVer(vfv, vysv, vmivmf, tem, 3, y0, vyuv);
	tsredvers = napMasEKTPVer(vfv, vysv, vmivmf, tem, 4, y0, vyuv);
    t=[0];
end

function t = poiskkb()
format longg;
t07=tem_07r(); ta=tem_air(); tc=tem_cen(); ht=24;
n=mean([length(t07),length(ta),length(tc)]);
for k=1:n
    lntzmt07=log(ta-t07);
    lntzmtc=log(ta-tc);
    tx(k)=k*ht;
end
m=max(lntzmt07); p=1;
for k=1:n
if (lntzmt07(k)==m)
    p=k;
    break;
end
end
tx1=tx(p:n);
ta1=ta(p:n);
lntzmt071=lntzmt07(p:n);
lntzmtc1=lntzmtc(p:n);
n=length(tx1)
for k=1:n
    tx2(k)=tx1(k)^2;
    xy1(k)=tx1(k)*lntzmt071(k);
    xy2(k)=tx1(k)*lntzmtc1(k);
end
sxy1=sum(xy1);
sxy2=sum(xy2);
stx1=sum(tx1);
stx2=sum(tx2);
sy1=sum(lntzmt071);
sy2=sum(lntzmtc1);
k1=(n*sxy1-stx1*sy1)/(n*stx2-stx1^2)
k2=(n*sxy2-stx1*sy2)/(n*stx2-stx1^2)
b1=(sy1-k1*stx1)/n
b2=(sy2-k2*stx1)/n
%t=postrGraph(tx1,lntzmt071,lntzmtc1);
t=0;
end

function p = postrGraph(x,y1,y2)
pl=plot(x,y1,'-b',x,y2,'-k');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Время, с'}); 
ylabel({'ln(t_ж-t_в)'}); 
title({'График для определения m'});
p=0;
end

function [ w ] = ktpkvi(vyb)
kktp(3)=0;
switch (vyb)
    case 0
        kktp(2) = 0.00015; kktp(1) = 0.068; %350
    case 1
        kktp(2) = 0.000125; kktp(1) = 0.082; %400
    case 2
        kktp(2) = 0.0001;   kktp(1) = 0.103; %500
    case 3
        kktp(2) = 0.00015;  kktp(1) = 0.116; %600
    case 4
        kktp(2) = 0.00017;  kktp(1) = 0.146; %700
    case 5
        kktp(2) = 0.00018;  kktp(1) = 0.156; %800
    case 6
        kktp(2) = 0.00019;  kktp(1) = 0.185; %900
    case 7
        kktp(2) = 0.00025;  kktp(1) = 0.246; %1000
end
n=length(kktp); lam=0;
dtem=1e2; temk=1e3; temn=0;
tem=temn:dtem:temk;
for j=1:length(tem)
    lam(j)=0;
for k=1:n
    lam(j)=lam(j)+kktp(k)*(tem(j)^(k-1));
end
%disp(tem(j)); disp(lam(j));
end
w=lam;
end

function w = stroki()
str10='Pryamougolnye_Parallelepipedy_';
str2={'Tr_','dpp_','Trs_','dnra_','razdT_','nKBr_','dl_vo_'};
str3='itom-620';
str5='.txt';
str4=str2(1);
str1=strcat(str10,str4);
str1=strcat(str1,str3);
str1=strcat(str1,str5)
str1=str10;
w=0;
end
        
function m = mfunct()
te0=273.15;
pkp=[1,0.9779,0.9412,0.8414,0.7524,0.7375,0.7299,0.7247,0.7215];
tpkp=(2e2:1e2:1e3)+te0;
te=(5e1:1e2:13e2)+te0;
n=length(te);
for k=1:n
    zte(k)=opredKTPTKTochVer(pkp, tpkp, te(k));
end
zte=zte'
m=0;
end

function t = opredKTPTKTochVer(ktptks, te, temp)
n=length(te); f=1; p=1;
for k=2:n
if ((te(k)>=temp) && (f>0)) 
        p=k; f=0;
end
end
    if ((f==1) && (p==1)) 
        p=n; f=0;
    end
    if (f==0)
	ko=(ktptks(p)-ktptks(p-1))/(te(p)-te(p-1)); 
    ktp=ktptks(p-1)+ko*(temp-te(p-1));
    if (temp>te(n))
        ktp=ktptks(p)+ko*(temp-te(p));
    end
    else
        ktp=0;
    end
    if (ktp<0)
        ktp=0;
    end
t=ktp;
end

function w = tmp22()
knuSha=0; knuSha=SredGrafSha;
Sp1=dlvoVer53101;
Sp2=dlvoSham1();
leSp1=length(Sp1);
leSp2=length(Sp2);
knuSham=0;
knuSham=preobMas();
format long g;
for k=1:leSp2
    if (Sp2(k)==Sp1(1)) 
        break; end;
end
t=k; knuSham1=0;
for k=1:leSp1
    if (rem(k,2)==0)
        knuSham1(k)=(knuSha(t)+knuSha(t+1))/2;
        t=t+1;
    else
        knuSham1(k)=knuSha(t);
    end
end
for k=1:leSp1
    r=abs(knuSham1(k)-knuSham(k));
    if (r>1e-3)
        %disp(knuSham1(k));
        disp(k);
    end
end
for k=1:leSp1
    Sp1(k)=1e4/Sp1(k);
end
%for k=1:po
    %disp(knuSha(k));
%end
%for k=1:po
   %disp(Sp1(k));
%end
%disp(length(knuSham));
%disp(length(knuSham1));
%disp(length(Sp1));
b=plot(Sp1,knuSham,'-b',Sp1,knuSham1,'-g');
%b=plot(Sp1,knuSham1,'-b');
set(b,'LineWidth',1);
xlabel('Длина волны, мкм');
ylabel('Коэффициент поглощения, мкм^-1');
title('График alpha(lambda)');

%x=0.01:0.01:10;
%p=length(x);
%E1(1)=0;
%for k=1:p
%        E1(k+1)=integroexpon(1,x(k));
%end
%E2(1)=1;
%   for k=1:p
%        E2(k+1)=(-x(k))*(E1(k)-exp(-x(k))/x(k));
    %E2(k+1)=integroexpon(2,x(k));
    %end
    %E3(1)=0.5;
    %for k=1:p
       %E3(k+1)=(-(x(k)^2))*(E2(k)/x(k)-exp(-x(k))/x(k)^2)/2;
    %end
    %disp(E3(1));
    %for k=1:p
    %   disp(x(k));
    %disp(E3(k+1));
    %end
    %E4(1)=1/3;
    %for k=1:p
        %E4(k+1)=(exp(-x(k))-x(k)*E3(k))/3;
    %end
    w=0;
end

function [ w ] = tem_07r()
w=[23.5
23.5
23.3
23.6
23.6
23.8
23.8
23.9
24
24
24.1
24.2
24.3
24.3
24.4
24.5
24.6
24.7
24.8
24.8
24.9
25
25.1
25.2
25.3
25.4
25.4
25.5
25.6
25.7
25.8
25.8
26
26
26.1
26.2
26.3
26.4
26.5
26.6
26.6
26.7
26.8
26.9
27
27
27.1
27.2
27.3
27.4
27.5
27.6
27.7
27.8
27.8
27.9
28
28.1
28.2
28.3
28.3
28.4
28.5
28.5
28.5
28.7
28.8
28.8
28.9
28.9
29
29.2
29.3
29.3
29.4
29.4
29.5
29.5
29.6
29.7
29.8
29.8
29.9
30
29.9
30
30.1
30.2
30.3
30.4
30.4
30.5
30.5
30.5
30.6
30.7
30.8
30.8
30.9
30.9
30.9
30.9
31.1
31.3
31.2
31.3
31.3
31.3
31.3
31.5
31.5
31.6
31.7
31.7
31.8
31.7
31.9
32
32.1
32.2
32.3
32.3
32.4
32.4
32.5
32.5
32.6
32.7
32.8
32.8
32.9
32.9
32.9
33
33
33.1
33.1
33.2
33.2
33.3
33.3
33.4
33.4
33.5
33.5
33.6
33.7
33.7
33.7
33.8
33.8
33.8
33.9
33.9
33.9
34
34
34
34.1
34.2
34.2
34.2
34.2
34.3
34.5
34.4
34.4
34.5
34.5
34.5
34.6
34.7
34.6
34.7
34.7
34.7
34.7
34.7
34.8
34.9
34.9
34.9
34.9
34.9
35
35
35
35.1
35.1
35.2
35.2
35.2
35.3
35.3
35.3
35.3
35.4
35.4
35.4
35.5
35.5
35.5
35.6
35.7
35.6
35.6
35.7
35.7
35.7
35.8
35.7
35.8
35.8
35.8
35.8
35.9
36
36
36
36
36
36
36.1
36.1
36.1
36.1
36.1
36.1
36.2
36.2
36.2
36.2
36.3
36.3
36.3
36.3
36.3
36.4
36.4
36.4
36.4
36.4
36.5
36.5
36.5
36.5
36.5
36.5
36.6
36.6
36.6
36.6
36.6
36.6
36.7
36.7
36.7
36.7
36.7
36.8
36.8
36.8
36.8
36.8
36.8
36.9
36.9
37
37
37
36.9
36.9
36.9
37
37
37
37.1
37
37.1
37.1
37.1
37.1
37.1
37.1
37.1
37.2
37.2
37.2
37.2
37.3
37.3
37.3
37.3
37.3
37.4
37.4
37.4
37.4
37.4
37.4
37.4
37.4
37.5
37.5
37.4
37.4
37.5
37.4
37.5
37.5
37.5
37.5
37.5
37.6
37.6
37.6
37.6
37.6
37.6
37.6
37.7
37.7
37.7
37.7
37.7
37.7
37.7
37.7
37.7
37.7
37.7
37.7
37.8
37.8
37.8
37.8
37.8
37.8
37.8];
end

function [ w ] = tem_air()
w=[23.6
24.4
25.1
25.6
26.1
26.7
26.8
27.2
27.4
27.9
28.3
28.5
29
29.2
29.5
29.7
30.1
30.1
30.3
30.6
30.7
30.9
31.1
31.3
31.7
31.7
31.8
32
32.2
32.2
32.3
32.4
32.5
32.6
32.8
32.9
33
33.1
33.1
33.3
33.5
33.5
33.5
33.6
33.6
33.9
33.9
33.9
33.9
34
34
34.1
34.1
34.2
34.4
34.5
34.5
34.5
34.6
34.6
34.6
34.7
34.8
34.7
34.7
34.8
34.8
34.8
35
35
35
35.1
35.3
35.2
35.3
35.3
35.5
35.5
35.5
35.5
35.6
35.6
35.7
35.6
35.7
35.8
35.8
35.8
35.8
35.8
35.8
35.7
35.9
35.8
35.8
36.1
36.1
36.1
35.9
36.1
36.2
36.1
36.2
36.3
36.3
36.3
36.3
36.3
36.3
36.4
36.4
36.4
36.4
36.4
36.6
36.4
36.6
36.6
36.6
36.7
36.7
36.6
36.7
36.7
36.8
36.9
36.9
36.9
36.9
36.9
36.9
36.9
36.9
36.9
36.8
36.9
36.9
36.9
37
36.9
37
37
37.3
37.3
37
37.3
37.2
37
37.3
37.3
37.3
37.3
37.4
37.5
37.5
37.5
37.5
37.5
37.5
37.5
37.7
37.5
37.5
37.7
37.5
37.7
37.7
37.5
37.7
37.5
37.7
37.7
37.7
37.7
37.7
37.7
37.5
37.7
37.7
37.5
37.7
37.7
37.7
37.8
37.9
37.8
37.9
37.8
37.8
37.8
37.8
37.9
37.8
37.9
38
38
38
37.9
38
38
38
38
37.9
38
37.9
38
38.1
38.1
38.1
38.1
38.1
38.1
38.1
38.3
38.3
38.1
38.3
38.1
38.1
38.1
38.1
38.3
38.3
38.3
38.3
38.3
38.3
38.4
38.4
38.3
38.5
38.3
38.4
38.3
38.4
38.4
38.3
38.4
38.4
38.4
38.5
38.5
38.6
38.4
38.5
38.4
38.4
38.4
38.5
38.4
38.4
38.5
38.4
38.4
38.4
38.5
38.5
38.5
38.4
38.4
38.4
38.4
38.4
38.5
38.4
38.4
38.4
38.5
38.5
38.5
38.6
38.7
38.5
38.5
38.6
38.6
38.6
38.5
38.6
38.6
38.6
38.6
38.6
38.6
38.9
38.7
38.9
38.7
38.6
38.7
38.9
38.6
38.7
38.6
38.7
38.7
38.7
38.7
38.6
39
38.9
39
38.7
39
38.9
38.7
38.7
38.7
38.7
39
38.9
39
39
39
39
39
39
38.7
38.9
38.9
38.9
39
39
39
39
39
39.1
39
39
39
39
39.1
39
38.9
39
39
39.1
39
39.1];
end

function [ w ] = tem_cen()
w=[23.2
23.3
23.3
23.4
23.5
23.6
23.7
23.7
23.8
23.9
24
24.1
24.2
24.3
24.3
24.4
24.5
24.6
24.7
24.8
24.8
24.9
25
25.1
25.2
25.3
25.4
25.5
25.5
25.6
25.7
25.8
25.8
25.9
26
26.1
26.2
26.2
26.4
26.4
26.5
26.6
26.6
26.7
26.8
26.9
27
27.1
27.2
27.3
27.4
27.3
27.6
27.7
27.8
27.8
27.9
28
28.1
28.1
28.1
28.3
28.4
28.4
28.5
28.6
28.7
28.7
29
29
29.1
29.1
29.3
29.2
29.3
29.3
29.4
29.4
29.5
29.5
29.6
29.7
29.7
29.9
29.8
29.9
30
30.1
30.1
30.2
30.3
30.3
30.2
30.3
30.5
30.6
30.5
30.7
30.7
30.7
30.8
31
30.8
31.3
31.1
31.1
31.1
31.1
31.2
31.3
31.2
31.4
31.4
31.6
31.6
31.6
31.7
31.8
31.9
32
32.1
32.1
32.2
32.2
32.3
32.4
32.4
32.5
32.5
32.5
32.6
32.6
32.6
32.6
32.7
32.7
32.7
32.8
32.8
32.9
32.9
33
33
33
33.1
33.1
33.2
33.2
33.2
33.3
33.3
33.4
33.4
33.4
33.4
33.5
33.5
33.5
33.5
33.6
33.6
33.7
33.7
33.7
33.9
33.8
33.8
33.9
34
33.9
33.9
34.1
34
34.1
34.1
34.1
34.1
34.1
34.2
34.2
34.2
34.2
34.3
34.3
34.3
34.3
34.4
34.4
34.4
34.5
34.5
34.5
34.6
34.6
34.6
34.7
34.7
34.7
35
34.8
34.8
34.8
34.9
35
35
35
35
35
35
35.1
35.1
35.1
35.1
35.1
35.2
35.2
35.3
35.3
35.3
35.3
35.3
35.3
35.4
35.4
35.4
35.4
35.4
35.5
35.5
35.5
35.5
35.5
35.6
35.6
35.6
35.6
35.6
35.6
35.7
35.7
35.7
35.7
35.7
35.8
35.8
35.8
35.8
35.8
35.9
35.9
35.9
35.9
35.9
35.9
35.9
36
36
36
36
36
36
36.1
36
36.1
36.2
36.2
36.2
36.2
36.2
36.2
36.2
36.2
36.2
36.2
36.3
36.3
36.3
36.3
36.4
36.4
36.4
36.4
36.5
36.4
36.4
36.5
36.6
36.5
36.5
36.6
36.6
36.6
36.6
36.6
36.7
36.7
36.7
36.7
36.7
36.7
36.8
36.8
36.8
36.7
36.8
36.8
36.8
36.8
36.7
36.8
36.8
36.8
36.9
36.8
36.9
36.9
36.9
36.9
36.9
37
37
37
37
37
37
37
37.1
37
37.1
37.1
37
37.1
37.1
37.1
37.1
37.1
37.1
37.1
37.1];
end

function [ vyma ] = vydelPol(temvcs, temvhs, qos, ektpvs, temvss, v, n)
tev0 = 273.15; templa = 134.0*1e1 + tev0; q = 1; 
for k = 1:n
    if ((temvcs(k)>0) && (temvhs(k) > 0) && (qos(k) > 0) && (ektpvs(k) > 0) && (temvss(k) > 0) && (temvhs(k)<templa)) 
        q=q+1; 
    end
end
qn = q; hf = 1e0; nf = 0;
        for k = 1:qn
            nf = nf + hf;
        end
	q = 1; 
for k=1:n
if ((temvcs(k)>0) && (temvhs(k)>0) && (qos(k)>0) && (ektpvs(k) > 0) && (temvss(k) > 0) && (temvhs(k) < templa)) 
                switch (v)
                case (0) 
                    vm(q) = temvcs(k);
                case (1) 
                    vm(q) = temvhs(k);
                case (2) 
                    vm(q) = qos(k);
                case (3) 
                    vm(q) = temvss(k);
                case (4) 
                    vm(q) = ektpvs(k);
                end
            q=q+1;
end
end
	vyma=vm;
end
