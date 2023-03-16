function [ nf ] = DulnevZernVer(por, te, laefm,wmg,wsi,wal,vfv,la_voz,vyb)
%por = 0.6635; %пористость вермикулита фракции 2-0,7 мм %disp('По методу Дульнева для зернистых систем для вермикулита');
dkoscvm=oprdopkoefoslabstchver(1,wal,wmg,wsi);
tdkoscvm=oprdopkoefoslabstchver(2,wal,wmg,wsi);
kuscvm=oprdopkoefoslabstchver(3,wal,wmg,wsi);
tkuscvm=oprdopkoefoslabstchver(4,wal,wmg,wsi);
up=1; uo=urovPod(por);
lete=length(te);
la_e=zeros(1, lete);
stchm=zeros(1, lete);
lam=zeros(1, lete);
for k=1:length(te)
    ts=te(k);
    la_e(k)=opredKTPTKTochSha(laefm, te, ts, length(te));
dkusct=opredKTPTKTochSha(dkoscvm, tdkoscvm, ts, length(tdkoscvm));
dkusct=ProvAdek(dkusct);
dkosce=opredKTPTKTochSha(kuscvm, tkuscvm, ts, length(tkuscvm));
dkosce=ProvAdek(dkosce);
ume=dkusct*dkosce;
eps=epsisred(ts)*ume;
stchm(k)=eps;
lam(k)=opredDulnLam1Ver(por,ts,eps,la_voz(k),la_e(k),vfv);
end
lam = proverkaVer(lam, la_e, la_voz, up, uo, length(lam));
%stchm=stchm'; f=ZapisFileOptio(stchm); disp('Степень черноты'); lam=lam'; laefm=laefm'; srk=lam';
lam=provnanuli(lam)';
switch (vyb)
    case (0)
        nf=lam;
    case (1)
        nf=stchm;
    otherwise
        nf=[0];
end
end

function u = urovPod(po)
porex=95e-2; pormakvi=53e-2; pomi=36e-2;
porist=28e-2; porg=5e-1; povi=6e-1; m2=po;
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

function urot = urovOtsechen(pro,ref,urot)
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
urot=fl;
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