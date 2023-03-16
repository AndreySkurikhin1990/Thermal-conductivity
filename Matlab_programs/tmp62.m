%Решение задачи по методу Висканта-Гроша для шамота
function t = tmp62()
format long g;
y0=65e-3;
salosha=39e-2;
vybsha=1;
vystsha=1;
delT=1e2; Tna=1e2; Tko=12e2; te0=273.15;
Te=Tna:delT:Tko; Te=Te+te0; n=length(Te);
koef=KoefEffecKTP(salosha,vybsha,vystsha);
tepv=EffecKTPSha(koef,Te-te0);
temvh=oprTholTgor(tepv,Te,y0,1);
temvc=oprTholTgor(tepv,Te,y0,2);
ts=(temvh+temvc)/2;
tepv=oprTholTgor(tepv,Te,y0,5);
for w=1:n
Tk=temvc(w)
T0=temvh(w)
qob=tepv(w)*abs(Tk-T0)/y0;
kektp=opredKTPTKTochSha(tepv,Te,ts(w));
f=opredBetaPBetaMSha(T0,Tk,y0,kektp,T0,Te,qob);
end
t=0;
end

function [ t ] = opredBetaPBetaMSha(T0,Tk,tol,koeftep,Tz,Tet,qo)
delT=1e0; Te=Tk:delT:T0; 
wmg=0.26; wal=39; wfe=1.49; wsi=wal+wfe+wmg; wsi=1e2-wsi; 
koal=rasKoefAlpha(Te,wmg,wsi,wal); s=0; n=length(Te);
for k=1:n
    s=s+koal(k);
end
s=s/n;
sig=5.67e-8;
%alf=sredkopoPatch(tol,Tna,Tko)
dv=dlinyvoln();
fileID = fopen('Koefficient_pogloscheniya_shamota.txt','r'); formatSpec='%f'; knuSha=fscanf(fileID,formatSpec); fclose(fileID);
dko=0.647528; dko2=65e-2; dko3=0.665; s=s*dko*dko2*dko3;
knuSha=s*knuSha;
fileID = fopen('Pokazatel_prelomleniya_shamota.txt','r'); formatSpec='%f'; npp=fscanf(fileID,formatSpec); fclose(fileID);
alsr=knusreddvvsrSha(Tz,knuSha,dv,npp)
crp=ProvAd(koeftep*alsr/4/sig/(Tz^3))/pi;
kctp=PoiskChisKTP(qo,tol,T0,Tk,sig,alsr)
koeftep=koeftep'
tau0=alsr*tol; Nc=2*ceil(tau0); tau00=tau0/Nc
tol1=tau00/alsr
Nt=2e1; taus=0:tau00/Nt:tau00; ksi=taus/alsr; n=length(taus);
kna=(Tk-T0)/tau0/Tz; kb=T0/Tz;
for k=1:n
    tetataus(k)=kna*taus(k)+kb;
    tetataus4(k)=tetataus(k)^4;
end
teta004=tetataus4(1);
tetatau004=tetataus4(n);
for k=1:length(npp)
    npp2(k)=npp(k)^2;
end
%----------слой шамота тощиной tau00--------
npr=knusreddvvsrSha(Tz,npp,dv,npp);
npr2=knusreddvvsrSha(Tz,npp2,dv,npp);
for k=1:length(taus)
E01m(k)=npr2*tetataus4(k)*ProvAd(integroexpon(1,abs(taus(k))));
E11m(k)=npr2*tetataus4(k)*ProvAd(integroexpon(1,abs(tau00-taus(k))));
end
E21=pi*ProvAd(integroexpon(2,tau00));
E01i=ProvAd(trapz(taus,E01m));
E11i=ProvAd(trapz(taus,E11m));
En2t40=2*npr2*teta004;
En2t41=2*npr2*tetatau004;
%1--------
disp('First');
bpbm=reshenbpbm(pi, E21, E21, pi, En2t40-E01i,  En2t41-E11i);
bp=ProvAd(bpbm(1)); bm=ProvAd(bpbm(2));
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau00,npr2,tetataus4,taus))
%-------
%reiz=DruMetOprqr(tetataus4,taus,tau00,npr2,Tz,sig,E21,En2t41,E11i);
reiz=opredphimalqbol(tau00,sig,Tk,T0);
%reiz=resuIzlu(tol/Nc,T0,Tk,koal,dv,knuSha/s,npp,tetataus*Tz,bp*2*pi*sig*(Tz^4),bm*2*pi*sig*(Tz^4),tau00,npr,length(tetataus),Tet,ksi)
Nt=2; taus=0:tau00/Nt:tau00; n=length(taus);
kna=(Tk-T0)/tau0/Tz; kb=T0/Tz;
for k=1:n
    tetataus(k)=kna*taus(k)+kb;
    tetataus4(k)=tetataus(k)^4;
end
%bplbmi=opredIpImNSham(npp,dv,tol,knuSha,tetataus*Tz,koal,length(ksi),ksi);
t=0;
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

function [ luchteplpot ] = rasqrlamtep(bp,bm,Tz,sig,tau0,n2,teta4,taus)
n=length(taus);
for k=1:n
    E20m(k)=ProvAd(integroexpon(2,taus(k)))*n2*teta4(k);
    E21m(k)=ProvAd(integroexpon(2,tau0-taus(k)))*n2*teta4(k);
end
E3tau0=ProvAd(integroexpon(3,tau0));
E2itaus=trapz(taus,E20m)/pi;
E2tau0taus=trapz(taus,E21m)/pi;
qr0=2*pi*sig*(Tz^4)*(bp/2-bm*E3tau0-E2itaus);
qrtau0=2*pi*sig*(Tz^4)*(bp*E3tau0-bm/2+E2tau0taus);
luchteplpot=[qr0,qrtau0];
end

function [ dl ] = dlinyvoln()
Sp=(dlvoSham1()+dlvoSham2())/2;
leSp=length(Sp);
for k=1:leSp
	Sp(k)=1e-2/Sp(k);
end
dl=Sp;
end

function epsi = epssr(tem,dv,npp,alf)
n=length(npp);
np=0; hi=0; ep=0;
for k=1:n
hi=dv(k)*alf(k)/4/pi;
np=(npp(k)^2+hi^2)^0.5;
ep(k)=emdiel(np);
%ep(k)=emmet(npp(k),hi);
end
epsi=epssred(dv,ep,tem,npp);
end

function [ res ] = reshenbpbm(a11,a12,a21,a22,b1,b2)
matko=[a11, a12; a21, a22]; matsvch=[b1,b2];
res=(inv(matko)*matsvch')' ;
end

function [ rezizl ] = resuIzlu(tol,T0,Tk,koal,dl,alfs,npp,temp,Iplus,Iminus,tau0,npr,Nt,Tet,ksi)
p=length(dl); 
sigma=5.67e-8;
dko=1e0; dko2=65e-2; koal=koal*dko*dko2;
for k=1:p
    dv(k)=dl(k)/npp(k);
end
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=2*pi*PP*(c0^2);
C2=PP*c0/PB;
%1----------
taut=0; I1g=0; I2g=0;
for j=1:Nt
    ksi1=0; 
    I11g=0;
for k=1:j
    eb1=0; it1=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
altek=alfs(m)*opredKTPTKTochSha(koal,Tet,temp(k));
altej=alfs(m)*opredKTPTKTochSha(koal,Tet,temp(j));
eb1(m)=altek*c1t/(lam^5)/me;
it1(m)=eb1(m)*ProvAd(integroexpon(2,abs(altej*ksi(j)-altek*ksi(k))));
    end
    ksi1(k)=ksi(k);
    I11g(k)=integ2ma(dv,it1);
end
I1g(j)=2*pi*integ2ma(ksi1,I11g);
%2----------
    ksi2=0; I22g=0; q=1;
for k=j:Nt
    eb2=0; 
    it2=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
altek=alfs(m)*opredKTPTKTochSha(koal,Tet,temp(k));
altej=alfs(m)*opredKTPTKTochSha(koal,Tet,temp(j));
eb2(m)=altek*c1t/(lam^5)/me;
taut=abs(altek*ksi(k)-altej*ksi(j));
it2(m)=eb2(m)*ProvAd(integroexpon(2,taut));
    end
    ksi2(q)=ksi(k);
    I22g(q)=integ2ma(dv,it2);
    q=q+1;
end
I2g(j)=2*pi*integ2ma(ksi2,I22g);
end
%3----------
for j=1:Nt
for k=1:p
    altek=alfs(k)*opredKTPTKTochSha(koal,Tet,T0);
    taut=altek*ksi(j); 
    E3pdv(k)=Iplus*ProvAd(integroexpon(3,taut));
    E3mdv(k)=Iminus*ProvAd(integroexpon(3,abs(tau0-taut)));
end
Iplusi(j)=2*pi*trapz(dv,E3pdv);
Iminusi(j)=2*pi*trapz(dv,E3mdv);
end
s=0;
for j=1:Nt
qre(j)=Iplusi(j)-Iminusi(j)+I1g(j)-I2g(j);
s=s+qre(j);
end
disp(sigma*(npr^2)*abs((T0^4)-(Tk^4)));
rezizl=s/Nt;
end

function inte = integ2ma(ar1,ar2)
p=length(ar1);
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end

function [ ko ] = rasKoefAlpha(te,wmg,wsi,wal)
Tna=3e2; Tko=16e2; delT=1e2; tem=Tna:delT:Tko;
mgo = [1 0.94 0.87 0.801 0.736 0.676 0.635 0.59 ...
    0.567 0.543 0.53 0.525 0.515 0.507];
al2o3 = [1 0.98 0.946 0.898 0.854 0.797 0.753 ...
    0.709 0.676 0.642 0.618 0.598 0.582 0.57];
sio2 = [1 0.984 0.953 0.917 0.854 0.808 0.756 ...
    0.711 0.578 0.523 0.495 0.468 0.448 0.429];
wo=(wmg+wsi+wal);
kumm=vychKoe(tem,te,mgo);
kuma=vychKoe(tem,te,al2o3);
kums=vychKoe(tem,te,sio2);
n=length(te);
for k=1:n
    koo(k)=kumm(k)*wmg/wo+kuma(k)*wal/wo+kums(k)*wsi/wo;
end
ko=koo';
end

function [ koe ] = vychKoe(tem,te,epra)
for q=1:length(te)
    f=1;
for k=2:length(tem)
    if (te(q)<tem(k))
        if (f>0)
        j=k;
        f=0;
        end
    end
end
ep(q)=epra(j-1)+(epra(j)-epra(j-1))*(te(q)-tem(j-1))/(tem(j)-tem(j-1));
end
koe=ep;
end

%определяет ПТП, КТП, температуры для шамота
function [ t ] = oprTholTgor(ektp,ete,h0,vy)

n=length(ete); nit=1e6; hf=1e0; tnoscv=3e2;
ep=1e-2; d=1e-3; Thna=3e2; Thnac=0; dt=1e0; 
qon=PoiskZavVelTemVer(0);
tena=PoiskZavVelTemVer(1);
koeq=koefPribSha(qon,tena);
l=length(koeq);
for k=1:n 
    ts=ete(k); g=0; p=0; 
for j=1:l 
    g=g+(ts^p)*koeq(j); 
    p=p+hf;
end
qob(k)=g;
end
g=0; 
for k=1:n
laef=ektp(k); p=0; 
Thnac=Thna+g*dt; del=1e0;
while ((del>ep) && (p<nit))
Th(k)=Thnac+p*d; 
Tg(k)=Th(k)+qob(k)*h0/laef;
del=(2*ete(k)-(Tg(k)+Th(k)));
p=p+hf;
end
g=g+hf;
end
for k=1:n
    ktps(k)=abs(Tg(k)-Th(k))/h0;
    ktps(k)=qob(k)/ktps(k);
end
tems=(Tg+Th)/2e0; 
switch (vy)
    case 1
t=Tg';
    case 2
t=Th';
    case 3
t=qob';
    case 4
t=tems';
    case 5
t=ktps'; 
end
end

function [ klaefm ] = KoefEffecKTP(salosha,vybsha,vystsha)
if ((salosha>=28e-2) && (salosha<38e-2))
if (vybsha==0) 
        if (vystsha==0) 
            kektp(4)=-0.435e-9; 
            kektp(3)=0.685e-6; 
            kektp(2)=0.134e-3; 
            kektp(1)=0.725; 
            porsha=(20+24)*1e-2/2;
        elseif (vystsha==1)
            kektp(4)=-0.867e-9; 
            kektp(3)=1.77e-6; 
            kektp(2)=-0.523e-3; 
            kektp(1)=0.806; 
            porsha=(24+30)*1e-2/2;
        end
end %задание коэффициентов - шамот средней пористости
if (vybsha==1)
    kektp(4)=-0.397e-9; 
    kektp(3)=0.71e-6; 
    kektp(2)=0.011e-3; 
    kektp(1)=0.851; 
    porsha=(16+20)*1e-2/2;
end %уплотненный шамот
if (vybsha==2) 
    kektp(4)=-0.377e-9; 
    kektp(3)=0.918e-6; 
    kektp(2)=-0.338e-3; 
    kektp(1)=0.77; 
    porsha=(30+33)*1e-2/2; 
end %низкоплотный шамот
if (vybsha==3) 
    kektp(4)=0; 
    kektp(3)=-0.607e-6; 
    kektp(2)=1.14e-3; 
    kektp(1)=0.641; 
    porsha=(10+16)*1e-2/2;
end
end %повышенной плотности шамот
if ((salosha>=38e-2) && (salosha<=45e-2))
if (vybsha==0) 
        if (vystsha==0) 
            kektp(4)=-0.124e-9; 
            kektp(3)=0.215e-6; 
            kektp(2)=0.125e-3; 
            kektp(1)=1.01; 
            porsha=(20+24)*1e-2/2;   
        elseif (vystsha==1) 
            kektp(4)=-0.333e-9; 
            kektp(3)=0.805e-6; 
            kektp(2)=-0.289e-3; 
            kektp(1)=0.903; 
            porsha=(24+30)*1e-2/2;
        end
end %задание коэффициентов - шамот средней пористости
if (vybsha==1) 
    kektp(4)=0; 
    kektp(3)=-0.154e-6; 
    kektp(2)=0.369e-3; 
    kektp(1)=1.03; 
    porsha=(16+20)*1e-2/2;
end %уплотненный шамот
if (vybsha==2) 
    kektp(4)=-0.377e-9; 
    kektp(3)=0.918e-6; 
    kektp(2)=-0.338e-3; 
    kektp(1)=0.77; 
    porsha=(30+33)*1e-2/2;
end %низкоплотный шамот
if (vybsha==3) 
    kektp(4)=0; 
    kektp(3)=-0.141e-6; 
    kektp(2)=0.437e-3; 
    kektp(1)=1.32; 
    porsha=(10+16)*1e-2/2;
end
end %повышенной плотности шамот
klaefm=kektp';
end

function knus = knusreddvvsrSha(tem,knuSha,dv,npp)
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

function [ t ] = EffecKTPSha(koef,Te)
ktp=0;
for k=1:length(Te)
    s=0;
    for j=1:length(koef)
        s=s+koef(j)*(Te(k)^(j-1));
    end
    ktp(k)=s;
end
t=ktp';
end

function [ t ] = koefPribSha(ktp, te)
yx2=0; yx=0; x4=0; x3=0; x2=0; x=0; y=0; le=length(ktp);
for k=1:le
        yx2=yx2+ktp(k)*(te(k)^2); 
        yx=yx+ktp(k)*te(k); 
        y=y+ktp(k); 
        x4=x4+(te(k)^4); 
        x3=x3+(te(k)^3); 
        x2=x2+(te(k)^2); 
        x=x+te(k);
end
	b=[yx2,yx,y]; 
    p=le; 
    A=[x4,x3,x2;x3,x2,x;x2,x,le];
	de=inv(A)*b'; ko=de;
	ko(3)=de(1); ko(1)=de(3);
    t=ko;
end

function t = opredKTPTKTochSha(ktptks,te,temp)
n=length(te); f=1; p=0;
for k=1:n
if ((te(k)>temp) && (f>0)) 
        p=k; f=0;
end
end
    if ((p==1) && (f==0)) 
        p=2;
    end
    if ((f>0) && (p==0)) 
        p=n; f=0;
    end
    if (f==0)
	ko=(ktptks(p)-ktptks(p-1))/(te(p)-te(p-1)); 
    ktp=ktptks(p-1)+ko*(temp-te(p-1));
    end
t=ktp;
end

function oipim = opredIpImNSham(npp,dl,tol,alfs,temp,koal,Nt,ksi)
p=length(dl);
for k=1:p
    dv(k)=dl(k)/npp(k);
end
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=2*pi*PP*(c0^2);
C2=PP*c0/PB;
%1----------
I1g=0; I2g=0;
for j=1:Nt
    ksi1=0; 
    I11g=0;
for k=1:j
    eb1=0; it1=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
altek=alfs(m)*koal(k);
altej=alfs(m)*koal(j);
eb1(m)=altek*c1t/(lam^5)/me;
it1(m)=eb1(m)*ProvAd(integroexpon(2,altej*ksi(j)-altek*ksi(k)));
    end
    ksi1(k)=ksi(k);
    I11g(k)=integ2ma(dv,it1);
end
I1g(j)=2*pi*integ2ma(ksi1,I11g);
%2----------
    ksi2=0; I22g=0; q=1;
for k=j:Nt
    eb2=0; 
    it2=0;
    tem=temp(k);
    for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
lam=dv(m);
me=exp(c2t/lam/tem)-1;
altek=alfs(m)*koal(k);
altej=alfs(m)*koal(j);
eb2(m)=altek*c1t/(lam^5)/me;
taut=altek*ksi(k)-altej*ksi(j);
it2(m)=eb2(m)*ProvAd(integroexpon(2,taut));
    end
    ksi2(q)=ksi(k);
    I22g(q)=integ2ma(dv,it2);
    q=q+1;
end
I2g(j)=-2*pi*integ2ma(ksi2,I22g);
end
%4-----------
I41gm=0;
for k=1:Nt
    it4=0;
    tem=temp(k);
    for m=1:p
        alte=alfs(m)*koal(k);
        it4(m)=ProvAd(integroexpon(3,alte*ksi(k)));
    end
    I41gm(k)=2*pi*integ2ma(dv,it4);
end
%5-----------
I51gp=0; 
for k=1:Nt
    it51=0; 
    for m=1:p
        alte=alfs(m)*koal(Nt);
        tau0=alte*tol;
        it52(m)=ProvAd(integroexpon(3,tau0));
        alte=alfs(m)*koal(k);
        taut=alte*ksi(k);
        it51(m)=ProvAd(integroexpon(3,tau0-taut));
    end
    I51gp(k)=-2*pi*integ2ma(dv,it51);
end
I1g=-(I1g+I2g);
%Total-----------
Ite=oprIpImSha(I41gm(1),I51gp(1),-1,I41gm(2),I51gp(2),-1,I41gm(3),I51gp(3),-1,I1g(1),I1g(2),I1g(3));
Ite=[ProvAd(Ite(1)),ProvAd(Ite(2)),ProvAd(Ite(3))]
oipim=0;
end

function [ mas ] = oprIpImSha(a11, a12, a13, a21, a22, a23, a31, a32, a33, b1, b2, b3)
m=[a11, a12, a13; a21, a22, a23; a31, a32, a33]; b=[b1, b2, b3]; x=inv(m)*b';
mas=[ProvAd(x(1)),ProvAd(x(2)),ProvAd(x(3))];
end

function t = DruMetOprqr(tetataus4,taus,tau00,npr2,Tz,sig,E21i,En2t41,E11i)
%2-------
n=length(taus);
disp('Second');
for j=1:n
tauss=0; I1=0;
ppo(j)=0; q=1;
for k=1:j
I1(q)=npr2*tetataus4(k);
tauss(q)=taus(k);
q=q+1;
end
if (j>1)
ppo(j)=trapz(tauss,I1);
else ppo(j)=0;
end
end
vpo=trapz(taus,ppo);
%-----
E4p=ProvAd(integroexpon(4,tau00));
E4pk=-pi*(tau00/2+E4p-1/3);
E3m=ProvAd(integroexpon(3,tau00));
E4mk=-pi*(1/3-E4p-tau00*E3m);
I1=0;
for k=1:n
I1(k)=npr2*tetataus4(k);
E3ma(k)=ProvAd(integroexpon(3,abs(tau00-taus(k))))*I1(k);
E2ma(k)=ProvAd(integroexpon(2,abs(tau00-taus(k))))*I1(k);
end
I1i=trapz(taus,I1);
E3mi=trapz(taus,E3ma);
E3mi=-(vpo-E3mi+I1i/2);
E3p1=1/2-ProvAd(integroexpon(3,tau00));
E2mi=trapz(taus,E2ma);
E2mi=-(I1i+E2mi);
E2mi=E2mi/E3p1;
bpbm=reshenbpbm(pi, -pi, E4pk, E4mk, E2mi, E3mi);
bp=ProvAd(bpbm(1)); bm=ProvAd(bpbm(2));
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau00,npr2,tetataus4,taus))
%3----------
disp('Third');
bpbm=reshenbpbm(pi, -pi, E21i, pi, E2mi,  En2t41-E11i);
bp=ProvAd(bpbm(1)); bm=ProvAd(bpbm(2));
qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau00,npr2,tetataus4,taus))
t=0;
end

function t = opredphimalqbol(tau0,sig,T1,T0)
nc=1e2; tau=0:tau0/nc:tau0; ep=1e-6;
n=length(tau); etek=1e2;
for k=1:n
phimal(k)=(2/3/tau0+1-tau(k)/tau0)/(1+4/3/tau0);
end
phimalt=phimal;
nit=1e2;
for j=1:nit
for k=1:n
    for q=1:n
    phiE1(q)=phimal(q)*ProvAd(integroexpon(1,abs(tau(k)-tau(q))));
    end
    phiE1i=trapz(tau,phiE1);
    E2t=ProvAd(integroexpon(2,tau(k)));
    phimalt(k)=(E2t+phiE1i)/2;
    taus=0; Q2=0;
    for q=1:k
        Q2(q)=phimalt(q)*ProvAd(integroexpon(2,tau(k)-tau(q)));
        taus(q)=tau(q);
    end
    if (k>1)
    Q2i=trapz(taus,Q2);
    else Q2i=0;
    end
    qq=1; taus=0; Q3=0;
    for q=k:n
        Q3(qq)=phimalt(q)*ProvAd(integroexpon(2,tau(q)-tau(k)));
        taus(qq)=tau(q);
        qq=qq+1;
    end
    if (k<n)
    Q3i=trapz(taus,Q3);
    else Q3i=0;
    end
    qbol(k)=2*(ProvAd(integroexpon(3,tau(k)))+Q2i-Q3i);
end
e=vychnevyaz(phimalt,phimal);
if (e<etek)
    etek=e;
    etek=etek'
end
phimalt=phimal;
if (e<ep)
    break;
end
end
s=0; st4=0; Q=0; p=0;
for k=1:n
Q(k)=sig*abs(T1^4-T0^4)*qbol(k);
st4(k)=(abs(T1^4-T0^4)*phimal(k)+(T1^4))^(1/4);
s=s+Q(k);
p=p+st4(k);
end
disp('q R sred');
s=s/n
disp('T sred');
p=p/n
t=0;
end

function t = vychnevyaz(a1,a2)
n=length(a1); s=0;
for k=1:n
    s=s+abs(a1(k)-a2(k));
end
t=s;
end

function t = PoiskChisKTP(qo,L,T1,T2,sig,beta)
a=0; b=1e9; ra=abs(b-a); ep=1e-4; q=0; ni=1e4;
while ((ra>ep) && (q<ni))
    c=(a+b)/2e0;
    fa=a*abs(T1-T2)/L+(4*sig/3/beta/L)*abs(T1^4-T2^4)-qo;
    fb=b*abs(T1-T2)/L+(4*sig/3/beta/L)*abs(T1^4-T2^4)-qo;
    fc=c*abs(T1-T2)/L+(4*sig/3/beta/L)*abs(T1^4-T2^4)-qo;
if ((fc*fb>0) && (fa*fc<0)) 
    b=c; 
end
    if ((fc*fa>0) && (fb*fc<0)) 
        a=c; 
    end
    q=q+1; ra=abs(b-a);
end
t=c;
end

function [ mk ] = napMasEKTPVer(vyfv,vysv,vmivmf,te,v,h,vyuv)
te0=273.15; F=13.85e-4; tc84=arrTemCold()+te0; tc207=arrTemCold207()+te0;
th84=arrTemHigh()+te0; th207=arrTemHigh207()+te0;
ts207=(th207+tc207)/2; ts84=(th84+tc84)/2;
switch (vyfv) %фракция 2-0,7 мм
    case (0)
        temvs=ts207;
        temvh=th207;
        temvc=tc207;
        tepv=arrTepPot207();
        tepv=tepv/F;
        for k=1:length(ts207)
            ktp(k)=abs(temvh(k)-temvc(k))/h;
            ktp(k)=tepv(k)/ktp(k);
        end
        switch (vysv)
            case (0) %исходный
                koeft=oprkoefKTPiskhchao(vmivmf,0,temvs,tepv,h,temvh,temvc,ktp);
                koefh=oprkoefKTPiskhchao(vmivmf,1,temvs,tepv,h,temvh,temvc,ktp);
                koefc=oprkoefKTPiskhchao(vmivmf,2,temvs,tepv,h,temvh,temvc,ktp);
                koefs=oprkoefKTPiskhchao(vmivmf,3,temvs,tepv,h,temvh,temvc,ktp);
                koefq=oprkoefKTPiskhchao(vmivmf,4,temvs,tepv,h,temvh,temvc,ktp);
            case (1) %после повторных измерений
                koefq=danPoTemTepl2072(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
                koeft=danPoTemTepl2072(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
                koefh=danPoTemH2072(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
                koefc=danPoTemC2072(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
                koefs=danPoTemC2072(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t; 
            case (2) %после обжига при 1000 °С
                koefq=danPoTemTepl2073(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
                koefh=danPoTemH2073(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
                koefc=danPoTemC2073(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
                koeft=danPoTemTepl2073(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
                koefs=danPoTemTepl2073(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t; 
            case (3) %после повторного обжига при 1000 °С
                koefq=danPoTemTepl2074(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
        end
    case (1) %фракция 8-4 мм
        temvs=ts84;
        temvh=th84;
        temvc=tc84;
        tepv=arrTepPot84();
        tepv=tepv/F;
        for k=1:length(ts84)
            ktp(k)=abs(th84(k)-tc84(k))/h;
            ktp(k)=tepv(k)/ktp(k);
        end
        switch (vyuv) 
            case (1) %плоско-параллельная засыпка
                switch (vysv)
                    case (0) %исходный
                        koefh=danPoTemH(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
                        koefc=danPoTemC(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
                        koeft=danPoTemTepl(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
                        koefq=danPoTemTepl(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
                        koefs=danPoTemTepl(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t; 
                    case (1) %после повторных измерений
                        koefh=danPoTemH3(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
                        koefc=danPoTemC3(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
                        koeft=danPoTemTepl3(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
                        koefq=danPoTemTepl3(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
                        koefs=danPoTemTepl3(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t; 
                    case (2) %после обжига
                        koefh=danPoTemH6(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
                        koefc=danPoTemC6(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
                        koeft=danPoTemTepl6(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
                        koefq=danPoTemTepl6(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
                        koefs=danPoTemTepl6(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t; 
                end
            case (2) %вертикальная засыпка
                switch (vysv)
                    case (0) %исходный
                        koefh=danPoTemH2(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
                        koefc=danPoTemC2(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
                        koefq=danPoTemTepl2(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
                        koeft=danPoTemTepl2(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
                        koefs=danPoTemTepl2(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t; 
                    case (1) %после повторных измерений
                        koefh=danPoTemH5(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
                        koefc=danPoTemC5(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
                        koeft=danPoTemTepl5(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
                        koefq=danPoTemTepl5(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
                        koefs=danPoTemTepl5(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t; 
                    case (2) %после обжига
                        koefh=danPoTemH4(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
                        koefc=danPoTemC4(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
                        koeft=danPoTemTepl4(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
                        koefq=danPoTemTepl4(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
                        koefs=danPoTemTepl4(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t; 
                end
        end
end
switch (v)
    case (0)
        vm=koefc;
    case (1)
        vm=koefh;
    case (2)
        vm=koefq;
    case (3)
        vm=koeft;
    case (4)
        vm=koefs;
end
for k=1:length(te)
    t=te(k); s=0;
    for j=1:length(vm)
        s=s+vm(j)*(t^(j-1));
    end
    ma(k)=s;
end
mk=ma;
end

function [ mk ] = oprkoefKTPiskhchao(vmiv,v,temvs,tepv,h,temvh,temvc,ktp)
te0=273.15;
switch (vmiv) %выбор метода измерений
    case (0) %установка Netzsch
        mt=arrTem_VVF2(); mt=mt+te0;
        ktpn=arrKTP_VVF2();
        tsredver=opredTempHolGorVer(ktpn, mt, h, 3);
        tgorv=opredTempHolGorVer(ktpn, mt, h, 0);
		tholv=opredTempHolGorVer(ktpn, mt, h, 1);
		qobv=opredTempHolGorVer(ktpn, mt, h, 2);
        thv=vydelPol(tholv,tgorv,qobv,ktpn,tsredver,0);
        tgv=vydelPol(tholv,tgorv,qobv,ktpn,tsredver,1);
        qov=vydelPol(tholv,tgorv,qobv,ktpn,tsredver,2);
        tsv=vydelPol(tholv,tgorv,qobv,ktpn,tsredver,3);
        kt=vydelPol(tholv,tgorv,qobv,ktpn,tsredver,4);
koeft=koefPribSha(kt, tsv); koefh=koefPribSha(tgv, tsv); koefc=koefPribSha(thv, tsv); 
koefs=koefPribSha(tsv, tsv); koefq=koefPribSha(qov, tsv); 
    case (1) %данные 2020 года
        ktp1=arrKTP_VVF1(); ktpn=[(ktp1(1)+ktp1(2)+ktp1(3))/3,ktp1(4),(ktp1(5)+ktp1(6))/2];
        t10=arrTem1VVF1(); t1=[(t10(1)+t10(2)+t10(3))/3,t10(4),(t10(5)+t10(6))/2]+te0;
        t20=arrTem2VVF1(); t2=[(t20(1)+t20(2)+t20(3))/3,t20(4),(t20(5)+t20(6))/2]+te0;
        t30=arrTem3VVF1(); t3=[(t30(1)+t30(2)+t30(3))/3,t30(4),(t30(5)+t30(6))/2]+te0;
        ts=(t1+t2)/2;
koeft=koefPribSha(ktpn,ts); 
koefh=koefPribSha(t1,ts); 
koefc=koefPribSha(t2,ts); 
koefs=koefPribSha(t3,ts); 
for k=1:length(ktpn)
    s=0;
    for j=1:length(koeft)
        s=s+koeft(j)*(ts(j)^(j-1));
    end
    qv(k)=s*abs(t1(k)-t2(k))/h;
end
koefq=koefPribSha(qv,ts); 
    case (2) %данные 2019 года
        koefq=danPoTemTepl2071(temvs,tepv); t=koefq(1); koefq(1)=koefq(2); koefq(2)=t; 
        koefh=danPoTemH2071(temvs,temvh); t=koefh(1); koefh(1)=koefh(2); koefh(2)=t; 
        koefc=danPoTemC2071(temvs,temvc); t=koefc(1); koefc(1)=koefc(2); koefc(2)=t; 
        koeft=danPoTemTepl2071(temvs,ktp); t=koeft(1); koeft(1)=koeft(2); koeft(2)=t; 
        koefs=danPoTemTepl2071(temvs,temvs); t=koefs(1); koefs(1)=koefs(2); koefs(2)=t;
end
switch (v)
    case (0)
        p=koeft;
    case (1)
        p=koefh;
    case (2)
        p=koefc;
    case (3)
        p=koefs;
    case (4)
        p=koefq;
    otherwise
        p=[0,0];
end
mk=p;
end

function [ vm ] = PoiskZavVelTemVer(v)
format longg; 
te0=273.15; delT=1e2; Tna=2e2+te0; n=8; Te(1)=Tna;
for k=2:n
    Te(k)=Te(k-1)+delT;
end
y0=30e-3; 
for k=1:n
    temvc(k)=0; temvh(k)=0; temvs(k)=0; qo(k)=0; ektpv(k)=0;
end
vyfv=0; vysv=0; vmivmf=1; vyuv=1; c=0;
nvyfv=[0,1]; nvysv=[0,1,2]; nvmivmf=[1,2]; nvyuv=[1,2];
for kvf=1:length(nvyfv)
    vyfv=nvyfv(kvf);
    for jvsv=1:length(nvysv)
        vysv=nvysv(jvsv);
        if (vyfv==0)
            for qvmi=1:length(nvmivmf)
                vmivmf=nvmivmf(qvmi);
                temvct=napMasEKTPVer(vyfv,vysv,vmivmf,Te,0,y0,vyuv);
                temvc=temvct+temvc;
                temvht=napMasEKTPVer(vyfv,vysv,vmivmf,Te,1,y0,vyuv);
                temvh=temvht+temvh;
                qot=napMasEKTPVer(vyfv,vysv,vmivmf,Te,2,y0,vyuv);
                qo=qo+qot;
                ektpvt=napMasEKTPVer(vyfv,vysv,vmivmf,Te,3,y0,vyuv);
                ektpvt=ektpvt+ektpv;
                c=c+1;
            end
        elseif (vyfv==1)
            for qvuv=1:length(nvyuv)
                vyuv=nvyuv(qvuv);
                temvct=napMasEKTPVer(vyfv,vysv,vmivmf,Te,0,y0,vyuv);
                temvc=temvct+temvc;
                temvht=napMasEKTPVer(vyfv,vysv,vmivmf,Te,1,y0,vyuv);
                temvh=temvht+temvh;
                qot=napMasEKTPVer(vyfv,vysv,vmivmf,Te,2,y0,vyuv);
                qo=qo+qot;
                ektpvt=napMasEKTPVer(vyfv,vysv,vmivmf,Te,3,y0,vyuv);
                ektpv=ektpvt+ektpv;
                c=c+1;
            end
        end
    end
end 
c=c';
temvc=temvc'/c;
temvh=temvh'/c;
qo=qo'/c;
temvs=(temvc+temvh)/2;
ektpv=ektpv'/c;
switch (v)
    case (0)
        vm=qo;
    case (1)
        vm=temvs;
end
end