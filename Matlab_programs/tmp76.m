%Определение распределения температур и средней температуры путем
%разложения по базису по методу Галеркина
function t = tmp76()
format longg;
delT=1e2; Tna=1e2; te0=273.15; cemv=11; Te=Tna:delT:((cemv-1)*delT+Tna); Te=Te+te0; y0=30e-3; 
%qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); kqt=polyfit(te,qo,2); r=polyval(kqt,858.34); qo=0;
vyfv=0
vysv=1
vmivmf=2
vyuv=0
temvci=napMasEKTPVer(vyfv,vysv,vmivmf,Te,0,y0,vyuv);
temvhi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,1,y0,vyuv);
qoi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,2,y0,vyuv);
ektpvi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,3,y0,vyuv);
temvsi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,4,y0,vyuv);
%---
temvc=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,0);
temvh=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,1);
qo=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,2);
temvs=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,3)'-te0
ektpv=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,4)'
n=length(temvc);
%---
for w=1:n
    Tk=temvc(w);
    T0=temvh(w);
    qob=qo(w);
    ektp=ektpv(w);
    T3=temvs(w);
    Ts=(Tk+T0)/2;
%f=opredBetaPBetaMVer(T0,Tk,y0,ektp,T0,qob,T3,te0);
end
t=0;
end

function p = int2ma(x,y)
n=length(x);
s=0;
for k=1:n-1
    ys=(y(k)+y(k+1))/2;
    dx=(x(k+1)-x(k));
    %dx=abs(dx);
    s=s+ys*dx;
end
p=s;
end

function [ mc ] = zadMasKoefc(c)
c0=c(1);
c1=c(2);
c2=c(3);
mct(1)=c0^5;
mct(2)=(c0^4)*(5*c1+3*c2);
mct(3)=9*(c0^3)*(c1^2)+6*(c0^3)*c1*c2+2*(c0^4)*c2;
mct(4)=14*(c0^3)*c1*c2+7*(c0^2)*(c1^3)+6*(c0^3)*(c2^2);
mct(5)=4*(c0^3)*(c2^2)+27*(c0^2)*(c1^2)*c2+2*c0*(c1^4)+6*(c0^2)*c1*(c2^2);
mct(6)=24*(c0^2)*c1*(c2^2)+14*c0*(c1^3)*c2+3*(c0^2)*c1*(c2^2)+3*(c0^2)*(c2^3);
mct(7)=7*(c0^2)*(c2^3)+27*c0*(c1^2)*(c2^2)+2*(c1^4)*c2;
mct(8)=14*c0*c1*(c2^3)+7*(c1^3)*(c2^2)+6*c0*c1*(c2^3);
mct(9)=5*c0*(c2^4)+9*(c1^2)*(c2^3);
mct(10)=5*c1*(c2^4);
mct(11)=c2^5;
mc=mct;
end

function [ vy ] = VybCAB(fa,fb,fc,xa,xb,xc)
if (((fa*fc)<0) && ((fb*fc)>0))
    xb=xc; 
end
if (((fa*fc)>0) && ((fb*fc)<0))
    xa=xc; 
end
vy = [xa,xb,xc];
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

function [ t ] = koefPribSha(ktp, te)
yx2=0; yx=0; le=length(ktp);
x4=0; x3=0; x2=0; x=0; y=0; de=0; de1=0; de2=0; de3=0; 	
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
	A=[x4,x3,x2; x3,x2,x; x2,x,le];
	de=A(1,1)*(A(3,3)*A(2,2)-A(3,2)*A(2,3))-A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2));
	de1=b(1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))-A(1,2)*(b(2)*A(3,3)-b(3)*A(2,3))+A(1,3)*(b(2)*A(3,2)-b(3)*A(2,2));
	de2=A(1,1)*(b(2)*A(3,3)-b(3)*A(2,3))-b(1)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))+A(1,3)*(A(2,1)*b(3)-A(3,1)*b(2));
	de3=A(1,1)*(A(2,2)*b(3)-A(3,2)*b(2))-A(1,2)*(A(2,1)*b(3)-A(3,1)*b(2))+b(1)*(A(2,1)*A(3,2)-A(3,1)*A(2,2));
	t=[de3/de,de2/de,de1/de];
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

function [ t ] = opredBetaPBetaMVer(T0,Tk,tol,koeftep,Tz,qo,T3,te0)
kt=polyfit([0,tol/2,tol],[T0,T3,Tk],2);
dko2m=[4.68,4.69,5.656,13.17,20.2,27.81]/1e2; 
dko2t=6e2:2e2:16e2; 
for k=1:length(dko2m)
    dko2m(k)=1-dko2m(k);
end
Nt=2e1+1; 
y=0:tol/Nt:tol; 
n=length(y);
for k=1:n
    Te(k)=polyval(kt,y(k));
end
Te=Te';
pkp=[1,0.971325109,0.937848799,0.836713218,0.757732344,0.735010914,0.725526014,0.717549581,0.70984116];
tpkp=(2e2:1e2:1e3)+te0;
wmg=217e-3; 
wal=132e-3; 
wsi=368e-3;
koep=rasKoefOslabEps(Te,wmg,wsi,wal); 
sig=5.67e-8; %alf=sredkopoPatch(tol,Tna,Tko)
dv=RasshDiapDlinVoln();
knuVer=RasMasKoAbs();
dko=0.444905; %из статистического моделирования
dko3=1.3002913762; %из пористой структуры самого вермикулита
knuVer=dko*dko3*knuVer;
betam=knuVer/dko; %коэффициент ослабления
npp=Kramers_n();
for k=1:length(npp)
    npp2(k)=npp(k)^2;
end
alsr=knusreddvvsrVer(Tz,knuVer,dv,npp,pkp,tpkp);
eps=epssr(Te,dv,npp,knuVer,koep,dko2m,dko2t,pkp,tpkp); %падение температуры и экспериментальные данные
f=1; e=1e-6;
for k=1:length(eps)
    if (eps(k)<e)
        f=-1;
        t=0;
        return;
    end
end
bp=(T0/Tz)^4; bp=bp/pi; 
bm=(Tk/Tz)^4; bm=bm/pi; 
n=length(Te);
for k=1:n
    als(k)=knusreddvvsrVer(Te(k),knuVer,dv,npp,pkp,tpkp);
    betas(k)=knusreddvvsrVer(Te(k),betam,dv,npp,pkp,tpkp);
    tetataus(k)=Te(k)/Tz;
    tetataus4(k)=(tetataus(k))^4;
    npr(k)=knusreddvvsrVer(Te(k),npp,dv,npp,pkp,tpkp);
    npr2(k)=knusreddvvsrVer(Te(k),npp2,dv,npp,pkp,tpkp);
end
als=als';
betas=betas';
npr2=npr2';
if (f>0)
kctp=PoiskChisKTP(qo,tol,T0,Tk,sig,betas,eps,pkp,tpkp,Te);
if (koeftep<kctp)
    t=0;
    return;
end
alss=int2ma(Te,als)/(Tk-T0)
npr2s=int2ma(Te,npr2)/(Tk-T0)
%koal=vychKoe(tpkp,Te,pkp);
taus=alss*y;
n=length(taus);
tau0=taus(n);
crp=ProvAd(kctp*alsr/4/sig/(Tz^3))/pi
%c=opredmatrc([T0,T3,Tk]/Tz,[0,tau0/2,tau0]);
c1a=-1e1; c1b=2e1; nit=1e3; k=0; ra=1e2; e=1e-8;
while ((ra>e) && (k<nit))
c1c=(c1a+c1b)/2;
ca=opredmatrc2([T0,Tk]/Tz,[0,tau0],c1a);
fa=naytia1(T0,Tk,Tz,tau0,npr2s,T3,ca,bp,bm,crp,0);
cb=opredmatrc2([T0,Tk]/Tz,[0,tau0],c1b);
fb=naytia1(T0,Tk,Tz,tau0,npr2s,T3,cb,bp,bm,crp,0);
cc=opredmatrc2([T0,Tk]/Tz,[0,tau0],c1c);
fc=naytia1(T0,Tk,Tz,tau0,npr2s,T3,cc,bp,bm,crp,0);
if ((fc*fb>0) && (fa*fc<0)) 
    c1b=c1c; 
end
    if ((fc*fa>0) && (fb*fc<0)) 
        c1a=c1c; 
    end
ra=abs(c1b-c1a);
k=k+1;
end
a1=naytia1(T0,Tk,Tz,tau0,npr2s,T3,cc,bp,bm,crp,1);
cc=cc'
m=length(cc);
n=length(taus);
    for j=1:n
        s=0; 
        d=taus(j);
    for k=1:m
    s=s+(d^(k-1))*cc(k);
    end
    teta(j)=s;
    end
    for j=1:m
        s=0; 
        d=tau0*(j-1)/2;
    for k=1:m
    s=s+(d^(k-1))*cc(k);
    end
    nb(j)=s;
    end
    T3=T3';
    T3r=nb(2)*a1*Tz;
    nb=[T0,T3,Tk]/Tz;
    n=length(teta);
    teta=teta*Tz*a1;
    teta0=teta(1);
    tetako=teta(n);
    dTp=abs(tetako-teta0);
    dTn=T0-Tk;
    Tsredk=int2ma(taus,teta)/tau0
    Tsredn=int2ma(taus,Te)/tau0
    oot=abs(Tsredk-Tsredn)*1e2/Tsredn
t=0;
end
end

function f = naytia1(T0,Tk,Tz,tau0,npr2s,T3,c,bp,bm,crp,v)        
sg=zadMasKoefc(c);
n=length(sg);
pc1=0;
for k=1:n
    pc10=sg(k)*(tau0^k)/k;
    pc1=pc1+pc10;
end
    pc1=npr2s*pc1/pi;
    pc2=opredpc2(tau0,bp,c);
    pc3=opredpc3(tau0,bm,c);
    pc4=RasFunKvadMark(tau0,c,npr2s);
    c0=c(1);
    c1=c(2);
    c2=c(3);   
    pc5=-2*crp*c2*tau0*(c0+c1*tau0/2+c2*(tau0^2)/3);
    a1=opreda1(pc1+pc4,pc2+pc3,pc5);
    m=length(c);
        s=0;
        taus=tau0/2;
        for j=1:m
            s=s+c(j)*a1*(taus^(j-1));
        end
        if (v==0)
        f=s-T3/Tz;
        elseif (v==1)
            f=a1;
        end
end

function t = opreda1(a14,a23,a5)
a1a=0; a1b=1e12; 
q=0; nq=1e5; eps=1e-4; ra=1e1;
while ((q<nq) && (ra>eps))
    a1c=(a1a+a1b)/2;
    fa=a14*(a1a^4)+a5*a1a+a23;
    fb=a14*(a1b^4)+a5*a1b+a23;
    fc=a14*(a1c^4)+a5*a1c+a23;
    la = VybCAB(fa,fb,fc,a1a,a1b,a1c);
    a1a=la(1); a1b=la(2); 
    ra=abs(a1b-a1a);
    q=q+1;
end
%if (a1c>1e2)
    %a1c=0;
%end
%if (a1c<0)
    %a1c=0;
%end
t=a1c;
end

function r = opredpc2(tau0,bp,c)
    c0=c(1);
    c1=c(2);
    c2=c(3);
    pc21=ProvAd(integroexpon(3,tau0));
    pc22=1/2;
    pc2=-c0*(pc21-pc22);
    pc21=2*ProvAd(integroexpon(4,tau0));
    pc21=pc21-exp(-tau0);
    pc22=2/3-1;
    pc2=pc2+c1*(pc21-pc22);
    pc21=-6*ProvAd(integroexpon(5,tau0));
    pc21=pc21+2*exp(-tau0);
    pc21=pc21-(tau0+1)*exp(-tau0);
    pc22=-6/4+2-1;
    pc2=pc2+c2*(pc21-pc22);
    r=-bp*pc2/2;
end

function r = opredpc3(tau0,bm,c)
    c0=c(1);
    c1=c(2);
    c2=c(3);
    pc31=ProvAd(integroexpon(3,tau0));
    pc32=1/2;
    pc3=c0*(pc31-pc32);
    pc32=-2*ProvAd(integroexpon(4,tau0));
    pc32=pc32+exp(-tau0);
    pc31=2/3-1;
    pc3=pc3+c1*(pc31+pc32);
    pc31=6*ProvAd(integroexpon(5,tau0));
    pc31=pc31-exp(-tau0);
    pc32=-1/2-tau0;
    pc3=pc3+c2*(pc31+pc32);
    r=bm*pc3/2;
end

function [ mc ] = opredmatrc(teta,x)
n=length(x);
for k=1:n
    for j=1:n
        matr(k,j)=(x(k))^(j-1);
    end
    b(k)=teta(k);
end
mc=inv(matr)*b';
end

function [ mc ] = opredmatrc2(teta,x,c1)
c0=teta(1);
tau0=x(2);
c2=(teta(2)-c0-c1*tau0);
c2=c2/(tau0^2);
mc=[c0,c1,c2];
end

function f = postgraf(koox,tetaxx,tetaxy)
pl=plot(koox,tetaxx,'-b',koox,tetaxy,'-k');
%legend('Невязка','Решение','location','best');
set(pl,'LineWidth',3);
hold on;
grid on;
%xlabel({'Координата, м'});
%ylabel({'Невязка'});
%title({'График зависимости невязки от координаты'});
f=0;
end

function mk = ProvAd(m)
if (isnan(m))   
    m=0; 
end; 
if (isinf(m))  
    m=0;  
end; 
mk=real(abs(m));
end

function t = PoiskChisKTP(qo,L,T1,T2,sig,betam,eps,pkp,tpkp,tem)
a=0; b=1e9; ra=abs(b-a);
ep=1e-4; q=0; ni=1e5;
eps1=eps(1); n=length(eps); eps2=eps(n);
Ts=(T1+T2)/2;
dko=opredKTPTKTochVer(pkp, tpkp, Ts);
betae=opredKTPTKTochVer(betam, tem, Ts)*dko;
while ((ra>ep) && (q<ni))
    c=(a+b)/2e0;
    fa=a*abs(T1-T2)/L+(sig/(3*betae*L/4+1/eps1+1/eps2-1))*abs(T1^4-T2^4)-qo;
    fb=b*abs(T1-T2)/L+(sig/(3*betae*L/4+1/eps1+1/eps2-1))*abs(T1^4-T2^4)-qo;
    fc=c*abs(T1-T2)/L+(sig/(3*betae*L/4+1/eps1+1/eps2-1))*abs(T1^4-T2^4)-qo;
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

function knus = knusreddvvsrVer(tem,knuSha,dv,npp,pkp,tpkp)
dko2=opredKTPTKTochVer(pkp, tpkp, tem);
knuSha=dko2*knuSha; %падение КП с ростом температуры
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
nc=int2ma(dl,Ibc);
nz=int2ma(dl,Ib);
knus=nc/nz;
end

function [ epsi ] = epssr(tem,dv,npp,alf,koos,dko2m,dko2t,pkp,tpkp)
n=length(npp); p=length(tem);
for j=1:p
    dko1=opredKTPTKTochVer(pkp, tpkp, tem(j));
    alf=alf*dko1;
for k=1:n
hi=dv(k)*alf(k)/4/pi/npp(k);
np=sqrt(npp(k)^2+hi^2);
ep(k)=emdiel(np);
%ep(k)=emmet(npp(k),hi);
end
dko2=opredKTPTKTochVer(dko2m, dko2t, tem(j)); %отличие от экспериментальных данных
eps(j)=epssred(dv,ep,tem(j),npp)*koos(j)*dko2;
alf=alf/dko1;
end
epsi=eps;
end

function [ m ] = zadKoef(v)
%---
x(9)=0.095012509837637;
x(10)=0.281603550779259;
x(11)=0.458016777657227;
x(12)=0.617876244402644;
x(13)=0.755404408355003;
x(14)=0.865631202387832;
x(15)=0.944575023073233;
x(16)=0.989400934991650;
%---
A(9)=0.189450610455069;
A(10)=0.182603415044924;
A(11)=0.169156519395003;
A(12)=0.149595988816577;
A(13)=0.124628971255534;
A(14)=0.095158511682493;
A(15)=0.062253523938648;
A(16)=0.027152459411754;
%---
for k=1:8
    f=9-k;
    x(f)=-x(8+k);
    A(f)=-A(8+k);
end
if (v==0)
m=x;
else 
    m=A;
end
end

function m = RasFunKvadMark(tau0,mc,n2)
tauk=zadKoef(0);
A=zadKoef(1);
a=0;
b=tau0; %delta=2;
n=length(tauk);
I=0;
c0=mc(1);
c1=mc(2);
c2=mc(3);
for k=1:n
    tau=(a+b)+(b-a)*tauk(k);
    tau=tau/2;
    An(k)=A(k); %An(k)=(b-a)*An(k)/delta;
    I=I+An(k)*funcRas(tau,tau0,c0,c1,c2);
end
m=-I*n2/2/pi;
end

function f = funcRas(tauk,tau0,c0,c1,c2)
a=tau0-tauk;
inex21=ProvAd(integroexpon(2,a));
inex22=ProvAd(integroexpon(2,tauk));
p=-c0*(inex21-inex22);
p=p+c1*(1-exp(-tau0));
inex21=ProvAd(integroexpon(3,a));
inex22=ProvAd(integroexpon(3,tauk));
p=p+c1*(inex21-inex22);
inex21=-2*ProvAd(integroexpon(4,a));
inex22=2*ProvAd(integroexpon(4,tauk));
inex23=-tau0*exp(-tau0);
p=p+c2*(inex21+inex22+inex23);
d=c0+c1*tauk+c2*(tauk^2);
d=d^4;
f=p*d;
end

function [ vm ] = opredTempHolGorVer(ektp, ete, h0, v) 
n=length(ete); nit=1e6; hf=1e0; tnoscv=3e2;
ep=1e-2; d=1e-3; Thna=tnoscv; Thnac=0; dt=1e0; 
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
qobv(k)=g;
end
g=0; 
for k=1:n
    if (qobv(k)>0)
laef=ektp(k); p=0; Thnac=Thna+g*dt; del=1e0;
while ((del>ep) && (p<nit))
tholv(k)=Thnac+p*d; 
tgorv(k)=tholv(k)+qobv(k)*h0/laef;
del=(2e0*ete(k)-(tholv(k)+tgorv(k)));
p=p+hf;
end
g=g+hf;
    else 
        tholv(k)=0; tgorv(k)=0;
    end
end
tsredver=(tgorv+tholv)/2e0; tgc=tgorv; thc=tholv; qobc=qobv; tsc=tsredver;
if (v==0) 
    vm=tgc;
elseif (v==1) 
    vm=thc;
elseif (v==2) 
    vm=qobc;
elseif (v==3) 
    vm=tsc;
end
end

function [ vm ] = vydelPol(temvc,temvh,qo,ektpv,temvs,v)
q=1; templa=17e2;
for k=1:length(qo)
    if ((temvc(k)>0) && (temvh(k)>0) && (qo(k)>0) && (ektpv(k)>0) && (temvs(k)>0) && (temvh(k)<templa))
        temvcn(q)=temvc(k);
        temvhn(q)=temvh(k);
        temvsn(q)=temvs(k);
        qon(q)=qo(k);
        ektpvn(q)=ektpv(k);
        q=q+1;
    end
end
if (v==0)
        vm=temvcn;
elseif (v==1)
        vm=temvhn;
elseif (v==2)
        vm=qon;
elseif (v==3)
        vm=temvsn;
elseif (v==4)
        vm=ektpvn;
end
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

function [ ko ] = rasKoefOslabEps(te,wmg,wsi,wal)
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
n=length(tem);
m=length(te);
for q=1:m
    f=1; j=1;
for k=1:n
    if (te(q)<tem(k))
        if (f>0)
        j=k;
        f=0;
        end
    end
end
if ((j==1) && (f==0)) 
    j=2;
end
if ((j==1) && (f==1))
    j=n;
end
ep(q)=epra(j-1)+(epra(j)-epra(j-1))*(te(q)-tem(j-1))/(tem(j)-tem(j-1));
end
koe=ep;
end