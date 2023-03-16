%КВИ. Келлет. Случай 2 - непрерывные температуры
function t = Kellet_cas_2_kvi()
sig=5.668e-8; te0=273.15; y0=3e1*1e-3; d=y0;
%npp=Kramers_n_kvi400();
%npp=Kramers_n_kvi500();
%npp=Kramers_n_kvi600();
%npp=Kramers_n_kvi700();
%npp=Kramers_n_kvi800();
%npp=Kramers_n_kvi900();
npp=Kramers_n_kvi1000();
%dv=dlinyvoln();
dv=dlinyvolnkvi();
fileID=fopen('Koefficient_pogloscheniya_kvi1000.txt','r'); 
formatSpec='%f';
knukvi=fscanf(fileID,formatSpec); 
fclose(fileID);
vybkvi=10; %4 - 400, 5 - 500, 6 - 600, 7 - 700, 8 - 800, 9 - 900, 10 - 1000
wsam=SodSiOAlOMgO(vybkvi, 0);
wsi=wsam(1);
wal=wsam(2);
wmg=wsam(3);
dko2m=SodSiOAlOMgO(vybkvi, 1);
tnd=6e2; dtd=2e2; dkoscil=6; 
dkoscit=tnd:dtd:((dkoscil-1)*dtd+tnd);
delT=1e2; Tna=1e2; Tko=12e2;
Te=Tna:delT:Tko; Te=Te+te0; n=length(Te);
koef=KoefEffecKTPkvi(vybkvi); 
ektpvk=EffecKTPkvi(koef,Te);
temvhk=opredTempHolGorVer(ektpvk, Te, y0, 0);
temvck=opredTempHolGorVer(ektpvk, Te, y0, 1);
qok=opredTempHolGorVer(ektpvk, Te, y0, 2);
temvsk=opredTempHolGorVer(ektpvk, Te, y0, 3);
temvc=vydelPol(temvck,temvhk,qok,ektpvk,temvsk,0);
temvh=vydelPol(temvck,temvhk,qok,ektpvk,temvsk,1);
qo=vydelPol(temvck,temvhk,qok,ektpvk,temvsk,2);
temvs=vydelPol(temvck,temvhk,qok,ektpvk,temvsk,3);
ektpk=vydelPol(temvck,temvhk,qok,ektpvk,temvsk,4);
for vtem=1:length(qo) 
ektp=ektpk(vtem)
Tz=temvh(vtem); T0=Tz;
Tk=temvc(vtem);
Tsr=temvs(vtem);
tem=Tk:1:Tz;
%dko=0.644057; dko3=0.778676; 
for k=1:length(dv)
    hikvi(k)=dv(k)*knukvi(k)/4/pi;
end
alsr=knusreddvvsrVer(Tz,knukvi,dv,npp)
hisr=knusreddvvsrVer(Tz,hikvi,dv,npp)
ronukvi=opredronu(npp, hikvi);
rosr=knusreddvvsrVer(Tz,ronukvi,dv,npp)
eps=epssr([T0,Tk],dv,npp,knukvi,vybkvi);
lamb=PoiskChisKTP(qo(vtem),d,T0,Tk,sig,alsr,eps)

lambm(vtem)=lamb;
alsrm(vtem)=alsr;
hisrm(vtem)=hisr;
rosrm(vtem)=rosr;

tC1=temvh(vtem); teta1=tC1
tC2=temvc(vtem); teta2=tC2
tetas=(teta1+teta2)/2;

m=(2*alsr*sig*(teta1^4-teta2^4)/(lamb*(teta1-teta2))+alsr^2)^(1/2); m=abs(m)
H1=2*sig*((teta1^4)-(teta2^4));
H1=H1+alsr*lamb*(teta1-teta2);
H2=m*lamb*(teta1-teta2)*(1+rosr)*(1+exp(-m*d));
H2=H2/((1-rosr)*(1-exp(-m*d)));
H1=H1+H2;
H2=2*(1+rosr)/(1-rosr);
H2=H2+alsr*d;
H3=2*(1+rosr)*(alsr^2);
H3=H3/((m^2)*(1-rosr));
H2=H2-H3;
H3=d*(1+rosr)*(1+exp(-m*d))*(alsr^2);
H3=H3/(m*(1-rosr)*(1-exp(-m*d)));
H2=H2+H3;
H=H1/H2

ptp(vtem)=H;
mm(vtem)=m;

rop=(1+rosr)/(1-rosr);
B=(1/2)*(teta1-teta2-(H*d*(alsr^2))/(lamb*(m^2)))/(1-exp(-m*d))
A=B*(1+exp(-m*d))*lamb*m+(alsr^2)*H/(m^2)-H;
A1=rop*A+alsr*lamb*teta1+2*sig*(teta1^4);
A=2*sig*(teta2^4);
A=A+alsr*lamb*teta2+alsr*H*d;
A=A-m*lamb*B*(1+exp(-m*d))*rop;
A2=A-rop*((alsr/m)^2)*H;
A=(A1+A2)/2
Am(vtem)=A;
Nx=1e2;
tetax=0;
koo=0:d/Nx:d;
C=-B*exp(-m*d);
for k=1:length(koo)
tetax(k)=B*exp(-m*koo(k))+C*exp(m*koo(k))-(koo(k)*H*(alsr^2)-alsr*A-6*alsr*sig*(tetas^4))/(lamb*(m^2));
end
%knukvi=knukvi/s;
end
H=ZapisFileTg(temvh);
H=ZapisFileTkh(temvc);
H=ZapisFileptp(ptp);
H=ZapisFilelamb(lambm);
H=ZapisFilealsr(alsrm);
H=ZapisFilerosr(rosrm);
H=ZapisFilehisr(hisrm);
H=ZapisFilemm(mm);
H=ZapisFileAm(Am);
H=ZapisFileektp(ektpk);
tetax=tetax-te0;
p=plot(koo*1e3,tetax,'-b');
set(p,'LineWidth',3);
hold on;
grid on;
xlabel({'Координата, мм'});
ylabel({'Температура, °С'});
title({'График зависимости T(x) по модели Келлета (случай 2)'});
t=0;
end

function t = opredKTPTKTochVer(ktptks, te, temp)
n=length(te); f=1; p=0;
for k=1:n
if ((te(k)>=temp) && (f>0)) 
        p=k; f=0;
end
end
    if ((p==1) && (f==0)) 
        p=2; 
    end
    if ((f==1) && (p==1)) 
        p=n; f=0;
    end
    if (f==0)
	ko=(ktptks(p)-ktptks(p-1))/(te(p)-te(p-1)); 
    ktp=ktptks(p-1)+ko*(temp-te(p-1));
    else ktp=0;
    end
    if (ktp<0)
        ktp=0;
    end
t=ktp;
end

function t = opreddko2(T0,Tk)
dko2m=[5,5,6,13,20,28]/1e2; dko2t=6e2:2e2:16e2; 
for k=1:length(dko2m)
    dko2m(k)=1-dko2m(k);
end
t=opredKTPTKTochVer(dko2m, dko2t, (T0+Tk)/2);
end

function knus = knusreddvvsrVer(tem,knuSha,dv,npp)
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

function [ t ] = opredronu(npp, hi)
        p=length(npp);
for k=1:p      
    nsha(k)=(npp(k)^2+hi(k)^2); 
    nsha(k)=abs(nsha(k));
    ronu(k)=(1/2)+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2);
    ronu(k)=ronu(k)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); 
    podln=(nsha(k)-1)/(nsha(k)+1);
    podln=abs(podln);
    ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));
    ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log(podln)/((nsha(k)^2+1)^3);
end
t=ronu;
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
        mt=arrTem_VVF2(); mt=mt'+te0;
        ktpn=arrKTP_VVF2();
        tsredver=opredTempHolGorVer(ktpn, mt, h, 3);
        tgorv=opredTempHolGorVer(ktpn, mt, h, 0);
		tholv=opredTempHolGorVer(ktpn, mt, h, 1);
		qobv=opredTempHolGorVer(ktpn, mt, h, 2);
        q=1;
        for k=1:length(qobv)
        if ((tsredver(k)>0) && (tgorv(k)>0) && (tholv(k)>0) && (qobv(k)>0) && (ktpn(k)>0))
            tsv(q)=tsredver(k); tgv(q)=tgorv(k); thv(q)=tholv(k); qov(q)=qobv(k); kt(q)=ktpn(k); q=q+1;
        end
        end
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

function [ vm ] = opredTempHolGorVer(ektp, ete, h0, v) 
del=1e0; nit=1e6; hf=1e0; ep=1e-2; d=1e-3; Thna=295; Thnac=0; dt=1e0; 
qon=PoiskZavVelTemVer(0); tena=PoiskZavVelTemVer(1);
koeq=koefPribSha(qon,tena); l=length(koeq); n=length(ete);
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
laef=ektp(k); p=0; Thnac=Thna+g*dt; del=1e0;
while ((del>ep) && (p<nit))
tholv(k)=Thnac+p*d; 
tgorv(k)=tholv(k)+qobv(k)*h0/laef;
del=(2e0*ete(k)-(tholv(k)+tgorv(k)));
p=p+hf; 
end
g=g+hf; 
end
tsredver=(tgorv+tholv)/2e0;  
switch (v)
    case (0)
        vm=tgorv;
    case (1)
        vm=tholv;
    case (2)
        vm=qobv;
    case (3)
        vm=tsredver;
end
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

function [ vm ] = PoiskZavVelTemVer(v)
format longg; 
delT=1e2; Tna=2e2; Tko=9e2; te0=273.15;
Te=Tna:delT:Tko; Te=Te+te0; n=length(Te);
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

function [ vm ] = vydelPol(temvc,temvh,qo,ektpv,temvs,v)
q=1;
for k=1:length(qo)
    if ((temvc(k)>0) && (temvh(k)>0) && (qo(k)>0) && (ektpv(k)>0) && (temvs(k)>0))
        temvcn(q)=temvc(k);
        temvhn(q)=temvh(k);
        temvsn(q)=temvs(k);
        qon(q)=qo(k);
        ektpvn(q)=ektpv(k);
        q=q+1;
    end
end
switch (v)
    case 0
        vm=temvcn;
    case 1
        vm=temvhn;
    case 2
        vm=qon;
    case 3
        vm=temvsn;
    case 4
        vm=ektpvn;
end
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

function t = PoiskChisKTP(qo,L,T1,T2,sig,beta,eps)
a=0; b=1e9; ra=abs(b-a); ep=1e-4; q=0; ni=1e4;
while ((ra>ep) && (q<ni))
    c=(a+b)/2e0;
    fa=a*abs(T1-T2)/L+(sig/(3*beta*L/4+1/eps(1)+1/eps(2)-1))*abs(T1^4-T2^4)-qo;
    fb=b*abs(T1-T2)/L+(sig/(3*beta*L/4+1/eps(1)+1/eps(2)-1))*abs(T1^4-T2^4)-qo;
    fc=c*abs(T1-T2)/L+(sig/(3*beta*L/4+1/eps(1)+1/eps(2)-1))*abs(T1^4-T2^4)-qo;
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

function [ epsi ] = epssr(tem,dv,npp,alf,vybkvi)
n=length(npp); p=length(tem);
wsam=SodSiOAlOMgO(vybkvi, 0);
wsi=wsam(1);
wal=wsam(2);
wmg=wsam(3);
dko2m=SodSiOAlOMgO(vybkvi, 1);
tnd=6e2; dtd=2e2; dkoscil=6; dkoscit=tnd:dtd:((dkoscil-1)*dtd+tnd);
koal=rasKoefAlpha(tem,wmg,wsi,wal); 
for j=1:p
    s=koal(j);
    dko1=0.647528; 
    dko3=0.6238223;
    dko2=opredKTPTKTochVer(dko2m, dkoscit, tem(j)); %поправка на экспериментальные значения
    s=s*dko1*dko2*dko3;
    alf=s*alf;
for k=1:n
hi=dv(k)*alf(k)/(4*pi*npp(k));
np=(npp(k)^2+hi^2)^0.5;
ep(k)=emdiel(np); %ep(k)=emmet(npp(k),hi);
end
eps(j)=epssred(dv,ep,tem(j),npp);
alf=alf/s;
end
epsi=eps;
end

function t = ZapisFileTg(massi)
fid = fopen('Kellet_kvi_sluchay_2_Temperatura_goryachaya.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileTkh(massi)
fid = fopen('Kellet_kvi_sluchay_2_Temperatura_kholodnaya.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileptp(massi)
fid = fopen('Kellet_kvi_sluchay_2_PTP.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFilelamb(massi)
fid = fopen('Kellet_kvi_sluchay_2_CHKTP.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFilealsr(massi)
fid = fopen('Kellet_kvi_sluchay_2_Alpha.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFilehisr(massi)
fid = fopen('Kellet_kvi_sluchay_2_Pokazatel_Pogloscheniya.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFilerosr(massi)
fid = fopen('Kellet_kvi_sluchay_2_Koefficient_Otrazheniya.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFilemm(massi)
fid = fopen('Kellet_kvi_sluchay_2_mm.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileAm(massi)
fid = fopen('Kellet_kvi_sluchay_2_Am.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileektp(massi)
fid = fopen('Kellet_kvi_sluchay_2_EKTP.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function [ dl ] = dlinyvolnkvi()
fileID=fopen('Dliny_voln_kvi.txt','r'); 
formatSpec='%f';
dv=fscanf(fileID,formatSpec); 
fclose(fileID);
le=length(dv);
for k=1:le
    dv(k)=1e-2/dv(k);
end
dl=dv;
end

function [ wsam ] = SodSiOAlOMgO(vybkvi, vyvz)
ko=1e-2;
porkvi400=52*ko;
porkvi500=52*ko;
porkvi600=52*ko;
porkvi700=40*ko;
porkvi800=41*ko;
porkvi900=39*ko;
porkvi1000=36*ko;
if (vybkvi==4) 
salok = 33; smgok = 15; ssiok = 52; porkvi = porkvi400; 
wal = 25e0; wsi = 11e0; wmg = 4e1; %КВИ-400
elseif (vybkvi==5) 
salok = 34; smgok = 11; ssiok = 54;
wal = 28; wmg = 8; wsi = 44; porkvi = porkvi500; %КВИ-500
elseif (vybkvi==6) 
salok = 36; smgok = 9; ssiok = 55;
wal = 3e1; wsi = 7; wmg = 45; porkvi = porkvi600; %КВИ-600
elseif (vybkvi==7) 
salok = 37e0; smgok = 8e0; ssiok = 55e0; 
wal = 31e0; wmg = 6e0; wsi = 45e0; porkvi = porkvi700; %КВИ-700
elseif (vybkvi==8) 
salok = 38e0; smgok = 7e0; ssiok = 55e0; 
wal = 3e1; wmg = 5e0; wsi = 45e0; porkvi = porkvi800; %КВИ-800
elseif (vybkvi==9) 
salok = 39e0; smgok = 6e0; ssiok = 55e0; 
wal = 32e0; wmg = 5e0; wsi = 45e0; porkvi = porkvi900; %КВИ-900
elseif (vybkvi == 10) 
salok = 39e0; smgok = 6e0; ssiok = 55e0; 
wal = 32e0; wmg = 4e0; wsi = 45e0; porkvi = porkvi1000; %КВИ-1000
end
salok = salok*ko; smgok = smgok*ko; ssiok = ssiok*ko;
wal = wal*ko; wsi = wsi*ko; wmg = wmg*ko;
k = 1; dkosckm(k) = 3e0; k=k+1; 
dkosckm(k) = 6e0; k=k+1; 
dkosckm(k) = 7e0; k=k+1; 
dkosckm(k) = 15e0; k=k+1; 
dkosckm(k) = 2e1; k=k+1; 
dkosckm(k) = 26e0;
if (vybkvi>6) 
    k=1; dkosckm(k)=7e0;
end
if (vybkvi>9)
    k=2; dkosckm(k)=8e0;
end
dkosckl=length(dkosckm);
for k = 1:dkosckl
tm = dkosckm(k)*ko; 
dkosckm(k) = 1e0 - tm; 
end
if (vyvz==0)
wsam=[wsi,wal,wmg];
elseif (vyvz==1) 
wsam=dkosckm;
else wsam=[0,0,0];
end
end

function p = PraCha10(phi, phis, n1, n2)
p=n1*sin(phi)-n2*sin(phis);
end

function p = provUgla(ugol)
ugo=ugol;
if (ugo>(pi/2e0)) 
    ugo=pi-ugo;
end
if (ugo<0) 
    ugo=abs(ugo);
end
p=ugo;
end

function prf = PoiskReflPhi(phi, phis, n1, n2)
rpa=(n2*cos(phi)-n1*cos(phis))/(n2*cos(phi)+n1*cos(phis)); %коэффициент отражения параллельный плоскости падения
rpe=(n1*cos(phi)-n2*cos(phis))/(n1*cos(phi)+n2*cos(phis)); %коэффициент отражения перпендикулярный к плоскости падения
prf=(abs((rpa^2e0))+abs((rpe^2e0)))/2e0;
end

function [ klaefm ] = KoefEffecKTPkvi(no)
tei0=273.15; te200=2e2+tei0; te380=38e1+tei0;
if (no==4) 
    k=0.00016;
    b=0.14;
	kkvi400(2)=k; 
    kkvi400(1)=b;
	kti=kkvi400;
elseif (no==5)
    k=0.000126;
    b=0.101842;
	kkvi00(2)=k; 
    kkvi00(1)=b;
	kti=kkvi00;
elseif (no==6)
	k=0.00015;
    b=0.116;
	kkvi600(2)=k; 
    kkvi600(1)=b;
	kti=kkvi600;
elseif (no==7) 
	k=0.000139;
    b=0.146526;
	kkvi700(2)=k; 
    kkvi700(1)=b;
	kti=kkvi700;
elseif (no==8) 
    k=0.000139;
    b=0.156526;
	kkvi800(2)=k; 
    kkvi800(1)=b;
	kti=kkvi800;
elseif (no==9) 
    k=7.36842e-5;
    b=0.193157895;
    kkvi900(2)=k; 
    kkvi900(1)=b;
	kti=kkvi900;
elseif (no==10) 
    k=7.78947e-5;
    b=0.248052632;
    kkvi1000(2)=k; 
    kkvi1000(1)=b;
	kti=kkvi1000;
end
kti(3)=0;
klaefm=kti';
end

function [ t ] = EffecKTPkvi(koef,Te)
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