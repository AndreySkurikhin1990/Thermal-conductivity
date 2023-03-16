%Серое приближение для вермикулита
%function t = tmp66(vyfv,vysv,vyuv,vmivmf)
function t = tmp66()
delT=1e2; Tna=1e2; te0=273.15; cemv=11; Te=Tna:delT:((cemv-1)*delT+Tna);
Te=Te+te0;
y0=30e-3; 
%qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); kqt=polyfit(te,qo,2); r=polyval(kqt,858.34); qo=0;
vyfv=0
vysv=0
vyuv=0
%vmivmf=0;
for vmivmf=0:2
vmivmf=vmivmf'
temvci=napMasEKTPVer(vyfv,vysv,vmivmf,Te,0,y0,vyuv);
temvhi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,1,y0,vyuv);
qoi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,2,y0,vyuv);
ektpvi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,3,y0,vyuv);
temvsi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,4,y0,vyuv);
%---
temvc=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,0);
temvh=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,1);
qo=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,2);
temvs=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,3);
ektpv=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,4);
n=length(temvc);
%---
for w=1:n
    Tk=temvc(w);
    T0=temvh(w);
    qob=qo(w);
    ektp=ektpv(w);
    T3=temvs(w);
    Ts=(Tk+T0)/2
f=opredBetaPBetaMVer(T0,Tk,y0,ektp,T0,qob,T3);
end
tem=[28,94,193,342,491,641,791,991]; 
tem=tem+te0; n=length(tem);
for w=1:n
    f(w)=opredKTPTKTochVer(ektpv, temvs, tem(w));
end
end
ektpv=ektpv';
Te=Te';
f=f';
tem=tem';
t=0;
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

function [ t ] = opredBetaPBetaMVer(T0,Tk,tol,koeftep,Tz,qo,T3)
kt=polyfit([0,tol/2,tol],[T0,T3,Tk],2);
dko2m=[5,5,6,13,20,28]/1e2; dko2t=6e2:2e2:16e2; 
for k=1:length(dko2m)
    dko2m(k)=1-dko2m(k);
end
dko2=opredKTPTKTochVer(dko2m, dko2t, (T0+Tk)/2);
Nt=2e1+1; y=0:tol/Nt:tol; n=length(y);
for k=1:n
    Te(k)=polyval(kt,y(k));
end
Te=Te';
wmg=217e-3; wal=132e-3; wsi=368e-3;
koal=rasKoefAlpha(Te,wmg,wsi,wal); 
s=mean(koal);
sig=5.67e-8; %alf=sredkopoPatch(tol,Tna,Tko)
dv=RasshDiapDlinVoln();
knuVer=RasMasKoAbs();
dko=0.444905; 
dko3=1.3002913762; 
s=s*dko*dko2*dko3;
knuVer=s*knuVer;
npp=Kramers_n();
for k=1:n
    alsrm(k)=knusreddvvsrVer(Te(k),knuVer,dv,npp);
end
alsr=trapz(Te,alsrm)/(Tk-T0); alsr=abs(alsr)
crp=ProvAd(koeftep*alsr/4/sig/(Tz^3))/pi
eps=epssr([T0,Tk],dv,npp,knuVer);
kctp=PoiskChisKTP(qo,tol,T0,Tk,sig,alsr/dko,eps)
koeftep=koeftep'
tau0=alsr*tol;
Nc=2*ceil(tau0); 
tau00=tau0/Nc;
%tau00=tau0
tol1=tau00/alsr;
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
npr=knusreddvvsrVer(Tz,npp,dv,npp);
npr2=knusreddvvsrVer(Tz,npp2,dv,npp);
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
%disp('First');
%bpbm=reshenbpbm(pi, E21, E21, pi, En2t40-E01i,  En2t41-E11i);
%bp=ProvAd(bpbm(1)); bm=ProvAd(bpbm(2));
%qr=abs(rasqrlamtep(bp,bm,Tz,sig,tau00,npr2,tetataus4,taus));
%-------
%reiz_drm=DruMetOprqr(tetataus4,taus,tau00,npr2,Tz,sig,E21,En2t41,E11i)
ri=DruMetOprqr2(tau0,((T0/Tz)^4)/pi,((Tk/Tz)^4)/pi,Tk/Tz,T0/Tz,crp,T3/Tz,npp2,Tz,npr2,dv,npp);
%reiz=opredphimalqbol(tau00,sig,Tk,T0);
%bpl=bp*sig*(Tz^4); bmi=bm*sig*(Tz^4); tetatausz=tetataus*Tz; tolt=tol/Nt;
%reizm=resuIzlu(tolt,T0,Tk,koal,dv,knuVer/s,npp,tetatausz,bpl,bmi,tau00,npr,length(tetataus),Te,ksi,0,dko,dko3)
%reiz=reizm(1); Iplu=reizm(2); Imin=reizm(3);
%reiz2=resuIzlu2(tol,T0,Tk,koal,dv,knuVer/s,npp,tetatausz,sig*(T0^4),sig*(Tk^4),tau0,npr,length(tetataus),Te,ksi,0,dko,dko3,dko2m,dko2t)
Nt=2; taus=0:tau00/Nt:tau00; n=length(taus); 
kna=(Tk-T0)/tau0/Tz; kb=T0/Tz;
for k=1:n
    tetataus(k)=kna*taus(k)+kb;
    tetataus4(k)=tetataus(k)^4;
end
%bplbmi=opredqrezVerm(npp,dv,tol,knuVer,Te,koal,length(y),y,T0,Tk,T3)
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

function [ epsi ] = epssr(tem,dv,npp,alf)
n=length(npp); p=length(tem);
np=0; hi=0; ep=0;
for j=1:p
    eps(j)=0;
for k=1:n
hi=dv(k)*alf(k)/4/pi/npp(k);
np=(npp(k)^2+hi^2)^0.5;
ep(k)=emdiel(np);
%ep(k)=emmet(npp(k),hi);
end
eps(j)=epssred(dv,ep,tem(j),npp);
end
epsi=eps;
end

function [ res ] = reshenbpbm(a11,a12,a21,a22,b1,b2)
matko=[a11, a12; a21, a22]; matsvch=[b1,b2];
res=(inv(matko)*matsvch')' ;
end

function [ rezizl ] = resuIzlu(tol,T0,Tk,koal,dl,alfs,npp,temp,Iplus,Iminus,tau0,npr,Nt,Tet,ksi,vybor,dko,dko2)
p=length(dl); 
sigma=5.67e-8; %dko=1e0; %учет пористой структуры %dko2=65e-2; %учет КП от КО
koal=koal*dko*dko2;
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
I1g(j)=2*integ2ma(ksi1,I11g);
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
I2g(j)=2*integ2ma(ksi2,I22g);
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
s=mean(qre);
sp=mean(Iplusi); 
sm=mean(Iminusi);
FKhN=sigma*(npr^2)*abs((T0^4)-(Tk^4))
rezizl=[s,sp,sm];
end

function [ rezizl ] = resuIzlu2(tol,T0,Tk,koal,dl,alfs,npp,temp,Iplus,Iminus,tau0,npr,Nt,Tet,ksi,vybor,dko,dko2,dko2m,dko2t)
p=length(dl); 
sigma=5.67e-8; %dko=1e0; %учет пористой структуры %dko2=65e-2; %учет КП от КО
koal=koal*dko*dko2;
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
altek=alfs(m)*opredKTPTKTochSha(koal,Tet,temp(k))*opredKTPTKTochSha(dko2m,dko2t,temp(k));
altej=alfs(m)*opredKTPTKTochSha(koal,Tet,temp(j))*opredKTPTKTochSha(dko2m,dko2t,temp(k));
eb1(m)=altek*c1t/(lam^5)/me;
it1(m)=eb1(m)*ProvAd(integroexpon(2,abs(altej*ksi(j)-altek*ksi(k))));
    end
    ksi1(k)=ksi(k);
    I11g(k)=integ2ma(dv,it1);
end
I1g(j)=2*integ2ma(ksi1,I11g);
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
altek=alfs(m)*opredKTPTKTochSha(koal,Tet,temp(k))*opredKTPTKTochSha(dko2m,dko2t,temp(k));
altej=alfs(m)*opredKTPTKTochSha(koal,Tet,temp(j))*opredKTPTKTochSha(dko2m,dko2t,temp(j));
eb2(m)=altek*c1t/(lam^5)/me;
taut=abs(altek*ksi(k)-altej*ksi(j));
it2(m)=eb2(m)*ProvAd(integroexpon(2,taut));
    end
    ksi2(q)=ksi(k);
    I22g(q)=integ2ma(dv,it2);
    q=q+1;
end
I2g(j)=2*integ2ma(ksi2,I22g);
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
s=0; sp=0; sm=0;
for j=1:Nt
qre(j)=Iplusi(j)-Iminusi(j)+I1g(j)-I2g(j);
s=s+qre(j);
sp=sp+Iplusi(j);
sm=sm+Iminusi(j);
end
FKhN=sigma*(npr^2)*abs((T0^4)-(Tk^4))
if (vybor==0)
rezizl=s/Nt;
elseif (vybor==1)
rezizl=sp/Nt;
elseif (vybor==2)
rezizl=sm/Nt;
end
end

function inte = integ2ma(ar1,ar2)
p=length(ar1);
if (p==1)
    inte=0;
else
su=0;
for k=2:p
    su=su+(ar2(k)+ar2(k-1))*(ar1(k)-ar1(k-1))/2;
end
inte=su;
end
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

function oipim = opredqrezVerm(npp,dl,tol,alfs,temp,koal,Nt,ksi,T1,T2,T3)
p=length(dl); nk=length(koal); nk2=nk/2; Tz=T1;
for k=1:p
    dv(k)=dl(k)/npp(k);
end
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=2*pi*PP*(c0^2);
C2=PP*c0/PB;
sig=5.67e-8;
%1----------q_r(0), q_r(tau0nu)
for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
        lam=dv(m);
        tau0nu=(alfs(m)*koal(nk))*ksi(nk);
        tau0nu2=(alfs(m)*koal(nk2))*ksi(nk2);
        tau00nu=(alfs(m)*koal(2))*ksi(2); %на горячей стороне
        tau01nu=tau0nu-tau00nu; %на холодной стороне
        tau02nu=(alfs(m)*koal(nk2+1))*ksi(nk2+1); %в середине образца
        tau=0; it1=0; it2=0; it22=0; it21=0;
        for j=1:Nt
tem=temp(j);
me=exp(c2t/lam/tem)-1;
taus=(alfs(m)*koal(j))*ksi(j);
tau(j)=taus;
it10(j)=c1t/(lam^5)/me;
it1(j)=it10(j)*ProvAd(integroexpon(2,taus));
if (tau0nu>=taus)
it2(j)=it10(j)*ProvAd(integroexpon(2,tau0nu-taus));
else it2(j)=0;
end
if (tau00nu<=taus)
it22(j)=it10(j)*ProvAd(integroexpon(2,taus-tau00nu));
it21(j)=0;
else it22(j)=0;
it21(j)=it10(j)*ProvAd(integroexpon(2,tau00nu-taus));
end
        end
        it31m(m)=ProvAd(integroexpon(3,tau00nu));
        if (tau01nu>=0)
        it32m(m)=ProvAd(integroexpon(3,tau01nu));
        else it32m(m)=0;
        end
        E3m(m)=ProvAd(integroexpon(3,tau0nu));
        E32m(m)=ProvAd(integroexpon(3,tau0nu2));
    I11g(m)=integ2ma(tau,it1);
    I12g(m)=integ2ma(tau,it2);
    I21g(m)=integ2ma(tau,it21);
    I22g(m)=integ2ma(tau,it22);
end
E3i=2*pi*integ2ma(dv,E3m);
E32i=2*pi*integ2ma(dv,E32m);    
It31i=2*pi*integ2ma(dv,it31m);
It32i=-2*pi*integ2ma(dv,it32m);
I1g=-2*integ2ma(dv,I11g);
I2g=2*integ2ma(dv,I12g);
I21gi=2*integ2ma(dv,I21g);
I22gi=-2*integ2ma(dv,I22g);
I212gi=I21gi+I22gi;
%2----------
for m=1:p
        c1t=C1/(npp(m)^2);
        c2t=C2/npp(m);
        lam=dv(m);
        tau0nu=(alfs(m)*koal(nk))*ksi(nk);
        tau0nu2=(alfs(m)*koal(nk2))*ksi(nk2);
        tau=0; it10=0; it1=0; I11g=0; I12g=0;
for j=1:nk2
    tem=temp(j);
    me=exp(c2t/lam/tem)-1;
    taus=alfs(m)*koal(j)*ksi(j);
    tau(j)=taus;
    it10(j)=c1t/(lam^5)/me;
if (tau0nu2>=taus)
it1(j)=it10(j)*ProvAd(integroexpon(2,tau0nu2-taus));
else it1(j)=0;
end
end
    I11g(m)=integ2ma(tau,it1);
    tau=0; q=1; it10=0; it1=0;
for j=nk2:Nt
    tem=temp(j);
    me=exp(c2t/lam/tem)-1;
    taus=alfs(m)*koal(j)*ksi(j);
    tau(q)=taus;
    it10(q)=c1t/(lam^5)/me;
if (tau0nu2<=taus)
it1(q)=it10(q)*ProvAd(integroexpon(2,-tau0nu2+taus));
else it1(q)=0;
end
q=q+1;
end
    I12g(m)=integ2ma(tau,it1);
end
I3g=2*integ2ma(dv,I11g);
I4g=-2*integ2ma(dv,I12g);
I34g=I3g+I4g;
Nt=length(ksi);
Ite=oprIpImVer(dv,alfs,temp,koal,npp,tol,ksi,Nt,T1,T2,T3,dl,E32i,E3i,I212gi,I1g/4/sig/(Tz^4),I2g/4/sig/(Tz^4),I34g/4/sig/(Tz^4))
Ip=Ite(1); Im=Ite(2);
q0=pi*Ip-E3i*Ip+I1g
q2=E3i*Ip-pi*Im+I2g
q1=E32i*Ip-E32i*Im+I34g
Ip=Ite(3); Im=Ite(4);
q0=pi*Ip-E3i*Ip+I1g
q2=E3i*Ip-pi*Im+I2g
q1=E32i*Ip-E32i*Im+I34g
%Total-----------
%tau0=opredtau0nu(koal,alfs,tol,dl,T1,T2,temp)
oipim=0;
end

function [ Ite ] = oprIpImVer(dv,alfs,temp,koal,npp,tol,ksi,Nt,T1,T2,T3,dl,e32i,e3i,i212gi,it1is,it2is,i34gs)
Tz=T1; teta2=T2/Tz; teta3=T3/Tz; sig=5.67e-8; e32i=e32i/sig/(Tz^4)/4; e3i=e3i/sig/(Tz^4)/4;
p=length(dv); nk=length(koal); nk2=nk/2; %disp(length(ksi));
for j=1:Nt
    n2m(j)=n2sredVer(dl,npp,temp(j));
    teta(j)=temp(j)/Tz;
    teta4(j)=(teta(j))^4;
    n2t4(j)=n2m(j)*teta4(j);
end
n2m=n2m';
n2i=integ2ma(temp,n2m)/abs(T2-T1)
n2i=abs(n2i);
tau0=integ2ma(dv,alfs)*koal(nk)*tol;
tau02=integ2ma(dv,alfs)*koal(nk2)*tol/2;
%------
for j=1:Nt
taum=0;
for m=1:p
        taum(m)=alfs(m)*koal(j)*ksi(j);
end
tau(j)=integ2ma(dv,taum);
E2m(j)=ProvAd(integroexpon(2,tau(j)));
if (tau0>=tau(j))
E20m(j)=ProvAd(integroexpon(2,tau0-tau(j)));
else E20m(j)=0;
end
it1(j)=E2m(j)*n2m(j)*(teta4(j));
it2(j)=E20m(j)*n2m(j)*(teta4(j));
end
it1i=-2*integ2ma(tau,it1)/4;
%it1i=it1is;
it2i=2*integ2ma(tau,it2)/4;
%it2i=it2is;
%------
it1=0; taun=0;
for j=1:Nt/2
    taus=tau(j); taun(j)=taus;
if (tau02>=taus)
it1(j)=n2m(j)*ProvAd(integroexpon(2,tau02-taus))*teta4(j);
else it1(j)=0;
end
end
I3g=2*integ2ma(taun,it1)/4;
q=1; it1=0; taun=0;
for j=Nt/2:Nt
    taus=tau(j); taun(q)=taus;
if (tau02<=taus)
it1(q)=n2m(j)*ProvAd(integroexpon(2,taus-tau02))*(teta4(j));
else it1(q)=0;
end
q=q+1;
end %disp(length(taun)); disp(length(it1));
I4g=-2*integ2ma(taun,it1)/4;
I34g=I3g+I4g;
%I34g=i34gs;
%------
for j=1:Nt
E10m(j)=n2m(j)*teta4(j)*ProvAd(integroexpon(1,tau(j)));
E11m(j)=n2m(j)*teta4(j)*ProvAd(integroexpon(1,abs(tau02-tau(j))));
E12m(j)=n2m(j)*teta4(j)*ProvAd(integroexpon(1,abs(tau0-tau(j))));
end
E10i=integ2ma(tau,E10m);
E11i=integ2ma(tau,E11m);
E12i=integ2ma(tau,E12m);
ep0=pi/sig/(Tz^4)/4;
em0=4*ep0*ProvAd(integroexpon(2,tau0));
em1=4*ep0*ProvAd(integroexpon(2,tau02));
ep1=(4*ep0+em1)/2;
em2=(em0+em1)/2;
n2t4s0=(n2t4(1)+n2t4(nk2))/2;
n2t4s1=(n2t4(nk)+n2t4(nk2))/2;
E101i=n2t4s0-(E10i+E11i)/2/2;
E112i=n2t4s1-(E11i+E12i)/2/2;
e332i=(e3i-e32i)/(tau02);
ep032i=(e32i-ep0)/(tau02);
I2i34=(it2i-I34g)/(tau02);
I341i=(I34g-it1i)/(tau02);
e332i=e332i+em2/2;
ep032i=ep032i+ep1/2;
E101i=E101i-I341i;
E112i=E112i-I2i34;
disp('Boleye tochnoye resheniye');
Ipm=BolyeTochReshIpIm(n2t4(1),n2t4(nk2),n2t4(nk),ep0,em0,em1,E10i,E11i,e32i,I34g,e3i,it1i,it2i,tau02,T1,T2)
Ipm=[ProvAd(Ipm(1)),ProvAd(Ipm(2))]; Ipmn=[Ipm(1),Ipm(2)];
disp('Priblizhennoye resheniye');
Ipm=reshSysUrav(ep032i,e332i,e332i,ep032i,E101i,E112i);
Ipm=[ProvAd(Ipm(1)),ProvAd(Ipm(2))];  
Ipmn(3)=Ipm(1); Ipmn(4)=Ipm(2); Ipl=Ipm(1); Imi=Ipm(2);
dqdt0=n2m(1)*teta4(1)-(ep0*Ipl+em0*Imi+E10i)/2
dqdt2=n2m(nk2)*teta4(nk2)-(em1*Ipl+em1*Imi+E11i)/2
dqdt1=n2m(nk)*teta4(nk)-(em0*Ipl+ep0*Imi+E12i)/2
Ite=Ipmn;
end

function [ x ] = reshSysUrav(a11,a12,a21,a22,b1,b2)
mas=[a11,a12;a21,a22]; b=[b1,b2];
xn=inv(mas)*b';
xn=[ProvAd(xn(1)),ProvAd(xn(2))];
x=xn;
end

function ns = n2sredVer(dv,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458; vl=c0; dl=0;
for k=1:length(npp)
    vl=vl/npp(k);
    c1=2*pi*PP*(vl^2);
    c2=PP*vl/PB;
    lambda=dv(k)/npp(k);
    Ib(k)=c1/(lambda^5)/(exp(c2/lambda/tem)-1);
    Ibn(k)=(npp(k)^2)*Ib(k);
    vl=c0;
    dl(k)=lambda;
end
nc=trapz(dl,Ibn);
nz=trapz(dl,Ib);
ns=abs(real(nc/nz));
end

function t = opredtau0nu(koal,alfs,tol,dv,T0,Tk,tem)
tau0a=0; tau0b=1e12; q=0; nq=1e4; eps=1e-4; nit=length(koal); ra=1e1;
while ((q<nq) && (ra>eps))
    tau0c=(tau0a+tau0b)/2;
    tau0ma=0:tau0a/nit:tau0a; ktaa=(Tk-T0)/tau0a; ktba=T0;
    tau0mb=0:tau0b/nit:tau0b; ktab=(Tk-T0)/tau0b; ktbb=T0;
    tau0mc=0:tau0c/nit:tau0c; ktac=(Tk-T0)/tau0c; ktbc=T0;
    fua=0;
    for k=1:length(tau0ma)
        fum=0;
        tema=ktaa*tau0ma(k)+ktba;
        koala=opredKTPTKTochVer(koal, tem, tema);
        for m=1:length(dv)
        fum(m)=1/(koala*alfs(m));
        end
        fua(k)=integ2ma(dv,fum);
    end
    fub=0;
    for k=1:length(tau0mb)
        fum=0; 
        temb=ktab*tau0mb(k)+ktbb;
        koalb=opredKTPTKTochVer(koal, tem, temb);
        for m=1:length(dv)
        fum(m)=1/(koalb*alfs(m));
        end
        fub(k)=integ2ma(dv,fum);
    end
    fuc=0;
    for k=1:length(tau0mc)
        fum=0;
        temc=ktac*tau0mc(k)+ktbc;
        koalc=opredKTPTKTochVer(koal, tem, temc);
        for m=1:length(dv)
        fum(m)=1/(koalc*alfs(m));
        end
        fuc(k)=integ2ma(dv,fum);
    end
    fa=integ2ma(tau0ma,fua)-tol;
    fb=integ2ma(tau0mb,fub)-tol;
    fc=integ2ma(tau0mc,fuc)-tol;
    la = VybCAB(fa,fb,fc,tau0a,tau0b,tau0c);
    tau0a=la(1); tau0b=la(2); 
    ra=abs(tau0b-tau0a);
    q=q+1;
end
t=0;
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
ppo(j)=integ2ma(tauss,I1);
else ppo(j)=0;
end
end
vpo=integ2ma(taus,ppo);
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

function t = DruMetOprqr2(tau0,betap,betam,teta2,teta0,Ncrp,teta3,npp2,Tz,npr2,dv,npp)
npr2=npr2';
Nt=1e2;
taus=0:tau0/Nt:tau0;
kt=polyfit([0,tau0/2,tau0],[teta0,teta3,teta2],2);
n=length(taus);
n2=round(n/2);
for k=1:n
    teta(k)=polyval(kt,taus(k));
    tetataus4(k)=(teta(k))^4;
    np2m(k)=knusreddvvsrVer(teta(k)*Tz,npp2,dv,npp);
end
npr2=integ2ma(teta,np2m)/(teta2-teta0); npr2=abs(npr2)
%2-------
r=1;
for j=1:n
tauss=0; I1=0;
for k=1:j
I1(k)=npr2*tetataus4(k);
tauss(k)=taus(k);
end
q=1;
if (j>=n2)
for k=n2:j
    tauss1(q)=taus(k);
    I11(q)=npr2*tetataus4(k);
    q=q+1;
end
end
if (j>1)
ppo(j)=integ2ma(tauss,I1);
else ppo(j)=0;
end
if (j>n2)
    ppo1(j)=integ2ma(tauss1,I11);
else ppo1(j)=0;
end
if (j<=n2)
    taus3(j)=taus(j);
    ppo3(j)=ppo(j);
else
        taus3k(r)=taus(j);
        ppo3k(r)=ppo1(j);
        r=r+1;
end
end
vpo=integ2ma(taus,ppo);
vpo3=integ2ma(taus3,ppo3);
vpo3k=integ2ma(taus3k,ppo3k);
%-----
E4p=ProvAd(integroexpon(4,tau0));
E4pk=-(tau0/2+E4p-1/3)*betap/2; %I2
E3m=ProvAd(integroexpon(3,tau0));
E4mk=-(1/3-E4p-tau0*E3m)*betam/2; %I3
I1=0; q=1;
for k=1:n
I1(k)=npr2*tetataus4(k);
E3ma(k)=ProvAd(integroexpon(3,abs(tau0-taus(k))))*I1(k);
E2ma(k)=ProvAd(integroexpon(2,abs(tau0-taus(k))))*I1(k);
if (k<=n2)
        I13(k)=I1(k);
        E3ma30(k)=ProvAd(integroexpon(3,tau0/2-taus(k)))*I1(k);
        E2ma30(k)=ProvAd(integroexpon(2,tau0/2-taus(k)))*I1(k);
else
    taus31(q)=taus(k);
    I131(q)=I1(k);
    E2ma3(q)=ProvAd(integroexpon(2,taus(k)))*I1(k);
    E2ma31(q)=ProvAd(integroexpon(2,taus(k)-tau0/2))*I1(k); 
    E2ma31k(q)=E2ma31(q);
    E2ma31(q)=E2ma31(q)-ProvAd(integroexpon(2,taus(k)))*I1(k);
    E3ma3(q)=ProvAd(integroexpon(3,abs(taus(k)-tau0/2)))*I1(k);
    E3ma31k(q)=E2ma31(q);
    E3ma3(q)=E3ma3(q)-ProvAd(integroexpon(3,taus(k)))*I1(k);
    E3ma32(q)=ProvAd(integroexpon(3,tau0-taus(k)))*I1(k);
    q=q+1;
end
end
%ppok=ppo(length(ppo));
E3ma31ki=-integ2ma(taus31,E3ma31k);
E2ma31ki=-integ2ma(taus31,E2ma31k);
E3ma32ik=integ2ma(taus31,E3ma32);
I131i=integ2ma(taus31,I131);
I1i=integ2ma(taus,I1);
E3mi=integ2ma(taus,E3ma);
E2mi=integ2ma(taus,E2ma);
E3mik=E3mi; 
E2maik=E2mi; 
E3mi=(vpo-I1i/2+E3mi)/2/pi; %I4
E3mi=vpo/pi+E4pk+E4mk+E3mi;
disp('d_teta/d_tau(0)');
E3mi=-(E3mi-Ncrp*(teta2-teta0))/tau0/Ncrp; %d_teta/d_tau(0)
disp('d_teta/d_tau(tau0)');
E3mi=(I1i/2/pi-betap*(1/2-E3m)/2-betam*(1/2-E3m)/2+E2mi/2/pi)/Ncrp+E3mi; %d_teta/d_tau(tau0)
%----------
%betapl=betap;
%betami=betam;
I1i=-integ2ma(taus31,I13)/2;
I3i0=integ2ma(taus3,E3ma30); 
E3ma30i=-I3i0;
I3i1=integ2ma(taus31,E3ma3);
I2i=integ2ma(taus31,E2ma3);
I2i=-tau0*I2i/2;
E3mi=vpo3+I1i+I3i0+I3i1+I2i;
E3mi=-E3mi/2/pi;
disp('d_teta/d_tau(0)');
betam=(teta3^4)/pi;
E4p=ProvAd(integroexpon(4,tau0/2));
E4p0=ProvAd(integroexpon(4,tau0));
E4pk=-(tau0/2/2+E4p-1/3)*betap/2; %I2
E3m=ProvAd(integroexpon(3,tau0));
E4mk=-(E4p-E4p0-tau0*E3m/2)*betam/2; %I3
E3mi=vpo3/pi+E4pk+E4mk+E3mi;
E3mi=E3mi-Ncrp*(teta3-teta2);
E3mi=-E3mi/Ncrp/(tau0/2); %d_teta/d_tau(0)
%---
disp('d_teta/d_tau(tau0/2)');
I1i=trapz(taus3,I13)/2/pi;
E3p=ProvAd(integroexpon(3,tau0/2));
E3pk=E3p;
E3p=-(1/2-E3p)*betap/2;
E3mk=E3m;
E3m=E3p-E3m;
E3m=-E3m*betam/2;
I2i=integ2ma(taus3,E2ma30); 
E2ma30i=I2i; I2i=I2i/2/pi;
I2i1=integ2ma(taus31,E2ma31)/2/pi;
E3mi=(I2i+I2i1+E3m+E3p+I1i)+Ncrp*E3mi;
E3mi=E3mi/Ncrp; %d_teta/d_tau(tau0/2)
disp('d_teta/d_tau(tau0/2)');
%E3mi=E3mik+E3ma30i+E3ma31ki+(tau0/2)*(E2ma30i+E2ma31ki);
E3mi=E3ma32ik;
E3mi=-E3mi/2/pi;
I3i=1/3-tau0*E3pk/2-E4p;
I3i=-I3i*betam/2;
I2i=tau0*E3pk/2+E4p0-E4p;
I2i=-I2i*betap/2;
E3mi=vpo3k/pi+I2i+I3i+E3mi;
E3mi=E3mi-Ncrp*(teta2-teta3);
E3mi=-E3mi/Ncrp/(tau0/2) %d_teta/d_tau(tau0/2)
E3mik=E3mi;
disp('d_teta/d_tau(tau0)');
I2i=-(-E3mk+E3p)*betap/2;
I3i=-(1/2-E3p)*betam/2;
E3mi=E2ma31ki+E2maik;
E3mi=I131i/pi+I3i+I2i-E3mi/2/pi;
E3mi=E3mi/Ncrp+E3mik
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
Q(k)=sig*abs(T1^4-T0^4)*qbol(k); %qbol=2Q, q=(I+-I-)Q
st4(k)=(abs(T1^4-T0^4)*phimal(k)+(T1^4))^(1/4);
s=s+Q(k);
p=p+st4(k);
end
disp('q R sred');
s=s/n
disp('T sred');
p=p/n
KhrNus=((5.67e-8)*(n^2)*abs(T1^4-T0^4))
t=0;
end

function t = vychnevyaz(a1,a2)
n=length(a1); s=0;
for k=1:n
    s=s+abs(a1(k)-a2(k));
end
t=s;
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

%определяет ПТП, КТП, температуры для вермикулита
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

function [ vm ] = BolyeTochReshIpIm(n2t4n,n2t4s,n2t4k,ep0,em0,em1,E10i,E11i,e32i,I34g,e3i,it1i,it2i,tau02,Tn,Tk)
sig=5.67e-8;
Ipl=sig*(Tn^4);
Imi=sig*(Tk^4);
e332i=(e3i-e32i)/(tau02);
ep032i=(e32i-ep0)/(tau02);
I341i=(I34g-it1i)/(tau02);
I2i34=(it2i-I34g)/(tau02);
nit=1e2;
k=0;
eps=1e-4;
ra=1e2;
while ((k<nit) && (ra<eps))
arglc1=ep032i*Ipl+e332i*Imi+I341i;
argpc11=n2t4n-(4*ep0*Ipl+em0*Imi+E10i)/2;
argpc12=n2t4s-(em1*(Ipl+Imi)+E11i)/2;
arglc2=e332i*Ipl+ep032i*Imi+I2i34;
argpc21=n2t4k-(em0*Ipl+4*ep0*Imi+E12i)/2;
argpc22=argpc12;
alc1=1/sqrt(1+arglc1^2);
alc2=1/sqrt(1+arglc2^2);
apc11=1/sqrt(1+argpc11^2);
apc12=1/sqrt(1+argpc12^2);
apc21=1/sqrt(1+argpc21^2);
apc22=1/sqrt(1+argpc22^2);
dfdx11=ep032i*alc1+ep0*apc11+em1*apc12/4;
dfdx12=e332i*alc1+em0*apc11/4+em1*apc12/4;
dfdx21=e332i*alc2+em0*apc21/4+em1*apc22/4;
dfdx22=ep032i*alc2+ep0*apc21+em1*apc22/4;
A=[dfdx11,dfdx12;dfdx21,dfdx22];
f1=atan(arglc1);
f1=f1-atan(n2t4n-(4*ep0*Ipl+em0*Imi+E10i)/2)/2;
f1=f1-atan(n2t4s-(em1*(Ipl+Imi)+E11i)/2)/2;
f2=atan(arglc2);
f2=f2-atan(n2t4k-(em0*Ipl+4*ep0*Imi+E12i)/2)/2;
f2=f2-atan(n2t4s-(em1*(Ipl+Imi)+E11i)/2)/2;
f=[f1,f2];
dx=-inv(A)*f';
Ipl=Ipl+dx(1);
Imi=Imi+dx(2);
ra=sqrt(dx(1)^2+dx(2)^2);
k=k+1;
end
Ipm=[ProvAd(Ipl),ProvAd(Imi)]; 
Ipl=Ipm(1); Imi=Ipm(2); 
Ipm=[Ipl,Imi];
vm=Ipm;
end