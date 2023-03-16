%Определение tau0
function tt = tmp74()
format longg;
te0=273.15; delT=1e2; Tna=1e2+te0; n=13; Te(1)=Tna;
for k=2:n
    Te(k)=Te(k-1)+delT;
end
Te=Te'
y0=3e1*1e-3; vyuv=0;
nvyfv=[0,1]; nvysv=[0,1,2]; nvmivmf=[0,1,2]; nvyuv=[1,2];
for kvf=1:length(nvyfv)
    vyfv=nvyfv(kvf)
    for jvsv=1:length(nvysv)
        vysv=nvysv(jvsv)
        if (vyfv==0)
            for qvmi=1:length(nvmivmf)
                vmivmf=nvmivmf(qvmi)
                c=funcRaschTau0(Te,vyfv,vysv,vmivmf,y0,vyuv,te0);
            end
        elseif (vyfv==1)
            for qvuv=1:length(nvyuv)
                vyuv=nvyuv(qvuv)
                c=funcRaschTau0(Te,vyfv,vysv,vmivmf,y0,vyuv,te0);
            end
        end
    end
end
tt=c;
end

function tt = funcRaschTau0(Te,vyfv,vysv,vmivmf,y0,vyuv,te0)
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
    Ts(w)=(Tk+T0)/2;
    f(w)=opredBetaPBetaMVer(T0,Tk,y0,ektp,T0,qob,T3,te0);
end
Ts=Ts'
f=f'
tt=0;
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
    tetan(k)=Te(k)/Tz;
end
Te=Te';
pkp=[1,0.971325109,0.937848799,0.836713218,0.757732344,0.735010914,0.725526014,0.717549581,0.70984116];
tpkp=(2e2:1e2:1e3)+te0;
wmg=217e-3; 
wal=132e-3; 
wsi=368e-3;
koep=rasKoefAlpha(Te,wmg,wsi,wal); %sig=5.67e-8; %alf=sredkopoPatch(tol,Tna,Tko)
dv=RasshDiapDlinVoln();
knuVer=RasMasKoAbs();
dko=0.444905; %из статистического моделирования
dko3=1.3002913762; %из пористой структуры самого вермикулита
knuVer=dko*dko3*knuVer;
betam=knuVer/dko; %КО
npp=Kramers_n();
%for k=1:length(npp)
    %npp2(k)=npp(k)^2;
%end
%alsr=knusreddvvsrVer(Tz,knuVer,dv,npp,pkp,tpkp);
eps=epssr(Te,dv,npp,betam,koep,dko2m,dko2t,pkp,tpkp); %падение температуры и экспериментальные данные
f=1; e=1e-6;
for k=1:length(eps)
    if (eps(k)<e)
        f=-1;
        break; 
    end
end
%bp=(T0/Tz)^4;
%bp=bp/pi; 
%bm=(Tk/Tz)^4; 
%bm=bm/pi;
for k=1:n
    als(k)=knusreddvvsrVer(Te(k),knuVer,dv,npp,pkp,tpkp);
    %betas(k)=knusreddvvsrVer(Te(k),betam,dv,npp,pkp,tpkp);
    taus(k)=als(k)*y(k);
    %tetataus(k)=Te(k)/Tz;
    %tetataus4(k)=(tetataus(k))^4;
    %npr(k)=knusreddvvsrVer(Te(k),npp,dv,npp,pkp,tpkp);
    %npr2(k)=knusreddvvsrVer(Te(k),npp2,dv,npp,pkp,tpkp);
end
n=length(taus);
tau0=taus(n);
tau0=opredtau0(taus,als,tol,n-1,tau0);
%tau0=taus(n);
%koal=vychKoe(tpkp,Te,pkp)
%tau0=opredtau0uss(koal,knuVer,tol,dv,T0,Tk,Te,T3,tau0)
t=tau0;
return
if (f>0)
kctp=PoiskChisKTP(qo,tol,T0,Tk,sig,betas,eps,pkp,tpkp,Te);
if (koeftep<kctp)
    kctp=0;
end
crp=ProvAd(kctp*alsr/4/sig/(Tz^3))/pi;
for k=1:n
    pc2=taus(k);
    pc2=ProvAd(integroexpon(2,pc2));
    pc2=bp*pc2/2;
    pc3=tau0-taus(k);
    if (pc3>=0)
	pc3=ProvAd(integroexpon(2,pc3));
    pc3=bm*pc3/2;
    else
        pc3=0;
    end
    nev2(k)=npr2(k)*tetataus4(k)/pi-pc2-pc3-RasFunKvadMark(tau0,tetataus4(k),npr2(k),taus(k))-crp*2*kt(1);
end
nev2=nev2';
%t=postgraf(y,nev2,y,tetataus);
end
end

function f = postgraf(koox,tetaxx,kooy,tetaxy)
pl=plot(koox,tetaxx,'-b',kooy,tetaxy,'-k');
legend('Невязка','Решение','location','best');
set(pl,'LineWidth',3);
hold on;
grid on;
xlabel({'Координата, м'});
ylabel({'Невязка'});
title({'График зависимости невязки от координаты'});
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
nc=trapz(dl,Ibc);
nz=trapz(dl,Ib);
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

function m = RasFunKvadMark(tau0,t4,n2,ts)
x=zadKoef(0);
A=zadKoef(1);
a=0;
b=tau0; %delta=2;
n=length(x);
I=0;
for k=1:n
    tau=(a+b)+(b-a)*x(k);
    tau=tau/2;
    An(k)=A(k); %An(k)=(b-a)*An(k)/delta;
    I=I+An(k)*funcRas(t4,n2,tau,ts);
end
m=I/2/pi;
end

function f = funcRas(teta4,n2,tau,taus)
ab=abs(tau-taus);
inex=ProvAd(integroexpon(1,ab));
f=n2*teta4*inex;
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

function t = opredtau0(tau,als,L,nit,tau0)
tau0a=1e-1*tau0; tau0b=1e1*tau0; 
q=0; nq=1e4; eps=1e-2; ra=1e1;
while ((q<nq) && (ra>eps))
    tau0c=(tau0a+tau0b)/2;
    tau0ma=0:tau0a/nit:tau0a; 
    tau0mb=0:tau0b/nit:tau0b;
    tau0mc=0:tau0c/nit:tau0c; 
for k=1:length(tau0ma)
        alpa(k)=opredKTPTKTochVer(als, tau, tau0ma(k));
        fuma(k)=1/alpa(k);
end
fa=trapz(tau0ma,fuma)-L;
for k=1:length(tau0mb)
        alpb(k)=opredKTPTKTochVer(als, tau, tau0mb(k));
        fumb(k)=1/alpb(k);
end
fb=trapz(tau0mb,fumb)-L;
for k=1:length(tau0mc)
        alpc(k)=opredKTPTKTochVer(als, tau, tau0mc(k));
        fumc(k)=1/alpc(k);
end
fc=trapz(tau0mc,fumc)-L;
la=VybCAB(fa,fb,fc,tau0a,tau0b,tau0c);
tau0a=la(1); 
tau0b=la(2); 
ra=abs(tau0b-tau0a);
q=q+1;
end
t=tau0c;
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

function t = opredtau0uss(koal,alfs,tol,dv,T0,Tk,tem,T3,tau0)
tau0a=1e-2*tau0; tau0b=1e2*tau0; koalmin=min(koal);
q=0; nq=1e4; eps=1e-3; nit=length(koal); ra=1e1;
while ((q<nq) && (ra>eps))
    tau0c=(tau0a+tau0b)/2;
    tau0ma=0:tau0a/nit:tau0a; kta=polyfit([T0,T3,Tk],[0,tau0a/2,tau0a],2);
    tau0mb=0:tau0b/nit:tau0b; ktb=polyfit([T0,T3,Tk],[0,tau0b/2,tau0b],2);
    tau0mc=0:tau0c/nit:tau0c; ktc=polyfit([T0,T3,Tk],[0,tau0c/2,tau0c],2);
    fua=0;
    for k=1:length(tau0ma)
        fum=0;
        tema=polyval(kta,tau0ma(k));
        koala=opredKTPTKTochVer(koal, tem, tema);
        if (koala<=0)
            koala=koalmin;
        end
        for m=1:length(dv)
        fum(m)=1/(koala*alfs(m));
        end
        fua(k)=trapz(dv,fum);
    end
    fub=0;
    for k=1:length(tau0mb)
        fum=0; 
        temb=polyval(ktb,tau0mb(k));
        koalb=opredKTPTKTochVer(koal, tem, temb);
        if (koalb<=0)
            koalb=koalmin;
        end        
        for m=1:length(dv)
        fum(m)=1/(koalb*alfs(m));
        end
        fub(k)=trapz(dv,fum);
    end
    fuc=0;
    for k=1:length(tau0mc)
        fum=0;
        temc=polyval(ktc,tau0mc(k));
        koalc=opredKTPTKTochVer(koal, tem, temc);
        if (koalc<=0)
            koalc=koalmin;
        end        
        for m=1:length(dv)
        fum(m)=1/(koalc*alfs(m));
        end
        fuc(k)=trapz(dv,fum);
    end
    fa=trapz(tau0ma,fua)-tol;
    fb=trapz(tau0mb,fub)-tol;
    fc=trapz(tau0mc,fuc)-tol;
    la = VybCAB(fa,fb,fc,tau0a,tau0b,tau0c);
    tau0a=la(1); tau0b=la(2); 
    ra=abs(tau0b-tau0a);
    q=q+1
end
t=tau0c;
end