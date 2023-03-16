%Определение уменьшения коэффициента поглощения
function tt = tmp70()
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
                %c=poiskUmenshAlphy(Te,vyfv,vysv,vmivmf,y0,vyuv,Te);
            end
        elseif (vyfv==1)
            for qvuv=1:length(nvyuv)
                vyuv=nvyuv(qvuv)
                c=poiskUmenshAlphy(Te,vyfv,vysv,vmivmf,y0,vyuv,Te);
            end
        end
    end
end
tt=c;
end

function tt = poiskUmenshAlphy(Te,vyfv,vysv,vmivmf,y0,vyuv,Tem0)
temvci=napMasEKTPVer(vyfv,vysv,vmivmf,Te,0,y0,vyuv); temvhi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,1,y0,vyuv); qoi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,2,y0,vyuv);
ektpvi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,3,y0,vyuv); temvsi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,4,y0,vyuv);
%---
temvc=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,0); temvh=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,1); qo=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,2);
temvs=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,3); ektpv=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,4); n=length(temvc);
dv=RasshDiapDlinVoln();
knuVer=RasMasKoAbs();
npp=Kramers_n();
%srpv207=[(0.0000165002+0.0000159712+0.0000158110)/3, 0.0000157067, 0.0000160877];
srpv207=[5e1,1e2,15e1]*1e-6;
srpv84=[(0.0000187133+0.0000185602)/2; (0.0000172837+0.0000186844)/2; (0.0000183523+0.0000181617)/2];
vy=1;
%vy=vysv+1;
if (vyfv==0)
srpe=srpv207(vy)
else
srpe=srpv84(vy)
end
q=1;
for w=1:n
    Tk=temvc(w);
    T0=temvh(w);
    %qob=qo(w);
    %ektp=ektpv(w);
    T3=temvs(w);
    Te=opredelTempRas(y0,T0,T3,Tk);
    ko=koefoslablen(Te,0);
    koal=koefoslablen(Te,1);
    ko2=koefoslablen(Te,2);
    p=opredAlpha(dv,npp,knuVer,ko,srpe,Te,koal,ko2);
    f=1;
    for k=1:length(p)
        if (p(k)<0)
            f=-1;
        end
    end
    if (f>0)
        t(q)=mean(Te);
        r(q)=mean(p);
        if (q>1)
        %r(q)=r(q)/r(1);
        end
        q=q+1;
    end
end
%r(1)=1;
Tem0=Tem0';
rr=vychKoe(t,Tem0,r)';
t=t'
r=r'
tt=0;
end

function [ f ] = opredAlpha(dv,npp,alf,ko,d,Te,koal,ko2)
ko=ko';
koal=koal';
alpsr=knusreddvvsrVer(Te,ko,1);
eps=epssr(Te,dv,npp,alf,koal,ko2);
refl=Rasronu(npp,dv,alf,Te);
a=1e-2; b=3e1;
for j=1:length(Te)
        alsr=alpsr(j);
        ab=eps(j);
        rb=refl(j);
        ad=alsr*d;
        ra=1e2; ep=1e-3; 
        k=0; kit=1e2; p=0;
    while ((k<kit) && (ra>ep))
        c=(a+b)/2;
        fa=exp(-ad*a);
        fda=fa;
        ffa=((rb*fa)^2);
        fa=fa/(1-ffa);
        fda=rb*(((1-rb)*fda)^2);
        fda=fda/(1-ffa);
        fa=(1-rb)*(1-fa)-ab-fda;
        fb=exp(-ad*b);
        fdb=fb;
        ffb=((rb*fb)^2);
        fb=fb/(1-ffb);
        fdb=rb*(((1-rb)*fdb)^2);
        fdb=fdb/(1-ffb);
        fb=(1-rb)*(1-fb)-ab-fdb;
        fc=exp(-ad*c);
        fdc=fc;
        ffc=((rb*fc)^2);
        fc=fc/(1-ffc);
        fdc=rb*(((1-rb)*fdc)^2);
        fdc=fdc/(1-ffc);
        fc=(1-rb)*(1-fc)-ab-fdc;
        if ((fc*fb>0) && (fa*fc<0)) 
        b=c; 
        end
        if ((fc*fa>0) && (fb*fc<0)) 
        a=c; 
        end
        if ((fc*fa>0) && (fb*fc>0)) 
            if (p>5)
                c=-1;
                break;
            end
        p=p+1;
        end
        ra=abs(b-a);
        k=k+1;
    end
    cm(j)=c;
end
Te=Te';
f=cm;
end

function [ rs ] = Rasronu(nsh,dl,knuSha,tem)
ronu=0; p=length(nsh);
for k=1:p
    hiSha(k)=abs(dl(k)*knuSha(k))/(4*pi)/nsh(k);
    nsh(k)=abs(nsh(k)); 
    nsha(k)=sqrt((nsh(k)^2+hiSha(k)^2)); 
    nsha(k)=abs(nsha(k));
    ronu(k)=5e-1+((nsha(k)-1)*(3*nsha(k)+1))/(6*(nsha(k)+1)^2);
    ronu(k)=ronu(k)-(2*(nsha(k)^3)*(nsha(k)^2+2*nsha(k)-1))/((nsha(k)^2+1)*(nsha(k)^4-1)); 
    podln=(nsha(k)-1)/(nsha(k)+1);
    podln=abs(podln);
    ronu(k)=ronu(k)+(8*(nsha(k)^4)*((nsha(k)^4)+1)*log(nsha(k)))/((nsha(k)^2+1)*((nsha(k)^4-1)^2));
    ronu(k)=ronu(k)+(nsha(k)^2)*((nsha(k)^2-1)^2)*log(podln)/((nsha(k)^2+1)^3);
end;
rs=ReflSreddv(ronu,nsha,dl,tem);
end

function [ knus ] = ReflSreddv(ronu,npp,dv,temmas)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c10=PP*c0^2;
c20=PP*c0/PB;
ronus=0;
for j=1:length(temmas)
c1=c10;
c2=c20;
tem=temmas(j);
for k=1:length(npp)
    c1=c1/(npp(k)^2);
    c2=c2/npp(k);
    lambda=dv(k)/npp(k);
    dl(k)=lambda;
    me=(exp(c2/lambda/tem)-1);
    Ib(k)=2*pi*c1/(lambda^5)/me;
Ibc(k)=ronu(k)*Ib(k);
c1=c10;
c2=c20;
end
nc=trapz(dl,Ibc);
nz=trapz(dl,Ib);
ronus(j)=nc/nz;
end
knus=ronus;
end

function [ tt ] = koefoslablen(Te,v)
wmg=217*1e-3; 
wal=132*1e-3; 
wsi=368*1e-3;
koal=rasKoefAlpha(Te,wmg,wsi,wal); 
dko=0.444905; 
dko3=1.3002914; 
for k=1:length(koal)
dko2(k)=opredeldko2(Te(k));
s(k)=dko*dko2(k)*dko3;
end
if (v==0)
tt=s;
elseif (v==1)
tt=koal;
elseif (v==2)
tt=dko2;
end
end

function [ tt ] = opredelTempRas(tol,T0,T3,Tk)
    kt=polyfit([0,tol/2,tol],[T0,T3,Tk],2);
    p=1;
    Nt=1e1-1e0;
    %y00=1e-3/p;
    y00=tol/p;
    y=0:(y00/Nt):y00;
    m=length(y);
    for k=1:m
    Te(k)=polyval(kt,y(k));
    end
    tt=Te;
end

function od2 = opredeldko2(Ts)
dko2m=[4.68,4.69,5.65,13.17,20.2,27.81]/1e2; 
dko2t=6e2:2e2:16e2; 
for k=1:length(dko2m)
    dko2m(k)=1-dko2m(k);
end
dko2=opredKTPTKTochVer(dko2m, dko2t, Ts);
od2=dko2;
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
    koo(k)=kumm(k)*wmg+kuma(k)*wal+kums(k)*wsi;
end
ko=koo'/wo;
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
x4=0; x3=0; x2=0; x=0; y=0; 
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
    else ktp=0; ko=0; 
    end
    if ((ktp<0) || (ktp>1))
        ktp=1;
    end
t=ktp;
end

%определяет ПТП, КТП, температуры для вермикулита
function [ vm ] = opredTempHolGorVer(ektp, ete, h0, v) 
n=length(ete); nit=1e6; hf=1e0; tnoscv=3e2;
ep=1e-2; d=1e-3; Thna=tnoscv; dt=1e0; 
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

function [ epsi ] = epssr(tem,dv,npp,alf,ko,ko2)
n=length(npp); p=length(tem);
for j=1:p
for k=1:n
hi=dv(k)*alf(k)/4/pi/npp(k);
np=sqrt(npp(k)^2+hi^2);
ep(k)=emdiel(np);
%ep(k)=emmet(npp(k),hi);
end
e=epssred(dv,ep,tem(j),npp);
eps(j)=e*ko(j)*ko2(j);
end
epsi=eps;
end

function es = epssred(dv,eps,T,npp)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*c0^2;
c2=PP*c0/PB;
iz0=0; iz1=0; dl=0;
c1t=c1; c2t=c2;
for j=1:length(dv)
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
iz0(j)=ProvAdekv(c1t/((la^5)*(exp(c2t/(la*T))-1)));
iz1(j)=eps(j)*iz0(j);
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
es=ProvAdekv(chi)/ProvAdekv(zna);
end

function [ als ] = knusreddvvsrVer(arrtem,ko,izm)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
npp=Kramers_n();
dv=RasshDiapDlinVoln();
knu=RasMasKoAbs();    
c1=PP*c0^2;
c2=PP*c0/PB;
koiz=1;
for j=1:length(arrtem)
    if (izm==1)
        koiz=ko(j);
    end
    tem=arrtem(j);
ct1=c1; ct2=c2;
for k=1:length(npp)
    ct1=ct1/(npp(k)^2);
    ct2=ct2/npp(k);
    lambda=dv(k)/npp(k);
    dv(k)=lambda;
    me=(exp(ct2/(lambda*tem))-1);
    Ib(k)=2*pi*ct1/((lambda^5)*me);
    knu(k)=knu(k)*koiz;
    Ibc(k)=knu(k)*Ib(k);
    knu(k)=knu(k)/koiz;
ct1=c1;
ct2=c2;
end
nc=trapz(dv,Ibc);
nz=trapz(dv,Ib);
als(j)=nc/nz;
end
knus=als;
end