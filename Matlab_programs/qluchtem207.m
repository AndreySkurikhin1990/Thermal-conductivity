function qlu = qluchten207()
%slu=4; 
nac=1; kon=4; fr=0;
format long g; 
y0=3e1*1e-3; te0=273.15;
dko=0.444905; %из статистического моделирования
dko3=1.3002913762; %из пористой структуры самого вермикулита
alfs=RasMasKoAbs();
alfs=dko*dko3*alfs;
npp=Kramers_n(); 
dl=RasshDiapDlinVoln();
temvh=arrTemHigh207(); temvh=temvh+te0; 
temvc=arrTemCold207(); temvc=temvc+te0; 
tepo=arrTepPot207(); 
F=13.85*1e-4;
qv=tepo/F; pt=length(temvc); 
tepv=koeftepv(qv,y0,temvh,temvc,pt);
temvs=(temvc+temvh)/2e0;
%Tz=1e2:5e1:12e2; Tz=Tz+te0;
delT=1e2; Tna=1e2; Tko=12e2;
Tz=Tna:delT:Tko; Tz=Tz+te0;
for slu=nac:kon
if (slu==1)
        %исходная фракция - установка Netzsch
        disp('Исходная фракция - установка Netzsch');
tem0=arrTem_VVF2(); tem0=tem0+te0;
ktp0=arrKTP_VVF2(); ketp0=polyfit(tem0,ktp0,2);
vyuv=0; vyfv=0; vysv=0; vmivmf=0; Te=Tz;
temvci=napMasEKTPVer(vyfv,vysv,vmivmf,Te,0,y0,vyuv);
temvhi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,1,y0,vyuv);
qoi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,2,y0,vyuv);
ektpvi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,3,y0,vyuv);
temvsi=napMasEKTPVer(vyfv,vysv,vmivmf,Te,4,y0,vyuv);
temvc=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,0);
temvh=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,1);
qo=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,2);
temvs=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,3);
ektpvi=vydelPol(temvci,temvhi,qoi,ektpvi,temvsi,4);
kc=polyfit(temvc,temvs,2);
kh=polyfit(temvh,temvs,2);
elseif (slu==2)
        %исходная фракция - стационарный метод
        %Tz=[347,617]+te0;
        %temvs
        disp('Исходная фракция - стационарный метод');
km=danPoTemTepl2071(temvs,tepv); kt1=km(1); kt2=km(2);
kc=danPoTemC2071(temvs,temvc); %ktc1=kc(1); ktc2=kc(2);
kh=danPoTemH2071(temvs,temvh); %kth1=kh(1); kth2=kh(2);
elseif (slu==3)
        %повторные измерения - стационарный метод
%Tz=[627,882]
%Tz=[356+352,619+600]/2+te0
disp('Повторные измерения - стационарный метод');
km=danPoTemTepl2072(temvs,tepv); kt1=km(1); kt2=km(2);
kc=danPoTemC2072(temvs,temvc); %ktc1=kc(1); ktc2=kc(2);
kh=danPoTemH2072(temvs,temvh); %kth1=kh(1); kth2=kh(2);
elseif (slu==4)
%Tz=[347+348,633+631]/2+te0;
%        Tz=[621,905];
        %после обжига при 1000 град С - стационарный метод
        disp('После обжига при 1000 °С - стационарный метод');
km=danPoTemTepl2073(temvs,tepv); kt1=km(1); kt2=km(2);
kc=danPoTemC2073(temvs,temvc); %ktc1=kc(1); ktc2=kc(2);
kh=danPoTemH2073(temvs,temvh); %kth1=kh(1); kth2=kh(2);
%Tz=[676,777,878,978,1079];
end
wmg=217e-3; wal=132e-3; wsi=368e-3;
p=length(Tz); koal=rasKoefAlpha(Tz,wmg,wsi,wal); 
%koeftep=0; ptpo=0; teco=0; teho=0; GrT=0; npT=0;
for k=1:p
    teco(k)=polyval(kc,Tz(k));
    teho(k)=polyval(kh,Tz(k));
    if (slu==1)
        koeftep(k)=abs(polyval(ketp0,Tz(k)));
    else
        koeftep(k)=kt1*Tz(k)+kt2;
    end
    GrT(k)=abs(teho(k)-teco(k))/y0;
    ptpo(k)=koeftep(k)*GrT(k);
    npT(k)=nsreddvVer(Tz(k));
    dko2=opreddko2(Tz(k));
    dko4=opredKTPTKTochVer(koal, Tz, Tz(k));
    alfs=alfs*dko2*dko4;
    alf(k)=sredRosSieg207(Tz(k),npp,alfs,dl);
    alfs=alfs/dko2/dko4;
end
alf=alf'
npT=npT'
GnpT(p)=0;
for k=2:p
    GnpT(k-1)=(npT(k)-npT(k-1))/(Tz(k)-Tz(k-1));
end
sigma=5.67e-8; %lamizl=0; lamte=0;
for k=1:p
    t=Tz(k)^3;
    t=abs(8*npT(k)*sigma*t/(3*alf(k)));
    t=t*(2*npT(k)+Tz(k)*GnpT(k));
    lamizl(k)=t;
    %lamte(k)=koeftep(k)-t;
end
n=length(Tz);
slu=slu'
lamizl=lamizl'
Tz=Tz'
koeftep=koeftep';
for k=1:length(lamizl)
    if (slu==nac)
        fr(k)=0;
    end
fr(k)=fr(k)+lamizl(k);
end
lamizl=lamizl';
Tz=Tz';
%t=zapvfile_alpha(alf);
%t=zapvfile_lamizl(lamizl);
%t=zapvfile_tem(Tz);
%-----
end
fr=fr'/abs(-nac+kon+1)
%slu=6
%дополнительные измерения, исходная фракция - стационарный метод
%Tz=sluchay6(npp,alfs,dl,Tz);
%p=sluchay6();
%if (slu<kon)
    %pl=plot(Tz(2:p),lamizl(2:p),'-b',Tz(2:p),lamte(2:p),'-k',Tz(2:p),koeftep(2:p),'-r');
    %legend('Излучательная компонента', 'Кондуктивная часть КТП', 'Общий КТП','location','best');
%else
    %pl=plot(Tz(2:p),lamizl(2:p),'-b',Tz(2:p),lamte(2:p),'-k',Tz(2:p),koeftep(2:p),'-r',Tz(2:p),lamef(2:p),'-g');
    %legend('Излучательная компонента', 'Кондуктивная часть КТП', 'Общий КТП 1','Общий КТП 2','location','best');
%end
%set(pl,'LineWidth',3); hold on; grid on; 
%xlabel({'Температура, К'}); 
%ylabel({'Составляющие коэффициента теплопроводности, Вт/(м*К)'}); 
%title({'График зависимости составляющих коэффициента теплопроводности от температуры'});
qlu=0;
end

%function sl = sluchay6(npp,alfs,dl,Tz)
function sl = sluchay6()
tepv=arrKTP_VVF1(); te0=273.15; y0=3e1*1e-3;
tepvs1=(tepv(1)+tepv(2)+tepv(3))/3;
tepvs2=tepv(4); tepvs3=(tepv(5)+tepv(6))/2;
temvh=0; temvh=arrTem1VVF1(); temvh=temvh+te0;
temvc=0; temvc=arrTem2VVF1(); temvc=temvc+te0;
temv3=0; temv3=arrTem3VVF1(); temv3=temv3+te0;
temvs=0; temvs=(1/2)*(temvc+temvh);
temvh1=(temvh(1)+temvh(2)+temvh(3))/3;
temvh2=temvh(4); temvh3=(temvh(5)+temvh(6))/2;
temvc1=(temvc(1)+temvc(2)+temvc(3))/3;
temvc2=temvc(4); temvc3=(temvc(5)+temvc(6))/2;
temvss1=(temvs(1)+temvs(2)+temvs(3))/3;
temvss2=temvs(4); temvss3=(temvs(5)+temvs(6))/2;
tem1=0; tem1=[temvss1 temvss2 temvss3];
teps1=0; teps1=[tepvs1 tepvs2 tepvs3];
km=polyfit(tem1,teps1,2);
tem2=0; tem2=[temvh1 temvh2 temvh3];
kh=polyfit(tem1,tem2,2);
tem2=0; tem2=[temvc1 temvc2 temvc3];
kc=polyfit(tem1,tem2,2); %Tz=0; 
Tz=[temvss3 temvss2 temvss1]; 
%Tz=tem1; 
p=length(Tz); koeftep=0; ptpo=0; teco=0; teho=0; GrT=0; npT=0;
for k=1:p
    teco(k)=polyval(kc,Tz(k));
    teho(k)=polyval(kh,Tz(k));
    koeftep(k)=abs(polyval(km,Tz(k)));
    %ptpo(k)=koeftep(k)*abs(teho(k)-teco(k))/y0;
    %GrT(k)=ptpo(k)/koeftep(k);
    %npT(k)=nsreddvVer(Tz(k));
    %alf(k)=sredRosSieg207(Tz(k),npp,alfs,dl);
end
%GnpT=0; 
%for k=2:p
    %GnpT(k)=(npT(k)-npT(k-1))/(Tz(k)-Tz(k-1));
%end
%sigma=5.67e-8; lamizl=0; lamte=0;
%for k=2:p
    %t=0;
    %t=abs(8*npT(k)*sigma*(Tz(k)^3)/(3*alf(k)));
    %t=t*(2*npT(k)+Tz(k)*GnpT(k));
    %lamizl(k)=t;
    %lamte(k)=koeftep(k)-t;
%end
Tz=Tz'
koeftep=koeftep'
sl=0;
end

function t = zapvfile_alpha(dl)
p=length(dl); 
fid = fopen('Koefficient_pogloscheniya_po_Rosselandu_2-07_mm.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',dl(k));
end
fclose(fid);
t=0;
end

function t = zapvfile_tem(dl)
p=length(dl); 
fid = fopen('Temperatura_KP_po_Rosselandu_2-07_mm.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',dl(k));
end
fclose(fid);
t=0;
end

function t = zapvfile_lamizl(dl)
p=length(dl); 
fid = fopen('Koefficient_teploprovodnosti_Luchistyy_po_Rosselandu_2-07_mm.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',dl(k));
end
fclose(fid);
t=0;
end

function rs = sredRosSieg207(tem,npp,alfs,dl)
te0=273.15; 
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4));
np2=nsredPlank2(dl,npp,tem);
Ibc=0; Ibz=0;
for k=1:length(npp)
c1m(k)=C1/(npp(k)^2);
c2m(k)=C2/npp(k);
t=(pi/2)*(c1m(k)*c2m(k))/sig;
dlv(k)=dl(k)/npp(k);
chi=exp(c2m(k)/tem/dlv(k));
zna=(chi-1)^2;
Ibz(k)=t*chi/zna/(dlv(k)^6)/(tem^5)/np2; 
Ibc(k)=Ibz(k)/alfs(k);
end
chasc=trapz(dlv,Ibc);
chasz=trapz(dlv,Ibz);
dlsvprfo2=chasc/chasz;
rs=1/dlsvprfo2;
end

function ns = nsredPlank2(dv,npp,tem)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
for k=1:length(npp)
    vl=c0;
    vl=vl/npp(k);
    c1=PP*(vl^2);
    c2=PP*vl/PB;
    %lambda=dv(k);
    lambda=dv(k)/npp(k);
    Ib(k)=2*pi*c1/((lambda^5)*(exp(c2/(lambda*tem))-1));
    %Ib(k)=Ib(k)/(npp(k)^2);    
Ibn(k)=(npp(k)^2)*Ib(k);
dvm(k)=lambda;
end
nc=trapz(dvm,Ibn);
nz=trapz(dvm,Ib);
ns=nc/nz;
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
x4=0; x3=0; x2=0; x=0; y=0; %de=0; de1=0; de2=0; de3=0; 	
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

function t = opreddko2(Ts)
dko2m=[4.68,4.69,5.65,13.17,20.2,27.81]/1e2; dko2t=6e2:2e2:16e2; 
for k=1:length(dko2m)
    dko2m(k)=1-dko2m(k);
end
t=opredKTPTKTochVer(dko2m, dko2t, Ts);
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
    if ((ktp<0) || (ktp>1))
        ktp=1;
    end
t=ktp;
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