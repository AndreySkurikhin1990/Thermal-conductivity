%Расчет коэффициента ослабления для ИТОМ и КВИ
function t = tmp73()
m=2;
ko=1e-6; koe=1e-3; koo=1e2; rokbr=2.75;
%mv=[0.95,1.19]*ko;
tol=[0.64,0.64,0.65,0.66]*koe;
n=length(tol);
mkbr=zeros(1,n);
hsl=zeros(1,n);
vkbr=zeros(1,n);
D=13e-3;
for k=1:n
mkbr(k)=(pi*(D^2)/4)*tol(k);
vkbr(k)=mkbr(k)/rokbr;
hsl(k)=8*k*ko;
end
for w=4:10 
rov=w*koo; j=w-3;
for k=1:n
vkvi(j,k)=vkbr(k)*hsl(k)/(tol(k)-hsl(k));
mv(j,k)=vkvi(j,k)*rov;
%k=RaschetDlyakvi(w,m);
%k=obratZadOprMasTol(w);
end
end
mv=mv/ko
%k=RaschetDlyaItom();
t=0;
end

function [ t ] = RaschetDlyakvi(ident,nom)
te0=273.15;
ko=1e-2;
roKao=2.6e3;
tem=te0+23; %комнатная температура
%вермикулит

knuVer=1e6*SredGraf();
dvVer=dlinyvolnVer();
nVer=Kramers_n_uk();
alsrVer=usrednen(tem,knuVer,dvVer,nVer);

dvVer=RasshDiapDlinVoln();
knuVer=RasMasKoAbs();
nVer=Kramers_n();
alsrVer=usrednen(tem,knuVer,dvVer,nVer);

%шамот
nSha=Kramers_n_Sha_uk(); 
dvSha=dlinyvolnSha();
knuSha=1e6*SredGrafSha();
alsrSha=usrednen(tem,knuSha,dvSha,nSha);

knuSha=preobMasSha_n();
nSha=Kramers_n_Sha();
dvSha=dvVer;
alsrSha=usrednen(tem,knuSha,dvSha,nSha);

%каолин
knuKao=1e6*SredGrafKao();
nKao=Kramers_n_Kao();
dvKao=dlinyvolnKao();
alsrKao=usrednen(tem,knuKao,dvKao,nKao);

nKao=Kramers_n_Kao_n();
knuKao=preobMasKao_n();
dvKao=dvVer;
alsrKao=usrednen(tem,knuKao,dvKao,nKao);

roSha=mean([2.33,2.23,2.13,2.14,2.12])*1e3
porSha=mean([14.06,17.10,21.28,19.26,26.54])*1e-2;
dvkvi=dvVer;
%--------------------------------------
%Гранулометрический состав фракции 0-10 мм
ropV=(95+1e2)/2;
roKgf=(120+125)/2;

vspV=[5.1,5.8,29.6,33.5,14,1e1,2,0]; 
n=length(vspV);
svspV=sum(vspV(1:n-1));
vspV=vspV(1:n-1)*1e2/svspV*ko; vspV(n)=0;
svspV=sum(vspV);
vspV=vspV/svspV;

vsKgf=[1.1,3.8,4.7,9.8,10.6,26.4,37.6,6];
n=length(vsKgf);
svsKgf=sum(vsKgf(1:n-1));
vsKgf=vsKgf(1:n-1)*1e2/svsKgf*ko; vsKgf(n)=0;
svsKgf=sum(vsKgf);
vsKgf=vsKgf/svsKgf;

vspV2=[0,0,0,(2e1+4e1)/2, (25+5e1)/2, (25+45)/2,0,0];
vsKgf2=[0,0,0,(3e1+45)/2, (3e1+4e1)/2, (1e1+25)/2,0,0];
eps=1e-6;
for k=1:length(vsKgf2)
    if (vsKgf2(k)<eps)
        vsKgf2(k)=1e0/2e0; 
        vspV2(k)=1e0/2e0;
    end
end
svsKgf2=sum(vsKgf2);
vsKgf2=vsKgf2*1e2/svsKgf2*ko;
svsKgf2=sum(vsKgf2);
svspV2=sum(vspV2);
vspV2=vspV2*1e2/svspV2*ko;
svspV2=sum(vspV2);
vsKgf2=vsKgf2'/svsKgf2;
vspV2=vspV2'/svspV2;

roVer010=0;
for k=1:length(vsKgf2)
roVer010=roVer010+ropV*vspV(k)*vspV2(k)+roKgf*vsKgf(k)*vsKgf2(k);
end
roVer010=roVer010';
sVer0612=[(2e1+4e1)/2,(3e1+45)/2];
ssVer0612=sum(sVer0612);
sVer0612=sVer0612*1e2/ssVer0612;
ssVer0612=sum(sVer0612);
roVer0612=ropV*sVer0612(1)/ssVer0612+roKgf*sVer0612(2)/ssVer0612; %определение плотности

saloSha=mean([42.6,38.5,34.9,33.6,37.6])*ko;
smgoSha=mean([0.18,0.32,0.37,0.37,0.3])*ko;
ssioSha=mean([53.3,56.3,59.5,61.1,57])*ko;
sfeoSha2=mean([1.34,1.45,1.86,1.73,1.56])*ko;
stioSha=mean([1.17,2.04,1.95,1.98,2.09])*ko;

soSha=[saloSha,smgoSha,ssioSha];
ssoSha=sum(soSha);
soSha=soSha*1e2/ssoSha*ko;
ssoSha=sum(soSha);
soSha=soSha/ssoSha;

wkaol=mean([85,90])*ko;
shoKao=13.96*wkaol*ko;
saloKao=mean([30,45,33,41,39.5*wkaol])*ko;
ssioKao=mean([42,45,46.54*wkaol])*ko;
sfeoKao2=mean([1.2,1.3,0.4,5,1.2,1.5])*ko;
stioKao=mean([1.2,1.4])*ko;
smgoKao=mean([0.1,0.2])*ko;
scaoKao=mean([0.3,0.5,0.9])*ko;
snaokoKao=mean([0,0.015,0.55,2])*ko;

soKao=[saloKao,smgoKao,ssioKao];
ssoKao=sum(soKao);
soKao=soKao*1e2/ssoKao*ko;
ssoKao=sum(soKao);
soKao=soKao/ssoKao;

ssioVer=36.8*ko;
saloVer=13.2*ko;
smgoVer=21.7*ko;
sfeoVer2=8.3*ko; sfeoVer=1.2*ko; scaoVer=2*ko;
soVer=[saloVer,smgoVer,ssioVer]';
ssoVer=sum(soVer);
soVer=soVer*1e2/ssoVer*ko;
ssoVer=sum(soVer);
soVer=soVer/ssoVer;

if (nom==1)
    roVer=roVer010;
else
    roVer=roVer0612;
end

%КВИ
obsove=1.8;
%--------------------------------------
switch (ident)
    case (4) %400
disp('kvi-400');
vKao=(2e2+22e1)/2; 
vSha=0;
npp=Kramers_n_kvi400();
    case (5) %500
disp('kvi-500');
if (nom==1) %0-10 мм
vKao=(2e2+275)/2;
vSha=(1e2+21e1)/2;
else %0,6-1,2 мм
vKao=(2e2+23e1)/2; 
vSha=(1e2+15e1)/2; 
end
npp=Kramers_n_kvi500();
    case (6) %600
if (nom==1) %0-10 мм
vKao=(25e1+375)/2; 
vSha=(1e2+28e1)/2; 
else %0,6-1,2 мм
vKao=(23e1+32e1)/2; 
vSha=(15e1+2e2)/2; 
end
npp=Kramers_n_kvi600();
    case (7) %700
if (nom==1) %0-10 мм
vKao=(275+475)/2; 
vSha=(1e2+38e1)/2;
else %0,6-1,2 мм
vKao=(23e1+425)/2;
vSha=(1e2+35e1)/2;
end
npp=Kramers_n_kvi700();
    case (8) %800
if (nom==1) %0-10 мм
vKao=(37e1+575)/2;
vSha=(1e2+42e1)/2;
else %0,6-1,2 мм
vKao=(32e1+52e1)/2;
vSha=(1e2+375)/2; 
end
npp=Kramers_n_kvi800();
    case (9) %900
if (nom==1) %0-10 мм
vKao=(475+675)/2; 
vSha=(1e2+39e1)/2;
else %0,6-1,2 мм
vKao=(4e2+6e2)/2;
vSha=(1e2+34e1)/2;
end
npp=Kramers_n_kvi900();
    case (10) %1000
if (nom==1) %0-10 мм
vKao=(54e1+7e2)/2; 
vSha=(1e2+46e1)/2; 
else %0,6-1,2 мм
vKao=(48e1+7e2)/2;
vSha=(1e2+325)/2;
end
npp=Kramers_n_kvi1000();
end
vSha=vSha/roSha;
vKao=vKao/roKao;
%-------------------
Sp=RasshDiapDlinVoln();
knukvi=oprknukvi(knuVer, knuKao, knuSha, nVer, nKao, nSha, tem, dvkvi, obsove, vSha/roSha, vKao/roKao);
sokvi=oprsodoxkvi(soVer,soSha,soKao,obsove,roVer,vSha,roSha,vKao,roKao)';
disp(ident);
n2=opredn2(npp,tem,Sp)
srRoskvi=sredRosSiegel(tem,npp,Sp, knukvi, n2)
q=ZapisFileOptiokvi(knukvi,ident);
q=ZapisFileOptiokvi_dv(dvkvi);
t=0;
end

function  [ kkvi ] = oprknukvi(knuVer, knuKao, knuSha, nVer, nKao, nSha, tem, dvkvi, vVer, vSha, vKao)
pVer=vVer/(vVer+vKao+vSha);
pSha=vSha/(vVer+vKao+vSha);
pKao=vKao/(vVer+vKao+vSha);
knukvi=knuVer*pVer+knuSha*pSha+knuKao*pKao;
nkvi=nVer*pVer+nSha*pSha+nKao*pKao;
alsrkvi=usrednen(tem,knukvi,dvkvi,nkvi);
knukvi=knuVer*pVer+knuSha*(pSha+pKao);
alsrkvi=usrednen(tem,knukvi,dvkvi,nkvi)
kkvi=knukvi;
end

function [ r ] = oprsodoxkvi(soVer,soSha,soKao,vVer,roVer,vSha,roSha,vKao,roKao)
saloVer=soVer(1); smgoVer=soVer(2); ssioVer=soVer(3);
saloSha=soSha(1); smgoSha=soSha(2); ssioSha=soSha(3);
saloKao=soKao(1); smgoKao=soKao(2); ssioKao=soKao(3);
mSha=roSha*vSha;
mVer=roVer*vVer;
mKao=roKao*vKao;
mob=mSha+mVer+mKao;
mdSha=mSha/mob;
mdKao=mKao/mob;
mdVer=mVer/mob;
salokvi=saloSha*mdSha+saloVer*mdVer+saloKao*mdKao; 
smgokvi=smgoSha*mdSha+smgoVer*mdVer+smgoKao*mdKao; 
ssiokvi=ssioSha*mdSha+ssioVer*mdVer+ssioKao*mdKao; 
r=[salokvi,smgokvi,ssiokvi];
end

function t = ZapisFileOptiokvi_dv(massi)
p=length(massi);
fid = fopen('DlinaVolny_kvi.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptiokvi(massi,f)
p=length(massi);
switch (f)
    case (4)
        fid = fopen('Koefficient_pogloscheniya_kvi400.txt','w');
    case (5)
        fid = fopen('Koefficient_pogloscheniya_kvi500.txt','w');
    case (6)
        fid = fopen('Koefficient_pogloscheniya_kvi600.txt','w');
    case (7)
        fid = fopen('Koefficient_pogloscheniya_kvi700.txt','w');
    case (8)
        fid = fopen('Koefficient_pogloscheniya_kvi800.txt','w');
    case (9)
        fid = fopen('Koefficient_pogloscheniya_kvi900.txt','w');
    case (10)
        fid = fopen('Koefficient_pogloscheniya_kvi1000.txt','w');
end
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function [ t ] = RaschetDlyaItom()
te0=273.15;
roKao=2.6*1e3;
tem=te0+23; %комнатная температура
nitom=Kramers_n_Itom();
dvitom=dlinyvolnitom();
knuitom=SredGrafItom();
disp('Itom-440');
alsritom440=usrednen(tem,knuitom,dvitom,nitom)
n2=opredn2(nitom,tem,dvitom)
alsritomRos440=sredRosSiegel(tem,nitom,dvitom, knuitom, n2)

dvVer=RasshDiapDlinVoln();
q=ZapisFileOptio_dv(dvVer);
%dvVer=dlinyvolnVer();
knuVer=RasMasKoAbs();
%q=postgraf(dvVer*1e6,knuVer);

%knuVer=1e6*SredGrafKoAbsVer();
%knuVer=1e6*SredGraf();
nVer=Kramers_n();
alsrVer=usrednen(tem,knuVer,dvVer,nVer);

nSha=Kramers_n_Sha_uk(); 
dvSha=dlinyvolnSha();
knuSha=1e6*SredGrafSha();
alsrSha=usrednen(tem,knuSha,dvSha,nSha);

roVer=250;
roSha=2.39174*1e3*0.91;
sVer=6e1*1e-2;
sSha=1-sVer;
pV=sVer/roVer;
pS=sSha/roSha;
pVer=pV/(pV+pS);
pSha=pS/(pV+pS);
alsVS=alsrVer*pVer+alsrSha*pSha;

alsrKao=1e6*SredGrafKao();
nKao=Kramers_n_Kao();
dvKao=dlinyvolnKao();
alsrKao=usrednen(tem,alsrKao,dvKao,nKao);

roVer=2.5*1e2;
roSha=2.6*1e3;
sVer=6e1*1e-2;
sSha=1-sVer;
pV=sVer/roVer;
pS=sSha/roSha;
pVer=pV/(pV+pS);
pKao=pS/(pV+pS);

disp('Sham itom440');
alsVK=alsrVer*pVer+alsrSha*pKao;
disp('Kao itom-440');
alsVK=alsrVer*pVer+alsrKao*pKao;

pV=5e1*1e-2;
pK=1-pV;
pV=pV/roVer;
pK=pK/roSha;
pVer=pV/(pV+pK);
pKao=pK/(pV+pK);
pSha=0;

disp('Itom-620');
knuSha=preobMasSha_n();
%knuKao=preobMasKao_n();
%disp('Kao itom-620');
%alitom620=pVer*knuVer+pSha*knuSha+pKao*knuKao;
nitom620=Kramers_n_Itom620();
%alitom620sred=usrednen(tem,alitom620,dvVer,nitom620);
%q=ZapisFileOptio(alitom620,1);
%disp('Sham itom-620');
Sp=RasshDiapDlinVoln();
alitom620=pVer*knuVer+(pSha+pKao)*knuSha;
q=ZapisFileOptio(alitom620,1);
alitom620sred=usrednen(tem,alitom620,dvVer,nitom620)
n2=opredn2(nitom620,tem,Sp)
alsritomRos620=sredRosSiegel(tem,nitom620,Sp, alitom620, n2)

disp('itom-860');
pV=25*1e-2;
pS=3e1*1e-2;
pK=(1-pV-pS);
pK=pK/roKao;
pV=pV/roVer;
pS=pS/roSha;
pVer=pV/(pV+pS+pK);
pSha=pS/(pV+pS+pK);
pKao=pK/(pV+pS+pK);
%disp('Kao itom-860');
%alitom860=pVer*knuVer+pSha*knuSha+pKao*knuKao;
nitom860=Kramers_n_Itom860();
%alitom860sred=usrednen(tem,alitom860,dvVer,nitom860);
%q=ZapisFileOptio(alitom860,2);
%disp('Sham itom-860');
alitom860=pVer*knuVer+(pSha+pKao)*knuSha;
q=ZapisFileOptio(alitom860,2);
alitom860sred=usrednen(tem,alitom860,dvVer,nitom860)
n2=opredn2(nitom860,tem,Sp)
alsritomRos860=sredRosSiegel(tem,nitom860,Sp, alitom860, n2)

disp('itom-1000');
pV=2e1*1e-2;
pS=4e1*1e-2;
pK=(1-pV-pS);
pK=pK/roKao;
pV=pV/roVer;
pS=pS/roSha;
pVer=pV/(pV+pS+pK);
pSha=pS/(pV+pS+pK);
pKao=pK/(pV+pS+pK);
%disp('Kao itom-1000');
%alitom1000=pVer*knuVer+pSha*knuSha+pKao*knuKao;
nitom1000=Kramers_n_Itom1000();
%alitom1000sred=usrednen(tem,alitom1000,dvVer,nitom1000);
%q=ZapisFileOptio(alitom1000,3);
%disp('Sham itom-1000');
alitom1000=pVer*knuVer+(pSha+pKao)*knuSha;
q=ZapisFileOptio(alitom1000,3);
alitom1000sred=usrednen(tem,alitom1000,dvVer,nitom1000)
n2=opredn2(nitom1000,tem,Sp)
alsritomRos1000=sredRosSiegel(tem,nitom1000,Sp, alitom1000, n2)

pVK=(alsritom440-alsVK)*1e2/alsVK;
pSK=(alsritom440-alsVS)*1e2/alsVS

t=0;
end

function t = ZapisFileOptio_dv(massi)
p=length(massi);
fid = fopen('DlinaVolny_itomN.txt','w');
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function t = ZapisFileOptio(massi,f)
p=length(massi);
switch (f)
    case (1)
        fid = fopen('Koefficient_pogloscheniya_itom620.txt','w');
    case (2)
        fid = fopen('Koefficient_pogloscheniya_itom860.txt','w');
    case (3)
        fid = fopen('Koefficient_pogloscheniya_itom1000.txt','w');
end
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function nk = opredn2(PokPrel,Tsr,dv)
mn2=0;
for k=1:length(PokPrel)
mn2(k)=PokPrel(k)^2;
end
nk=usrednen(Tsr,mn2,dv,PokPrel);
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

function [ dl ] = dlinyvolnitom()
Sp=(dlvoItom1()+dlvoItom2())/2;
leSp=length(Sp);
for k=1:leSp
	Sp(k)=1e-2/Sp(k);
end
dl=Sp;
end

function [ dl ] = dlinyvolnSha()
Sp=(dlvoSham1()+dlvoSham2())/2;
leSp=length(Sp);
for k=1:leSp
	Sp(k)=1e-2/Sp(k);
end
dl=Sp;
end

function [ dl ] = dlinyvolnKao()
Sp=(dlvokao1()+dlvokao2())/2;
leSp=length(Sp);
for k=1:leSp
	Sp(k)=1e-2/Sp(k);
end
dl=Sp;
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

function [ dlv ] = dlinyvolnVer()
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end

function [ pre ] = preobMasKao_n()
format longg;
alp2=1e6*SredGrafKao();
dlv2=(dlvokao1()+dlvokao2())/2;
ledlv2=length(dlv2);
dlv1=RasshDiapDlinVoln();
ledlv1=length(dlv1);
dlv1=preobDlinVoln(ledlv1,ledlv2,dlv1,dlv2,1);
%dlv1=dlvoVer53101(); ledlv1=length(dlv1);
t=1;
for j=1:ledlv1
    a=dlv1(j);
    fa=-1;
for k=2:ledlv2
        if ((dlv2(k-1)>=a) && (dlv2(k)<a))
            a1a=dlv2(k-1);
            b1a=dlv2(k);
            f1a=alp2(k-1);
            f2a=alp2(k);
            fa=1;
            break;
        end;
end;
if (fa>0)
    tm(t)=f1a+(a-a1a)*(f2a-f1a)/(b1a-a1a);
    dvtm(t)=a;
    t=t+1;
end
end;
a1=dlv2(1);            
b1=dlv2(2);
f1=alp2(1);
f2=alp2(2);
tm(1)=f1+(f2-f1)*(dlv1(1)-a1)/(b1-a1);
n=ledlv1;
m=ledlv2;
a1=dlv2(m-1);            
b1=dlv2(m);
f1=alp2(m-1);
f2=alp2(m);
tm(n)=f2+(f2-f1)*(dlv1(n)-b1)/(b1-a1);
tm(n-1)=(tm(n)+tm(n-2))/2;
for k=1:n
    tm(k)=ProvAd(tm(k));
end
pre=tm;
end

function [ w ] = preobDlinVoln(ledlv1,ledlv2,dlv1,dlv2,v)
for k=1:ledlv1
    dlv1(k)=1e-2/dlv1(k);
end
for k=1:ledlv2
    dlv2(k)=1e-2/dlv2(k);
end
if (v==1)
    w=dlv1;
else
    w=dlv2;
end
end

function f = postgraf(koox,tetaxx)
pl=plot(koox,tetaxx,'-b');
%legend('КО(ДВ)','location','best');
set(pl,'LineWidth',3);
hold on;
grid on;
xlabel({'ДВ, м'});
ylabel({'КО, м-1'});
title({'График зависимости КП от ДВ'});
f=0;
end

function [ rs ] = sredRosSiegel(tem,npp,dl, alfs, n2)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
C1=PP*(c0^2);
C2=PP*c0/PB;
eps=1e-6;
%pi=3.1415926535897932;
sig=2*C1*(pi^5)/(15*(C2^4));
for k=1:length(npp)
c1m=C1/(npp(k)^2);
c2m=C2/npp(k);
Consk=(pi/2)*(c1m*c2m)/(sig*n2);
la=dl(k)/npp(k);
dlv(k)=la;
chi=exp(c2m/la/tem);
zna=(chi-1)^2;
Ibz(k)=Consk*chi/zna/(la^6)/(tem^5); 
if (abs(alfs(k))>eps)
Ibc(k)=Ibz(k)/alfs(k);
end
end
chasc=trapz(dlv,Ibc);
chasz=trapz(dlv,Ibz);
dlsvprfo2=chasc/chasz;
dlsvprfo2=1/dlsvprfo2;
rs=dlsvprfo2';
end

function t = obratZadOprMasTol(vkvi)
c=num2str(vkvi);
s=strcat('Koefficient_pogloscheniya_kvi',c,'00.txt');
fid=fopen(s,'r');
kpskvi=fscanf(fid,'%f');
sv='DlinaVolny_kvi.txt';
fidv=fopen(sv,'r');
dv=fscanf(fidv,'%f');
te0=273.15; tem=22; tem=tem+te0;
switch (vkvi)
case (4)
npp=Kramers_n_kvi400();
case (5)
npp=Kramers_n_kvi500();
case (6)
npp=Kramers_n_kvi600();
case (7)
npp=Kramers_n_kvi700();
case (8)
npp=Kramers_n_kvi800();
case (9)
npp=Kramers_n_kvi900();
case (10)
npp=Kramers_n_kvi1000();
end
al=usrednen(tem,kpskvi,dv,npp);
t=al;
end