%рассеяние только для шаров
function [ t ] = tmp42()
format long g;
npp=Kramers_ver()-1; npp=npp';
%npp=ZadMasPokPrel();
na=13.1; ko=27; de=1e-1;
dv=1.4:de:ko;
dvr=na:de:ko;
q=1; np=0;
for k=1:length(dv)
    if (dv(k)>dvr(q))
        np(q)=npp(k-1);
        q=q+1;
    end
end
fidTr=fopen('Shary_Tr.txt','w'); fiddpp=fopen('Shary_dpp.txt','w'); 
fidTrs=fopen('Shary_Trs.txt','w'); fiddnra=fopen('Shary_dnra.txt','w'); 
fidrazdT=fopen('Shary_razdT.txt','w'); fidn1=fopen('Shary_nKBr.txt','w'); 
fiddv=fopen('Shary_dl_vo.txt','w'); 
np(q)=npp(length(npp)); 
no=1 %no - номер фракции
nom=3 %nom - номер концентрации
%q=VyvodTr(dvr,no,nom);
for k=1:length(dvr)
    n1(k)=OprPokPreKBr(dvr(k));
end
n1=n1'; n1=0; tr=0; np=np'; npma=max(np); npmi=min(np);
for k=1:length(np)
trk=nacha(np(k),dvr(k),no,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
tr(k)=trk(1);
end
fclose(fidTr); fclose(fiddpp); fclose(fiddnra); fclose(fidrazdT); fclose(fidTrs); fclose(fidn1); fclose(fiddv);
t=tr;
end
function [ tm ] = nacha(popr,dv,no,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
format long g; %определяет коэффициент рассеяния для таблетки KBr + вермикулит
%t=2*PoisKorn()
t=2*erfcinv(1e-4)
switch (no) 
    case 1 %60 мкм
ti=RasFra60(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
    case 2 %60-100 мкм
ti=RasFra60_100(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
    case 3 %100-150 мкм
ti=RasFra100_150(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
    case 4 %150-200 мкм
ti=RasFra150_200(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
end
tm=ti;
end
function [ pk ] = PoisDn(vyfr,vyko,xv,vv,mo,si,mi,ma,dpp,dvr,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
TrRef=SredGrafRass(vyfr,vyko); dlv=(1e6)*dlinavolny; numele=length(dlv)-1; na=1; dvko=dvr
xv=xv'; xve=opredEffPar(xv,vyfr,vyko); TR=1e0;
n1=OprPokPreKBr(dvr)
vve=opredEffPar(vv,vyfr,vyko); dnre=0; q=1; eps=1e-2; chit=1e1; dnc=0; dpp=dpp'
while (dvr<=dvko)
    w=PoiskNomera(dlv,dvr); dlvo=dlv(w); TRs=TrRef(w-1)+(TrRef(w)-TrRef(w-1))*(dvr-dlv(w-1))/(dlv(w)-dlv(w-1))
    Tra=TRs; Trb=1e0; sc=1; razdT=abs(Tra-Trb); TR=(Tra+Trb)/2; dnc=0;
while ((razdT>eps) && (sc<chit))
    TR=(Tra+Trb)/2;
    dna=npmi; dnb=npma; dnc=(dna+dnb)/2; ep=1e-5; Nit=1e2; k=0; ra=abs(dna-dnb);
while ((ra>ep) && (k<Nit))
    dnc=(dna+dnb)/2;
    if (dna==0)
        fa=1;
    else
    fa=opreKoefProp(vve,mo,si,mi,ma,dlvo,dna);
    end
    fa=fa-TR;
    if (dnb==0)
        fa=1;
    else
    fb=opreKoefProp(vve,mo,si,mi,ma,dlvo,dnb);
    end
    fb=fb-TR;
    if (dnc==0)
        fc=1;
    else
    fc=opreKoefProp(vve,mo,si,mi,ma,dlvo,dnc);
    end
    fc=fc-TR;
    ro=vybPoKo(fa,fb,fc,dna,dnb,dnc); 
    dna=ro(1); dnb=ro(2); dnc=ro(3); 
    k=k+1;
    ra=abs(dna-dnb);
end
dnra=dnc
sc=sc+1;
if (dnra>dpp)
    Trb=TR;
else
    Tra=TR;
end
razdT=abs(Tra-Trb)
end
TR=TR'
dvr=dvr'
%w=w'
dvr=dvr+1e-1;
dnre(q)=dnc; q=q+1;
end
fprintf(fidTr,'%0.15f\n',TR);
fprintf(fiddpp,'%0.15f\n',dpp);
fprintf(fidTrs,'%0.15f\n',TRs);
fprintf(fiddnra,'%0.15f\n',dnre(1));
fprintf(fidrazdT,'%0.15f\n',razdT);
fprintf(fidn1,'%0.15f\n',n1);
fprintf(fiddv,'%0.15f\n',dvko);
pk=dnre;
end
%Данная функция должна возвращать вероятность попадания при выстрелах, mo - матем. ожид-е, si - СКО, mi, ma - мин. и макс. разм., no - выбор (цил. или прям. пар-д)
function mdn = opreKoefProp(vvk,mo,si,mi,ma,dlvo,dn)
Nr=1e0; raot=13e3; blis=0; ocp=0; sepoo=0; rosr=0; %raot - размер отверстия, blis - число удачных попаданий, ocp - общее число попыток
for g=1:Nr %усреднение по реализациям
suv=0; k=1; nk=1e9;
while ((suv<vvk) && (k<nk))
tetani=pi*rand(); %направление излучения - угол между направлением точки падения и осью Z (вверх)
phini=pi*rand(); %угол между проекцией точки падения на горизонтальную плоскость и осью X
ocp=ocp+1; %общее число выстрелов
ds=mo+si*randn(); ds=prov(ds,mi,ma); rs=ds/2; %диаметр шара
nalu=PoiskDeltaUgPad(dn,dlvo,phini,tetani,rs); deug=nalu(1); ko=nalu(2); deug=proverka(deug); ko=proverka(ko); rosr=rosr+ko;
vys=PopNePop(dlvo,deug,raot); suv=suv+4*pi*(rs^3)/3; blis=blis+vys; sepoo=sepoo+pi*(rs^2);
k=k+1;
end
end
w=blis/ocp; rosr=rosr/ocp;
p=sepoo/Nr/(pi*(raot^2)/4);
mdn=p*w*(1-rosr)+(1-p); %пропускание за счет рассеяния
end
function [ ug ] = oprUgl(phi,teta)
phk=1/cos(teta)/tan(phi); phk=atan(phk);
nalu=[0     -1      0]; nalu=PreobOtIskhKoorKNov(nalu,teta,phi,phk);
phe=asin(nalu(1));
if (phe<0)
    phe=abs(phe);
end
phe=provUgla(phe);
phet=acos(nalu(2));
if (phet>pi/2)
    phet=pi-phet;
end
phet=provUgla(phet);
phet=(phet+phe)/2;
ug=[phk    phet];
end
function [ ug ] = oprUglNaVykh(phi,teta,nalu)
phk=1/cos(teta)/tan(phi); phk=atan(phk);
nalu=PreobOtIskhKoorKNov(nalu,teta,phi,phk);
phe=asin(nalu(1));
if (phe<0)
    phe=abs(phe);
end
phe=provUgla(phe);
phet=acos(nalu(2));
if (phet>pi/2)
    phet=pi-phet;
end
phet=provUgla(phet);
phet=(phet+phe)/2;
ug=[phk    phet];
end
function pu = provUgla(phi)
if (phi<(-pi))
    phi=phi+pi;
end
if (phi<(-pi/2))
    phi=-phi-pi;
end
if (phi>pi)
    phi=phi-pi;
end
if (phi>(pi/2))
    phi=pi-phi;
end
if (phi<0)
    phi=abs(phi);
end
pu=phi;
end
%в плоскости ksis, etas
function pr10 = PraCha10(phi,phis,n1,n2)
pr10=n1*sin(phi)-n2*sin(phis);
end
function pr11 = PraCha11(phis,phiss,n1,n2)
pr11=n1*sin(phiss)-n2*sin(phis); %dn/n1
end
%показатель преломления бромида калия
function opn = OprPokPreKBr(la)
lam=[1 2 10 20 25]; pp=[1.487 1.48 1.52 1.48 1.453]; ko=polyfit(lam,pp,length(pp)-1); opn=polyval(ko,la);
end
function [ pesdp ] = PoiskPhiShtr(dnr,phieta,n1,n2)
phipre=0; pvo=0; a=0; b=pi/2; ep=1e-6; Nit=1e2; k=0; ra=abs(a-b); dephi=0; %находим угол преломления
if (dnr<0)
    phipre=asin(n2/n1); phipre=provUgla(phipre);
end
if (((phieta<phipre) && (dnr<0)) || (dnr>0))
while ((ra>ep) && (k<Nit))
    c=(a+b)/2;
    fa=PraCha10(phieta,a,n1,n2);
    fb=PraCha10(phieta,b,n1,n2);
    fc=PraCha10(phieta,c,n1,n2);
    ro=vybPoKo(fa,fb,fc,a,b,c); 
    a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
    phietas=provUgla(c);
else
    pvo=1; phietas=pi/2; dephi=2*pi;
end
pesdp=[phietas  dephi   pvo     phipre];
end
%поиск изменения угла после преломления
function rop = PoiskReflPerv(phieta,phietas,n1,n2)
rpa=(n2*cos(phieta)-n1*cos(phietas))/(n2*cos(phieta)+n1*cos(phietas)); %коэффициент отражения параллельный плоскости падения
rpe=(n1*cos(phieta)-n2*cos(phietas))/(n1*cos(phieta)+n2*cos(phietas)); %коэффициент отражения перпендикулярный к плоскости падения
rop=(abs(rpa^2)+abs(rpe^2))/2;
end
function rov = PoiskReflVtor(ugalp,phietass,n1,n2)
rpe=(n2*cos(ugalp)-n1*cos(phietass))/(n2*cos(ugalp)+n1*cos(phietass)); %коэффициент отражения перпендикулярный к плоскости падения 
rpa=(n1*cos(ugalp)-n2*cos(phietass))/(n1*cos(ugalp)+n2*cos(phietass)); %коэффициент отражения параллельный плоскости падения
rov=(abs(rpe^2)+abs(rpa^2))/2;
end
function [ ug ] = oprUglVykh(nalu,tovy,phi,teta,ugalp,dvnaplu,n1,n2)
phipre=pi/2; 
if (n2>n1)
    phipre=asin(n1/n2);
end
if (ugalp<=phipre)
    a=0; b=pi/2; Nit=1e2; k=0; ra=abs(a-b); ep=1e-6; %находим угол выхода из шара
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; fa=PraCha11(ugalp,a,n1,n2); %n1 - ПП KBr, n2 - ПП ВВ
    fb=PraCha11(ugalp,b,n1,n2); 
    fc=PraCha11(ugalp,c,n1,n2);
    ro=vybPoKo(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
phietass=c; phietass=provUgla(phietass);
nopl=NormKPloskVto(nalu,tovy); nopl=nopl'; tovy=tovy'; nalu=nalu'; dsc=(-nopl*tovy-nopl*nalu)/2;
pke=oprUglNaVykh(phi,teta,nalu); phiksi=pke(1); phie=pke(2); ugalp=(ugalp+phie)/2;
nalunk=PreobOtIskhKoorKNov(nalu,teta,phi,phiksi);
ksiss=sign(nalunk(1))*sin(phietass); etass=sign(nalunk(2))*cos(phietass); dzetass=0; 
naplu=[ksiss etass dzetass]; %новые оси ksi, eta, dzeta
naluskss=PreobOtKsiEtaDzetaKIskhKoor(naplu,teta,phi,phiksi); dvnaluskss=ModulVectora(naluskss);
dephi=acos((naplu*naluskss')/dvnaplu/dvnaluskss); dephi=proverka(dephi);
rov=PoiskReflVtor(ugalp,phietass,n1,n2);
else
    dephi=2*pi; rov=1;
end
ug=[dephi    rov];
end
function [ dpkov ] = PoiskUgVykhIzShara(n1,n2,nalusk,tovy,rs,dvnalusk)
dephi=0; rov=0; z=tovy(3); teta=acos(z/rs); x=tovy(1); 
phix=acos(x/rs/sin(teta)); y=tovy(2); phiy=asin(y/rs/sin(teta));
if ((imag(phix)==0) && (imag(phiy)==0) && (abs(y)<=rs) && (abs(x)<=rs))
    if (y<=0)
        phix=-phix+2*pi;
        if (x<=0)
        phiy=-phiy+pi;
        else
        phiy=phiy+2*pi;
        end
    else
        phix=phix+2*pi;
        if (x>=0)
        phiy=phiy+2*pi;
        else
        phiy=-phiy+3*pi;
        end
    end
    ugalp=acos((nalusk*tovy')/dvnalusk/rs);
    phi=(phix+phiy)/2;
    ugvy=oprUglVykh(nalusk,tovy,phi,teta,ugalp,dvnalusk,n1,n2); 
else
    dephi=pi; rov=1;
end
dpkov=[dephi     rov];
end
function [ pk ] = PoiskDeltaUgPad(dnr,dv,phi,teta,rs)
ep=1e-5; n1=OprPokPreKBr(dv); n2=n1+dnr;
ugphike=oprUgl(phi,teta); phiksi=ugphike(1); %угол поворота в плоскости раздела
phieta=ugphike(2); %угол падения
pesdp=PoiskPhiShtr(dnr,phieta,n1,n2);
phietas=pesdp(1); dephi=pesdp(2); pvo=pesdp(3); phpr=pesdp(4);
if (pvo<1)
rop=PoiskReflPerv(phieta,phietas,n1,n2);
nvsn=NormKPlosk(phi,teta,phiksi,phietas); nvs=PreobOtKsiEtaDzetaKIskhKoor(nvsn,teta,phi,phiksi); dnvs=ModulVectora(nvs); %нормаль к плоскости падения
dzetas=0; ksis=-sin(phietas); etas=-cos(phietas); nalus=[ksis   etas    dzetas]; %направление распространения луча в сечении после преломления в новых штрихованных координатах ksis, etas, dzetas
nalusks=PreobOtKsiEtaDzetaKIskhKoor(nalus,teta,phi,phiksi); dvnalu=ModulVectora(nalusks); 
topa=[rs*sin(teta)*cos(phi)     rs*sin(teta)*sin(phi)   rs*cos(teta)]; kto=[0 0 0]; %точка падения в старых координатах
dsc=-nvs*topa'; roo=abs((kto*nvs')+dsc)/dnvs; kck=-roo*nvs; %расстояние до плоскости
ron=abs((kck*nvs')+dsc); %если направление другое,
if (ron>ep)
    kck=-kck; %то меняются координаты центра круга
end %проверка, что центр окружности лежит в плоскости падения
rksp=-kck+topa;
drksp=ModulVectora(rksp); %радиус круга в плоскости сечения
nalu=[0     -1      0]; nalun=PreobOtIskhKoorKNov(nalu,teta,phi,phiksi);
alug=acos((-rksp*nalusks')/dvnalu/drksp);
if (alug>pi/2)
    alug=pi-alug;
end
tovy=2*cos(alug)*drksp*nalusks+topa;
dpvko=PoiskUgVykhIzShara(n1,n2,nalusks,tovy,rs,dvnalu); dephi=dpvko(1); rov=dpvko(2);
else
    ksisp=-sin(phieta); etasp=-cos(phieta);
    nalusp=[ksisp etasp 0]; naluskp=PreobOtKsiEtaDzetaKIskhKoor(nalusp,teta,phi,phiksi); dvnaluskp=ModulVectora(naluskp);
    ksiso=-sin(phieta); etaso=cos(phieta);
    naluso=[ksiso etaso 0]; nalusko=PreobOtKsiEtaDzetaKIskhKoor(naluso,teta,phi,phiksi); dvnalusko=ModulVectora(nalusko);
    dephi=acos(naluskp*nalusko')/dvnalusko/dvnaluskp;
    rop=1; rov=1;
end
pk=[dephi     rop*rov];
end
function [ nv ] = NormKPloskVto(nalu,tovy)
x(1)=nalu(1); x(2)=tovy(1); 
y(1)=nalu(2); y(2)=tovy(2); 
z(1)=nalu(3); z(2)=tovy(3);
C=1; b=-C*z;
a(1,1)=x(1); a(1,2)=y(1); a(2,1)=x(2); a(2,2)=y(2); a=inv(a);
nvs=a*b'; nvs(3)=C;
nv=nvs;
end
function [ nv ] = NormKPlosk(phi,teta,phiksi,phietas)
nvsk1=[0  -1  0]; nv1=PreobOtIskhKoorKNov(nvsk1,teta,phi,phiksi); %падающий луч - волновое число
nv2=[-sin(phietas)  -cos(phietas)   0]; 
a(1,1)=nv1(1); a(1,2)=nv1(2); 
a(2,1)=nv2(1); a(2,2)=nv2(2); 
C=1; b(1)=-C*nv1(3); b(2)=-C*nv2(3); %преломленный луч
a=inv(a); nvs=a*b'; nvs(3)=C;
nv=nvs;
end
%------------------пока не трогать--------------
function pd = PopNePop(lam,deltatetass,raot)
p=0; deko=lam/raot;
if (deltatetass>deko)
    p=0;
else p=1;
end
pd=p;
end
function [ v ] = vybPoKo(fa,fb,fc,a,b,c)
if ((fc*fb)>0) 
        if ((fa*fc)<0) 
            b=c; 
        end
end
    if ((fc*fa)>0) 
        if ((fb*fc)<0) 
            a=c; 
        end
    end
    v=[a b c];
end
function [ alsr ] = SredGrafRass(vyfr,vyko)
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
%vyfr - выбор фракции: 1 - 60, 2 - 60-100, 3 - 100-150, 4 - 150-200; vyko - выбор массовой доли: 1 - 0,1 %, 2 - 0,2 %, 3 - 0,3 %
switch (vyfr)
    case 1
Tal=privedkEdiPropus(TrkoVer5311()); Tal2=privedkEdiPropus(TrkoVer5312());Tal3=privedkEdiPropus(TrkoVer5313()); 
Tal4=privedkEdiPropus(TrkoVer53101()); Tal5=privedkEdiPropus(TrkoVer53102()); Tal6=privedkEdiPropus(TrkoVer53103()); 
alsre=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,vyko);
    case 4
Tal=privedkEdiPropus(TrkoVer5681()); Tal2=privedkEdiPropus(TrkoVer5682()); Tal3=privedkEdiPropus(TrkoVer5683()); 
Tal4=privedkEdiPropus(TrkoVer5691()); Tal5=privedkEdiPropus(TrkoVer5692()); Tal6=privedkEdiPropus(TrkoVer5693());
alsre=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,vyko);
    case 2
Tal=privedkEdiPropus(TrkoVer5421()); Tal2=privedkEdiPropus(TrkoVer5422()); Tal3=privedkEdiPropus(TrkoVer5423()); 
Tal4=privedkEdiPropus(TrkoVer5431()); Tal5=privedkEdiPropus(TrkoVer5432()); Tal6=privedkEdiPropus(TrkoVer5433());
Tal7=privedkEdiPropus(TrkoVer5441()); Tal8=privedkEdiPropus(TrkoVer5442()); Tal9=privedkEdiPropus(TrkoVer5443());
alsre=opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,vyko);
    case 3
Tal=privedkEdiPropus(TrkoVer5551()); Tal2=privedkEdiPropus(TrkoVer5552()); Tal3=privedkEdiPropus(TrkoVer5553()); 
Tal4=privedkEdiPropus(TrkoVer5561()); Tal5=privedkEdiPropus(TrkoVer5562()); Tal6=privedkEdiPropus(TrkoVer5563());
Tal7=privedkEdiPropus(TrkoVer5571()); Tal8=privedkEdiPropus(TrkoVer5572()); Tal9=privedkEdiPropus(TrkoVer5573());
alsre=opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,vyko);
end
alsr=alsre;
end
function [ prked ] = privedkEdiPropus(ar)
arr=ar;
p=length(arr);
for k=1:p
    if (arr(k)>1)
        arr(k)=1;
    end
    if (arr(k)<0)
        arr(k)=0;
    end
end
prked=arr;
end
function [ al ] = opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,vyko)
switch (vyko)
    case 1
alsr=opredAlphSred2(Tal,Tal4);
    case 2
alsr=opredAlphSred2(Tal2,Tal5);
    case 3
alsr=opredAlphSred2(Tal3,Tal6);
end
al=alsr;
end
function [ al ] =opredAlphSred2(Tal,Tal2)
alsr=0; nTal=length(Tal);
for k=1:nTal 
    alsr(k)=Tal(k)+Tal2(k);
end
al=alsr/2;
end
function [ al ] = opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,vyko)
switch (vyko)
    case 1
alsr=opredAlphSred3(Tal,Tal4,Tal7);
    case 2
alsr=opredAlphSred3(Tal2,Tal5,Tal8);
    case 3
alsr=opredAlphSred3(Tal3,Tal6,Tal9);
end
al=alsr;
end
function [ al ] =opredAlphSred3(Tal,Tal2,Tal3)
alsr=0; nTal=length(Tal);
for k=1:nTal 
    alsr(k)=Tal(k)+Tal2(k)+Tal3(k);
end
al=alsr/3;
end
function pk = PoisKorn()
ro=0; a=1e-5; b=1e5; ep=1e-6; Nit=3e1; k=0; ffv=1e-4; ra=abs(a-b); %ffv = 0,01 % - такая доля частиц не попадает в диапазон фракции
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; fa=erfc(a)-ffv; fb=erfc(b)-ffv; fc=erfc(c)-ffv; 
    ro=vybPoKo(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3); k=k+1;
    ra=abs(a-b);
end
pk=c;
end
function [ xv ] = EffectTols(mkbr,mv,tol,rov)
rokbr=2.75; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %в см3
vv(k)=mv(k)/(1e3*rov); %в см3
xvo(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %в мкм
end
xv=xvo;
end
function [ ras60 ] = RasFra60(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[250.239 249.740 249.223 250.395 250.336 249.55]; 
mv=[0.283 0.464 0.812 0.22 0.547 0.777]; 
tol=[0.73 0.72 0.72 0.72 0.71 0.7];
rokbr=2.75; rov=0.49; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %в см3
vv(k)=mv(k)/(1e3*rov); %в см3
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %в мкм
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
xv=EffectTols(mkbr,mv,tol,rov);
mini=0; maxi=6e1; mo=(maxi+mini)/2; %МО
si=abs((mini-maxi)/di); %СКО
disp('Фракция менее 60 мкм');
depp=PoisDn(1,nom,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras60=depp; %ras60=alphaRas(tol,length(tol),Tr);
end %di - для поиска СКО
function [ ras60_100 ] = RasFra60_100(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[250 249.629 249.294 249.706 249.510 249.307 250.328 249.604 249.206]; 
mv=[0.255 0.539 0.809 0.295 0.517 0.756 0.36 0.534 0.843]; 
tol=[0.72 0.71 0.7 0.7 0.73 0.72 0.74 0.7 0.76]; 
rokbr=2.75; rov=0.52; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
xv=EffectTols(mkbr,mv,tol,rov);
mini=6e1; maxi=1e2; mo=(maxi+mini)/2; %МО
si=abs((mini-maxi)/di); %СКО
disp('Фракция 60-100 мкм');
depp=PoisDn(2,nom,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras60_100=depp; %ras60_100=alphaRas(tol,length(tol),Tr);
end
function [ ras100_150 ] = RasFra100_150(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[249.913 249.607 249.218 249.929 249.695 249.306 250.405 249.625 249.348];
mv=[0.315 0.473 0.709 0.293 0.528 0.83 0.27 0.493 0.764]; 
tol=[0.74 0.74 0.72 0.72 0.71 0.7 0.78 0.73 0.76]; 
rokbr=2.75; rov=0.53; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
xv=EffectTols(mkbr,mv,tol,rov);
mini=1e2; maxi=15e1; mo=(maxi+mini)/2; %МО
si=abs((mini-maxi)/di); %СКО
disp('Фракция 100-150 мкм');
depp=PoisDn(3,nom,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras100_150=depp; %ras100_150=alphaRas(tol,length(tol),Tr);
end
function [ ras150_200 ] = RasFra150_200(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[250.882 249.590 249.213 250.299 249.441 249.365];
mv=[0.320 0.533 0.849 0.223 0.502 0.797]; 
tol=[0.76 0.72 0.69 0.73 0.73 0.73]; 
rokbr=2.75; rov=0.56; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); vv(k)=mv(k)/(1e3*rov); xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vv(k)=vv(k)*((1e4)^3); vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
end
xv=EffectTols(mkbr,mv,tol,rov);
mini=15e1; maxi=2e2; mo=(maxi+mini)/2; %МО
si=abs((mini-maxi)/di); %СКО
disp('Фракция 150-200 мкм');
depp=PoisDn(4,nom,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras150_200=depp; %ras150_200=alphaRas(tol,length(tol),Tr);
end
function p = prov(hk,mi,ma)
hp=hk;
if (hk<mi)
    hp=mi;
end
    if (hk>ma)
        hp=ma;
    end
    p=hp;
end
function oep = opredEffPar(mapa,vyfr,vyko)
switch (vyfr)
    case {1,4}
        ne=2; 
        switch (vyko)
            case 1
                no = [1 4];
            case 2
                no = [2 5];
            case 3
                no = [3 6];
        end
    case {2,3}
        ne=3;
        switch (vyko)
            case 1
                no = [1 4 7];
            case 2
                no = [2 5 8];
            case 3
                no = [3 6 9];
        end
end
s=0;
for k=1:ne
    s=s+mapa(no(k));
end
oep=s/ne;
end
function [ pik ] = PreobOtKsiEtaDzetaKIskhKoor(nalus,teta,phi,phiksi) %преобразование от новых координат к старым для шара
ksis=nalus(1); etas=nalus(2); dzetas=nalus(3); 
ksi=ksis*cos(phiksi)+dzetas*sin(phiksi); eta=etas; dzeta=-ksis*sin(phiksi)+dzetas*cos(phiksi);
xss=eta; yss=-dzeta; zss=-ksi;
xs=xss*sin(teta)-zss*cos(teta); zs=xss*cos(teta)+zss*sin(teta); ys=yss;
x=xs*cos(phi)-ys*sin(phi); y=xs*sin(phi)+ys*cos(phi); z=zs;
pik=[x  y   z];
end
function [ pik ] = PreobOtIskhKoorKNov(nalusk,teta,phi,phiksi) %преобразование от старых координат к новым для шара
x=nalusk(1); y=nalusk(2); z=nalusk(3);
xs=x*cos(phi)+y*sin(phi); ys=-x*sin(phi)+y*cos(phi); zs=z;
xss=xs*sin(teta)+zs*cos(teta); zss=-xs*cos(teta)+zs*sin(teta); yss=ys;
ksi=-zss; eta=xss; dzeta=-yss;
ksis=ksi*cos(phiksi)-dzeta*sin(phiksi); etas=eta; dzetas=ksi*sin(phiksi)+dzeta*cos(phiksi);
pik=[ksis      etas      dzetas];
end
function t = proverka(d)
tt=d;
if (isnan(d))
    tt=0;
end
if (isinf(d))
    tt=0;
end
t=tt;
end
function pn = PoiskNomera(dlvo,dvr)
n=length(dlvo); f=1; q=1;
for k=1:n
    if ((dlvo(k)>dvr) && (f>0))
        f=0; 
        q=k;
        break;
    end
end
pn=q;
end
function t = ModulVectora(ar)
s=0;
for k=1:length(ar)
    x=ar(k); s=s+(x^2);
end
t=sqrt(s);
end
function t = VyvodTr(dvr,nofr,noko)
dlv=(1e6)*dlinavolny();
TrRef=SredGrafRass(nofr,noko); 
for k=1:length(dvr)
    dvro=dvr(k);
    w=PoiskNomera(dlv,dvro); 
    TRs(k)=TrRef(w-1)+(TrRef(w)-TrRef(w-1))*(dvro-dlv(w-1))/(dlv(w)-dlv(w-1));
end
TRs=TRs'
t=0;
end
function [ dv ] = dlinavolny()
dl=0; dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dv = dl;
end
function [ npp ] = ZadMasPokPrel()
npp=[1.00011907361118
           1.0000610020059
          1.00010352435283
          1.00026544052075
          1.00046634658064
          1.00071079319126
          1.00084336308156
          1.00112273029087
          1.00146884236497
          1.00187606258167
          1.00234647526234
          1.00297509678711
           1.0038447561407
          1.00557561497049
          1.00524543264575
          1.00344825286348
          1.00028495907096
         0.998694221040409
          0.99767219693943
         0.997640928504681
         0.998126555209527
         0.999232025583488
          1.00020455247694
          1.00098820778839
          1.00173684191243
          1.00241634776438
           1.0029878596913
          1.00349978226967
          1.00391664257249
          1.00435682449077
          1.00470433098372
          1.00509531099981
          1.00555177014779
          1.00604223669635
          1.00652072987217
          1.00699457205259
           1.0075267886669
          1.00811060294143
          1.00880016332392
          1.00936827968911
          1.01010092854399
          1.01081180723323
          1.01164698873235
           1.0124997425883
          1.01352557321169
          1.01330258736592
          1.01218778806115
          1.01037276042359
          1.00898820254082
            1.010327920541
          1.01202057702625
          1.01221859110109
          1.01337306978747
          1.01451886643304
          1.01473130696709
          1.01600779373315
          1.01648482881737
          1.01716832550177
          1.01795032778285
          1.01872235309146
           1.0195149235251
          1.02030697918631
          1.02114964005253
          1.02200053284077
          1.02281424451782
          1.02361007180792
           1.0245399767928
          1.02536375752239
           1.0261534526837
          1.02717218584581
          1.02869209178493
          1.03013340681228
            1.031276930309
          1.03248900585036
          1.03363082182131
          1.03477593512202
          1.03657298409728
          1.03845632544771
           1.0378619001816
          1.03717667701603
          1.03652587078846
          1.03623793358092
           1.0362198068916
          1.03507477510147
          1.03236898100473
          1.02778566795049
          1.02088502932505
          1.01295042080023
          1.00591847888212
          1.00000042563761
         0.994926151953635
         0.991222020043967
         0.989126371130242
         0.988741030085959
         0.989188715347181
         0.989730582399168
         0.989857684750855
          0.99030614888852
         0.992263147821438
         0.994772072698691
         0.996829228393934
         0.998482362605873
           1.0002417873933
          1.00214864479417
          1.00470874912981
          1.00729797151487
          1.00971436653592
          1.01116022209667
           1.0121685106475
           1.0130303462452
          1.01448794553559
          1.01631085235256
          1.01830292416993
          1.02001775821042
          1.02127473987166
          1.02243811804975
          1.02330371076368
           1.0237619525922
          1.02415289792584
           1.0244964336966
          1.02477136340842
          1.02503510110572
          1.02535914832008
          1.02562349018686
          1.02613891876655
          1.02697314261601
          1.02794809240384
          1.02878013159331
          1.02975007809159
          1.03037635058781
          1.03048516134328
          1.02994479796484
          1.02922541520964
          1.02862804314824
          1.02799597625443
          1.02744005256444
          1.02659480564644
           1.0252840272806
          1.02427795075048
          1.02362557851852
          1.02335599258304
          1.02359220235741
          1.02420448976143
          1.02460592378073
          1.02516870534284
          1.02611206166743
          1.02698123957535
          1.02771612274801
          1.02863091565657
          1.02946224286162
          1.03070677562327
          1.03190371055802
          1.03293399699302
          1.03393759025396
          1.03521926689801
          1.03653389149038
          1.03789071961467
          1.03944936735797
          1.04127694107874
          1.04300286914829
          1.04477442081304
          1.04680538665475
          1.04878614375839
          1.05066800730621
          1.05253103201956
          1.05443936958158
          1.05639463591819
          1.05791592316764
           1.0599348230194
          1.06159182270925
          1.06299507498948
          1.06454948556952
          1.06557899995986
          1.06653602072092
          1.06737794718449
          1.06756981592448
            1.067563507831
          1.06659946360297
          1.06580400989599
          1.06532026889177
          1.06475133919382
          1.06514065302929
          1.06550263165024
          1.06639247221319
          1.06702386504619
          1.06799392673869
          1.06829184785964
          1.06905469842653
          1.06906297765933
          1.06925770487327
          1.06912281527762
          1.06918880636191
          1.07012635337116
           1.0706567704854
          1.07150069530132
           1.0719896291137
          1.07265757409144
          1.07251679524807
          1.07212667957055
          1.07000627614216
           1.0682133114926
          1.06557169653535
          1.06212447911007
          1.05788482317465
          1.05410207943269
          1.04968527167777
          1.04619995568765
          1.04201554820026
          1.03830147305285
          1.03500575691345
           1.0311239210541
          1.02810554234663
          1.02493636847548
          1.02101166923469
           1.0165016909263
          1.01330897165649
          1.00954951215232
          1.00767509852304
          1.00591131584428
          1.00439054853671
          1.00441327196132
          1.00351007039367
          1.00222744696439
          1.00137583480815
          1.00244648954098
          1.00200691668322
          1.00137178371663
           1.0023396259477
          1.00327661542804
          1.00241271382437
          1.00279698176169
           1.0031460010226
          1.00253063952328
          1.00326916224715
          1.00296229802362
          1.00295434223511
          1.00357012728321
          1.00575302283094
          1.00538882563976
          1.00344410692782
          1.00367175716039
          1.00595623880419
          1.00696347972079
          1.00592286292602
          1.00726223519338
          1.00858107813206
          1.01033547360824
          1.00839608085399
          1.00789787542337
           1.0115317230167
          1.01153608202299
          1.01216302446949
          1.01478883025654
          1.01333448449347
          1.01119062148113
          1.01222124010081
          1.01147091404249]-1;
end