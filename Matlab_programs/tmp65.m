%Шамот, рассеяние только для шаров
function [ t ] = tmp65()
format long g;
npp=Kramers_n_Shamot_uk()-1;
npp=npp';
%npp=ZadMasPokPrel();
na=19.8; ko=27; de=1e-1;
dv=1.4:de:ko;
dvr=na:de:ko;
q=1; np=0;
for k=1:length(dv)
    if (dv(k)>dvr(q))
        np(q)=npp(k-1);
        q=q+1;
    end
end
fidTr=fopen('Shary_Tr-Sha2.txt','w'); fiddpp=fopen('Shary_dpp-Sha2.txt','w'); 
fidTrs=fopen('Shary_Trs-Sha2.txt','w'); fiddnra=fopen('Shary_dnra-Sha2.txt','w'); 
fidrazdT=fopen('Shary_razdT-Sha2.txt','w'); fidn1=fopen('Shary_nKBr-Sha2.txt','w'); 
fiddv=fopen('Shary_dl_vo-Sha2.txt','w'); 
np(q)=npp(length(npp)); 
no=1 %no - номер фракции
nom=2 %nom - номер концентрации
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
%t=2*PoisKorn();
t=2*erfcinv(1e-4)
ti=RasFra60(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
tm=ti;
end
function [ pk ] = PoisDn(vyfr,vyko,xv,vv,mo,si,mi,ma,dpp,dvr,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
TrRef=SredGrafRass(vyfr,vyko,xv); dlv=(1e6)*dlinavolny(); 
numele=length(dlv)-1; na=1; dvko=dvr
xv=xv'; xve=opredEffPar(xv,vyfr,vyko); 
n1=OprPokPreKBr(dvr)
vve=opredEffPar(vv,vyfr,vyko); dnre=0; q=1; eps=1e-2; chit=1e1; dnc=0; ide=1; TR=0; dpp=dpp'
while (dvr<=dvko)
    w=PoiskNomera(dlv,dvr); dlvo=dlv(w); TRs=TrRef(w-1)+(TrRef(w)-TrRef(w-1))*(dvr-dlv(w-1))/(dlv(w)-dlv(w-1)); TRs=exp(-TRs*xve)
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
function [ alsr ] = SredGrafRass(vyfr,vyko,xv)
%vyfr - выбор фракции: 1; vyko - выбор массовой доли: 1 - 0,7 %, 2 - 0,5 %
alsre=0;
if (vyfr == 1)
Tal=TrkoSham1(); 
Tal2=TrkoSham2(); 
Tal=privedkEdiPropus(Tal); 
Tal2=privedkEdiPropus(Tal2);
alsre=opredAlphSred2(Tal,Tal2,xv,vyko);
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
function [ al ] = opredAlphSred2(Tal,Tal2,xv,vyko)
nTal=length(Tal);
for k=1:nTal
    Ta(k)=-log(Tal(k))/xv(1); 
    Ta2(k)=-log(Tal2(k))/xv(2);
end
if (vyko == 1)
alsr=Ta;
elseif (vyko == 2)
alsr=Ta2;
end
al=alsr;
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
function ras60 = RasFra60(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[238.57,227.973]; 
mv=[1.706,1.1];
tol=[0.7,0.68];
rokbr=2.75; rov=2.69; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %в см3
vv(k)=mv(k)/(1e3*rov); %в см3
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %в мкм
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
mini=0; maxi=63; mo=(maxi+mini)/2; %МО
si=abs(mini-maxi)/di; %СКО
disp('Фракция менее 60 мкм');
depp=PoisDn(1,nom,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras60=depp; %ras60=alphaRas(tol,length(tol),Tr);
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
if (vyfr == 1)
        if (vyko == 1)
                no = 1;
        elseif (vyko == 2)
                no = 2;
        end
end
oep=mapa(no);
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
dl=(dlvoSham1()+dlvoSham2())/2; 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dv = dl;
end