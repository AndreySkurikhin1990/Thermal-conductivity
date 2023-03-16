%рассеяние только для ПП
function [ t ] = tmp43()
npp=Kramers_ver()-1; 
%npp=PokazPrelMas();
na=1.4; ko=27; de=1e-1; ep=de/1e1;
dv=1.4:de:ko;
dvr=na:de:ko;
q=1; np=0;
for k=1:length(dv)
    n1(k)=OprPokPreKBr(dv(k));
    if (abs(dv(k)-dvr(q))<ep)
        np(q)=npp(k);
        q=q+1;
    end
end
n1=n1'
fidTr=fopen('Pryamougolnye_Parallelepipedy_Tr.txt','w'); fiddpp=fopen('Pryamougolnye_Parallelepipedy_dpp.txt','w'); 
fidTrs=fopen('Pryamougolnye_Parallelepipedy_Trs.txt','w'); fiddnra=fopen('Pryamougolnye_Parallelepipedy_dnra.txt','w'); 
fidrazdT=fopen('Pryamougolnye_Parallelepipedy_razdT.txt','w'); fidn1=fopen('Pryamougolnye_Parallelepipedy_nKBr.txt','w'); 
fiddv=fopen('Pryamougolnye_Parallelepipedy_dl_vo.txt','w'); 
npma=max(np); npmi=min(np); 
for k=1:length(np)
trk=nacha(np(k),dvr(k),npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
tr(k)=trk(1);
end
fclose(fidTr); fclose(fiddpp); fclose(fiddnra); 
fclose(fidrazdT); fclose(fidTrs); fclose(fidn1); fclose(fiddv);
t=0;
end
function tm = nacha(popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
format long g; %определяет коэффициент рассеяния для таблетки KBr + вермикулит
%t=2*PoisKorn();
t=2*erfcinv(1e-4)
no=3 %no - номер фракции
nom=1 %nom - номер концентрации
switch (no) 
    case 1
ti=RasFra60(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
    case 2
ti=RasFra60_100(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
    case 3
ti=RasFra100_150(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
    case 4
ti=RasFra150_200(t,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
end
tm=ti;
end
function pk = PoisDn(vyfr,vyko,xv,vv,mo,si,mi,ma,dpp,dvr,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
TrRef=SredGrafRass(vyfr,vyko,xv); dlv=(1e6)*dlinavolny(); 
dvko=dvr; xv=xv'
xve=opredEffPar(xv,vyfr,vyko); ide=1; n1=OprPokPreKBr(dvr); TR=0; dnre=0;
vve=opredEffPar(vv,vyfr,vyko); dnre=0; q=1; eps=1e-2; dpp=dpp'
while (dvr<=dvko)
    w=PoiskNomera(dlv,dvr); dlvo=dlv(w); chit=1e0; TRs=TrRef(w-1)+(TrRef(w)-TrRef(w-1))*(dvr-dlv(w-1))/(dlv(w)-dlv(w-1))
    Tra=TRs; Trb=1e0; sc=0; razdT=abs(Tra-Trb); TR=(Tra+Trb)/2; dnc=1e2;
while ((razdT>eps) && (sc<chit))
    TR=(Tra+Trb)/2;
    dnb=npma; dna=npmi; %dnb=1e3; 
    if (dna>0)
        dna=-dna;
    end
    if (abs(dna)>1)
        dna=-1;
    end
    if (abs(dnb)>1e2)
        dnb=1e2;
    end
    dna=dna';
    dnb=dnb';
    ep=1e-5; Nit=1e2; k=0; ra=abs(dna-dnb); dnc=(dna+dnb)/2; ddnc=dnc; kk=0; nikk=1e1;
while ((ra>ep) && (k<Nit))
    dnc=(dna+dnb)/2;
    if (dna==0)
        fa=1;
    else
    fa=opreKoefProp(vve,mo,si,mi,ma,dlvo,dna,ide);
    end
    fa=fa-TR;
    if (dnb==0)
        fa=1;
    else
    fb=opreKoefProp(vve,mo,si,mi,ma,dlvo,dnb,ide+1);
    end
    fb=fb-TR;
    if (dnc==0)
        fc=1;
    else
    fc=opreKoefProp(vve,mo,si,mi,ma,dlvo,dnc,ide+2);
    end
    fc=fc-TR;
    ro=vybPoKo(fa,fb,fc,dna,dnb,dnc); 
    dna=ro(1); dnb=ro(2); dnc=ro(3); 
    k=k+1;
    ra=abs(dna-dnb);
    if (abs(ddnc-dnc)<ep)
        kk=kk+1;
        if (kk>nikk)
        break;
        end
    end
    ddnc=dnc;
    if (dnc>dpp)
    dnb=dnc;
    else
    dna=dnc;
    end
end
dnra=dnc
sc=sc+1;
if (dnra>dpp)
    Trb=TR; 
else
    Tra=TR; 
end
razdT=abs(Tra-Trb);
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
pk=dnre(1);
end
function [ nk ] = PovorotIzStarkNov(vec,upvoz,teta,phi)
vec=PovorotVokrugOZ(vec,upvoz); 
vec=PovorotnaUgolTeta(vec,teta);
vec=PovorotnaUgolPhi(vec,phi);
nk=vec;
end
function [ nk ] = PovorotVokrugOZ(koord,upvoz)
x=koord(1); y=koord(2); z=koord(3);
xs=x*cos(upvoz)+y*sin(upvoz); ys=-x*sin(upvoz)+y*cos(upvoz); zs=z;
nk=[xs,ys,zs];
end
function [ nk ] = PovorotnaUgolPhi(koord,phi)
x=koord(1); y=koord(2); z=koord(3);
xs=x; ys=y*cos(phi)-z*sin(phi); zs=y*sin(phi)+z*cos(phi);
nk=[xs,ys,zs];
end
function [ nk ] = PovorotnaUgolTeta(koord,teta)
x=koord(1); y=koord(2); z=koord(3);
xs=x*sin(teta)-z*cos(teta); ys=y; zs=x*cos(teta)+z*sin(teta);
nk=[xs,ys,zs];
end
%Данная функция должна возвращать вероятность попадания при выстрелах
%mo - матем. ожид-е, si - СКО, mi, ma - мин. и макс. разм., no - выбор (цил. или прям. пар-д)
function mdn = opreKoefProp(vvk,mo,si,mi,ma,dlvo,dn,ide)
n1=OprPokPreKBr(dlvo); n2=n1+dn; rosr=0;
Nr=1e0; raot=13e3; blis=0; ocp=0; sepoo=0; %raot - размер отверстия в микронах, blis - число удачных попаданий
for g=1:Nr %усреднение по реализациям
suv=0; k=0; nk=1e12;
while ((suv<vvk) && (k<=nk))
hk=mo+si*randn(); hk=prov(hk,mi,ma); %определение длины направляющей
ak=mo+si*randn(); ak=prov(ak,mi,ma); bk=mo+si*randn(); bk=prov(bk,mi,ma); %определение длин основания для ПП
tetani=pi*rand(); %направление излучения - угол между плоскостью основания и осью Z (вверх)
phini=pi*rand(); %угол между проекцией оси цилиндра на горизонтальную ось и осью X - коллинеарно направлению излучения
upvoz=pi*rand(); %угол поворота вокруг оси Z
suv=suv+hk*ak*bk; %объем ПП
    vat=[ak,0,0]; vat=PovorotIzStarkNov(vat,upvoz,tetani,phini); dva=DlinaVectora(vat); %одно основание на вертикальной плоскости
    vbt=[0,bk,0]; vbt=PovorotIzStarkNov(vbt,upvoz,tetani,phini); dvb=DlinaVectora(vbt); %второе основание на вертикальной плоскости
    vabt=vat+vbt; dvab=DlinaVectora(vabt); %большая диагональ на вертикальной плоскости
    vab=vat-vbt; dvabm=DlinaVectora(vab); %меньшая диагональ на вертикальной плоскости
    vht=[0,0,-hk]; vht=PovorotIzStarkNov(vht,upvoz,tetani,phini); dvh=DlinaVectora(vht); %вектор - проекция образующей боковой поверхности на вертикальную плоскость, три поворота векторов a=(a,0,0), b=(0,b,0), h=(0,0,-h), XZ - в новых координатах вертикальная плоскость    
    vaht=vat+vht; vbht=vbt+vht; vabht=vat+vbt+vht; uab=acos(abs(vbt*vat'))/dva/dvb; %острый угол между проекциями векторов a и b
    uah=acos(abs(vat*vht'))/dva/dvh; %острый угол между проекциями a и h
    ubh=acos(abs(vbt*vht'))/dvb/dvh; %острый угол между проекциями b и h
    svo=dvb*dva*sin(uab); svb=dva*dvh*sin(uah)+dvb*dvh*sin(ubh); svb=svb+svo; %проекции основания и боковой поверхности
    mz=VysLuchaPP(svo,svb,tetani,hk,abs(ak),abs(bk),n1,n2,dn,dlvo,raot,upvoz,phini,dva,dvb,dvh,uab,uah,ubh,ide,vat,vbt,vabt,vht); 
    vys=mz(1); vys=proverka(vys); blis=blis+vys; %1 - успешно, 0 - нет
    ko=mz(2); ko=proverka(ko); rosr=rosr+ko; sepoo=sepoo+hk*(ak*cos(upvoz)+bk*sin(upvoz)); ocp=ocp+1;
k=k+1;
end
end
w=blis/ocp;
rosr=rosr/ocp;
p=sepoo/Nr/(pi*(raot^2)/4);
mdn=p*w*(1-rosr)+(1-p); %пропускание за счет рассеяния
end
function rop = PoiskReflPerv(phi,phis,n1,n2)
rpa=(n2*cos(phi)-n1*cos(phis))/(n2*cos(phi)+n1*cos(phis)); %коэффициент отражения параллельный плоскости падения
rpe=(n1*cos(phi)-n2*cos(phis))/(n1*cos(phi)+n2*cos(phis)); %коэффициент отражения перпендикулярный к плоскости падения
rop=(abs(rpa^2)+abs(rpe^2))/2;
end
function rov = PoiskReflVtor(phips,phiss,n1,n2)
rpe=(n2*cos(phips)-n1*cos(phiss))/(n2*cos(phips)+n1*cos(phiss)); %коэффициент отражения перпендикулярный к плоскости падения 
rpa=(n1*cos(phips)-n2*cos(phiss))/(n1*cos(phips)+n2*cos(phiss)); %коэффициент отражения параллельный плоскости падения
rov=(abs(rpe^2)+abs(rpa^2))/2;
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
function dl = DlinaVectora(vect)
s=0;
for k=1:length(vect)
    sv=vect(k);
    s=s+sv^2;
end
dl=sqrt(s);
end
function [ vy ] = VysLuchaPP(svo,svb,teta,hk,ak,bk,n1,n2,dn,dlv,razotv,phiz,phi,va,vb,vh,uab,uah,ubh,ide,vat,vbt,vabt,vht)
s=svb*rand(); po=0; rop=0; rov=0;
ugpalu=[0,0,-1]; ugpalu=PovorotizNovkStar(phi,teta,phiz,ugpalu); phipaluxy=atan(abs(ugpalu(2)/ugpalu(1)));
if (s<svo) %попали в основание
    phiup=acos(abs(ugpalu(3))); phiup=provUgla(abs(phiup)); %угол (направление) падения луча на основание
    mz=PoiskTetaShtrNach(dn,phiup,n1,n2); phiups=mz(1); pvon=mz(3); %phiups - угол преломления, phiup - угол падения
    if ((phiups>=0) && (pvon==0))
        phiups=provUgla(abs(phiups));
        rop=PoiskReflPerv(phiup,phiups,n1,n2);
        nalu=[ugpalu(1),ugpalu(2)]; %направление луча в плоскости XY (основание)
        topa=[ak*rand(),bk*rand(),0]; %точка падения на плоскость XY
        mz=VykhodIzOsnPP(nalu,topa,ak,bk,hk,phiups); v=mz(1); vpp=mz(4); plo=mz(2); ugbe=mz(3); %ищем, откуда будет выходить луч, в плоскости самого основания
    if (v==2) %1 - выход через основание, 2 - через бок. пов-ть
        phiups=pi/2-abs(phiups); %новый угол падения
        mz=PoiskVykhUglaSSKonech(dn,phiups,n1,n2); pvok=mz(3); 
    if (pvok==0)
    phiupss=pi/2-mz(1); %выходит через боковую поверхность
    rov=PoiskReflVtor(phiups,mz(1),n1,n2);
    po=PopNePop(dlv,abs(phiup-phiupss),razotv); %попадание в диапазон малых углов
    else %при полном внутреннем отражении даем второй шанс при отражении от боковой поверхности
        nlxyz=n1*[sign(ugpalu(1))*cos(phipaluxy)*sin(phiup),sign(ugpalu(2))*sin(phipaluxy)*sin(phiup),sign(ugpalu(3))*cos(phiup)]; %направление луча при падении
        nlxyzs=n2*[sign(nlxyz(1))*cos(phipaluxy)*sin(phiups),sign(nlxyz(2))*sin(phipaluxy)*sin(phiups),sign(nlxyz(3))*cos(phiups)]; %направление луча при преломлении
        nnrlo=[nlxyzs(1),nlxyzs(2)]; nnrlo=nnrlo/DlinaVectora(nnrlo); topaxy=[sign(nlxyzs(1))*abs(nnrlo(1)*plo)+topa(1),sign(nlxyzs(2))*abs(nnrlo(2)*plo)+topa(2)];
        switch (vpp)
            case 3
            nlxyzs(1)=-nlxyzs(1); veni=2; %отражение от правого, YZ
            case 4
            nlxyzs(1)=-nlxyzs(1); veni=2; %отражение от левого, YZ
            case 1 
            nlxyzs(2)=-nlxyzs(2); veni=1; %отражение от верхнего, XZ
            case 2 
            nlxyzs(2)=-nlxyzs(2); veni=1; %отражение от нижнего, XZ
        end
        topaz=DlinaVectora(topaxy)/tan(ugbe); tpl=[topaxy(1),topaxy(2),sign(nlxyz(3))*abs(topaz)];
        mz=VykhodIzOsnPPPVO(ak,bk,hk,nlxyzs,tpl,veni,phiups,vpp,3,1,0); %ищем, откуда будет выходить луч
        vkp=mz(5); ugobet=mz(2);
        if ((((vkp==1) || (vkp==2)) && ((vpp==3) || (vpp==4))) || (((vkp==3) || (vkp==4)) && ((vpp==1) || (vpp==2))))
        phiups=pi/2-abs(phiups); %новый угол падения после полного внутреннего отражения
        end
        mz=PoiskVykhUglaSSKonech(dn,phiups,n1,n2); pvokk=mz(3); phiupss=mz(1);
        if ((pvokk==0) && (ugobet>=0))
        rov=PoiskReflVtor(phiups,phiupss,n1,n2);
        po=PopNePop(dlv,abs(phiup-phiupss),razotv); %попадание в диапазон малых углов
        else
            po=0; rop=0; rov=0;
        end
    end
    else
        mz=PoiskVykhUglaSSKonech(dn,abs(phiups),n1,n2); pvok=mz(3);
        if (pvok==0)
            phiupss=mz(1); phiupss=provUgla(abs(phiupss));
            po=1; rov=PoiskReflVtor(phiups,phiupss,n1,n2); %через основание - никакого ухода луча нет
        else po=0; rov=0; rop=0;
        end
    end
    else rop=0; rov=0; po=0;
    end
    if (ide==9) %построение графика
        cem=1e2; mxa=0:(vat(1)/cem):vat(1); mya=0:(vat(2)/cem):vat(2); mxh=0:(vht(1)/cem):vht(1); myh=0:(vht(2)/cem):vht(2); %новые координаты %ska=0:(ak/cem):ak; skh=0:(hk/cem):hk; skh=-skh; skb=0:(bk/cem):bk; %старые координаты
        pl=plot(mxa+vbt(1),mya+vbt(2),':r',mxh+vbt(1),myh+vbt(2),':g',mxh+vabt(1),myh+vabt(2),':b',mxa+vbt(1)+vht(1),mya+vbt(2)+vht(2),':m'); %новые координаты %pl=plot(ska,skh,':b',ska,skh,':k'); %старые координаты
        set(pl,'LineWidth',3); hold on; grid on; xlabel({'x'}); ylabel({'y'}); title({'y(x)'}); hold off;
    end
else %попадание в боковую поверхность
    vyav=vh*va*sin(uah); vybv=vh*vb*sin(ubh); vypo=vyav+vybv;
    vxy=vypo*rand();
    if (vxy<=vyav) %1 - верхнее основание, XZ
        tpx=ak*rand(); 
        if (phiz>pi/2) 
        tpy=bk; 
        else tpy=0;
        end
        tpz=-hk*rand(); veni=1;
    else %2 - нижнее основание, YZ
        tpx=ak; tpy=bk*rand(); tpz=-hk*rand(); veni=2;
    end %координаты точки падения в старых координатах
       tpl=[tpx,tpy,tpz]; upp=opreUgPadPre(ugpalu,veni,dn,n1,n2); 
       phiup=upp(1); phiups=upp(2); narlu=[upp(3),upp(4)]; %ищем углы падения и преломления
       if (phiups>=0)
           rop=PoiskReflPerv(phiup,phiups,n1,n2);
           mz=VykhodIzOsnBokPovPP(ak,bk,hk,narlu,tpl,veni,phiups); v=mz(1); vpp=mz(4); plbp=mz(3); ugbe=mz(2); %ищем, откуда будет выходить луч, в плоскости самого основания
           phiups=provUgla(abs(phiups));
    switch (v)
        case 2 %боковая поверхность напротив
                mz=PoiskVykhUglaSSKonech(dn,phiups,n1,n2); pvok=mz(3);
        if (pvok==0)
                phiupss=abs(mz(1)); phiupss=provUgla(phiupss);
                po=1; rov=PoiskReflVtor(phiups,phiupss,n1,n2); %через основание - никакого ухода луча нет
        else po=0; rov=0;
        end
        case 1 %через верхнее или нижнее основание или боковую поверхность сбоку
                mz=PoiskVykhUglaSSKonech(dn,pi/2-phiups,n1,n2); pvok=mz(3);
                if (pvok==0)
                phiupss=pi/2-abs(mz(1)); phiupss=provUgla(phiupss);
                rov=PoiskReflVtor(pi/2-phiups,pi/2-phiupss,n1,n2);
                po=PopNePop(dlv,abs(phiup-phiupss),razotv); %попадание в диапазон малых углов
                else %при полном внутреннем отражении даем второй шанс при отражении от боковой поверхности и от оснований
                    narlu=ugpalu;
                    if (veni==1) %XZ: ak*hk
                    phixyz=abs(narlu(3)/narlu(1)); phixyz=atan(phixyz); phixyz=provUgla(phixyz);
                    nlxyz=n1*[sign(ugpalu(1))*cos(phixyz)*sin(phiup),sign(ugpalu(2))*cos(phiup),sign(ugpalu(3))*sin(phixyz)*sin(phiup)]; %направление луча при падении
                    nlxyzs=n2*[sign(nlxyz(1))*cos(phixyz)*sin(phiups),sign(nlxyz(2))*cos(phiups),sign(nlxyz(3))*sin(phixyz)*sin(phiups)]; %направление луча при преломлении
                    nnrlo=[nlxyzs(1),nlxyzs(3)]; nnrlo=nnrlo/DlinaVectora(nnrlo); 
                    topaxz=[sign(nlxyzs(1))*abs(nnrlo(1)*plbp)+tpl(1),sign(nlxyzs(3))*abs(nnrlo(2)*plbp)+tpl(3)]; %направление луча на боковой поверхности
                    switch (vpp)
            case 3
            nlxyzs(1)=-nlxyzs(1); venin=2; %отражение от правого, YZ
            case 4
            nlxyzs(1)=-nlxyzs(1); venin=2; %отражение от левого, YZ
            case 1 
            nlxyzs(3)=-nlxyzs(3); venin=1; %отражение от верхнего, XY
            case 2 
            nlxyzs(3)=-nlxyzs(3); venin=1; %отражение от нижнего, XY
                    end
                    topay=DlinaVectora(topaxz)/tan(ugbe); tpl=[topaxz(1),sign(nlxyzs(2))*abs(topay),topaxz(2)];
                    mz=VykhodIzOsnPPPVO(ak,hk,bk,nlxyzs,tpl,venin,phiups,vpp,2,0,1); %ищем, откуда будет выходить луч
                    elseif (veni==2) %YZ: bk*hk
                        phixyz=abs(narlu(3)/narlu(2)); phixyz=atan(phixyz); phixyz=provUgla(phixyz);
                        nlxyz=n1*[sign(ugpalu(1))*cos(phiup),sign(ugpalu(2))*cos(phixyz)*sin(phiup),sign(ugpalu(3))*sin(phixyz)*sin(phiup)]; %направление луча при падении
                        nlxyzs=n2*[sign(nlxyz(1))*cos(phiups),sign(nlxyz(2))*cos(phixyz)*sin(phiups),sign(nlxyz(3))*sin(phixyz)*sin(phiups)]; %направление луча при преломлении
                        nnrlo=[nlxyzs(2),nlxyzs(3)]; nnrlo=nnrlo/DlinaVectora(nnrlo); topayz=[sign(nlxyzs(2))*abs(nnrlo(1)*plbp)+tpl(2),sign(nlxyzs(3))*abs(nnrlo(2)*plbp)+tpl(3)];
                        switch (vpp)
            case 3
            nlxyzs(2)=-nlxyzs(2); venin=2; %отражение от правого, XZ
            case 4
            nlxyzs(2)=-nlxyzs(2); venin=2; %отражение от левого, XZ
            case 1 
            nlxyzs(3)=-nlxyzs(3); venin=1; %отражение от верхнего, XY
            case 2 
            nlxyzs(3)=-nlxyzs(3); venin=1; %отражение от нижнего, XY
                        end
                        topax=DlinaVectora(topayz)/tan(ugbe); tpl=[topax,topayz(1),topayz(2)];
                        mz=VykhodIzOsnPPPVO(bk,hk,ak,nlxyzs,tpl,venin,phiups,vpp,1,0,2); %ищем, откуда будет выходить луч
                    end
                    vkp=mz(5); ugobet=mz(2);
        if ((((vkp==1) || (vkp==2)) && ((vpp==3) || (vpp==4))) || (((vkp==3) || (vkp==4)) && ((vpp==1) || (vpp==2))))
        phiups=pi/2-abs(phiups); %новый угол падения после полного внутреннего отражения
        end
        mz=PoiskVykhUglaSSKonech(dn,phiups,n1,n2); pvokk=mz(3); phiupss=mz(1);
        if ((pvokk==0) && (ugobet>=0))
        rov=PoiskReflVtor(phiups,phiupss,n1,n2);
        po=PopNePop(dlv,abs(phiup-phiupss),razotv); %попадание в диапазон малых углов
        else
            po=0; rov=0; rop=0;
        end
                end
                    rov=0; po=0; rop=0;
        case 0
            po=0; rov=0; rop=0;
    end
       else rop=0; rov=0; po=0;
       end
end %точка падения луча на боковую поверхность
vy=[po,rop*rov];
end
function pr10 = PraCha10(phi,phis,n1,n2)
pr10=n1*sin(phi)-n2*sin(phis);
end
function pr20 = PraCha20(phis,phiss,n1,n2)
pr20=n2*sin(phis)-n1*sin(phiss);
end
%показатель преломления бромида калия
function opn = OprPokPreKBr(la)
lam=[1,2,10,20,25]; pp=[1.487,1.48,1.52,1.48,1.453]; 
ko=polyfit(lam,pp,length(pp)-1); opn=polyval(ko,la);
end
function [ vy ] = VykhodIzOsnPPPVO(ak,bk,hk,narlu,tpl,veni,phiups,vpp,pk,os,vn) %ищем, откуда будет выходить луч
vyhm=[0,0,0,0,0];
if (os==1) %основание - 1, боковая поверхность - 0
    vyhm=VybVykhOsnPVO(ak,bk,hk,tpl,narlu,phiups,veni,vpp,pk); %выбрана боковая со стороной b
elseif (os==0)
if (vn==1) %ak*hk 
    vyhm=VybVykhOsnPVO(ak,hk,bk,tpl,narlu,phiups,veni,vpp,pk); %выбрана боковая со стороной b
elseif (vn==2) %bk*hk
    vyhm=VybVykhOsnPVO(bk,hk,ak,tpl,narlu,phiups,veni,vpp,pk); %выбрана боковая со стороной a
end
end
    vy=vyhm; %1 - выход через основание, 2 - через бок. пов-ть
end
function vy = VybVykhOsnPVO(ak,bk,hk,tpl,nlxyz,phiups,veni,vpp,peko)
vyh=1;
uvav=abs(bk-tpl(2)); uvav=atan(abs(ak/uvav)); %угол, под которым видна верхняя часть стороны a
uvnk=abs(ak/tpl(2)); uvnk=atan(uvnk); %угол видимости начала координат, 2
uvnk1=abs(bk/tpl(1)); uvnk1=atan(uvnk1); %угол видимости начала координат, 1
uvav1=abs(bk/(ak-tpl(1))); uvav1=atan(uvav1);
unl=abs(nlxyz(2)/nlxyz(1)); unl=atan(unl); %направление распространения луча
mz=opredPutLuchaPVO(nlxyz,unl,uvav,uvnk,uvav1,uvnk1,tpl,ak,bk,veni,vpp); 
pulu=mz(1); vnp=mz(2); vkp=mz(3); ugolbeta=-1;
if (pulu>0)
ugolbeta=pulu/abs(hk-tpl(peko)); ugolbeta=atan(ugolbeta); ugolbeta=provUgla(ugolbeta);
if (ugolbeta>phiups)
    vyh=1;
else vyh=2; %1 - выход через основание, 2 - через бок. пов-ть
end
else
    vyh=0;
end
vy=[vyh,ugolbeta,pulu,vnp,vkp];
end
function [ pl ] = opredPutLuchaPVO(nlxyz,unl,uvav,uvnk,uvav1,uvnk1,tpl,ak,bk,veni,vpp)
pulu=-1; vnp=0; vkp=0;
if ((nlxyz(1)>=0) && (nlxyz(2)>=0))
    if ((vpp==4) && (veni==2)) %боковые стороны
        unl=pi/2-unl;
if ((unl>=0) && (unl<=uvav))
    rov=abs(bk-tpl(2)); pulu=rov/cos(unl); vkp=1;
end
if ((unl>uvav) && (unl<=(pi/2)))
    rov=abs(ak); pulu=rov/sin(unl); vkp=3;
end
vnp=vpp;
    elseif ((vpp==2) && (veni==1))
if ((unl>=0) && (unl<=uvav1))
    rov=abs(ak-tpl(1)); pulu=rov/cos(unl); vkp=3;
end
if ((unl>uvav1) && (unl<=(pi/2)))
    rov=abs(bk); pulu=rov/sin(unl); vkp=1;
end
vnp=vpp;
    end
end
if ((nlxyz(1)<0) && (nlxyz(2)<0))
    if ((veni==2) && (vpp==3)) %боковые стороны
        unl=pi/2-unl;
if ((unl>=uvnk) && (unl<=(pi/2)))
    rov=abs(ak); pulu=rov/sin(unl); vkp=4;
end
if ((unl<uvnk) && (unl>=0))
    rov=abs(bk-tpl(2)); pulu=rov/cos(unl); vkp=2;
end
vnp=vpp;
    elseif ((veni==1) && (vpp==1))
if ((unl>=uvnk1) && (unl<=(pi/2)))
    rov=abs(bk); pulu=rov/sin(unl); vkp=2;
end
if ((unl<uvnk1) && (unl>=0))
    rov=abs(ak-tpl(1)); pulu=rov/cos(unl); vkp=4;
end
vnp=vpp;
    end
end
if ((nlxyz(1)<0) && (nlxyz(2)>=0))
    if ((veni==2) && (vpp==3)) %боковые стороны
        unl=pi/2-unl;
    if ((unl<uvav) && (unl>=0))
        rov=abs(bk-tpl(2)); pulu=rov/cos(unl); vkp=1;
    end
    if ((unl>=uvav) && (unl<=(pi/2)))
    rov=abs(ak); pulu=rov/sin(unl); vkp=4;
    end
    vnp=vpp;
    elseif ((veni==1) && (vpp==2))
    if ((unl>=uvnk1) && (unl<=(pi/2)))
        rov=abs(bk); pulu=rov/sin(unl); vkp=1;
    end
    if ((unl<uvnk1) && (unl>=0))
    rov=abs(ak-tpl(1)); pulu=rov/cos(unl); vkp=4;
    end
    end
    vnp=vpp;
end
if ((nlxyz(1)>=0) && (nlxyz(2)<0))
    if ((veni==2) && (vpp==4)) %боковые стороны
        unl=pi/2-unl;
    if ((unl>=uvnk) && (unl<=(pi/2)))
        rov=abs(ak); pulu=rov/sin(unl); vkp=3;
    end
    if ((unl<uvnk) && (unl>=0))
    rov=abs(bk-tpl(2)); pulu=rov/cos(unl); vkp=2;
    end
    vnp=vpp;
    elseif ((veni==1) && (vpp==1))
    if ((unl>=uvav1) && (unl<=(pi/2)))
        rov=abs(bk); pulu=rov/sin(unl); vkp=2;
    end
    if ((unl<uvav1) && (unl>=0))
    rov=abs(ak-tpl(1)); pulu=rov/cos(unl); vkp=3;
    end
    end
    vnp=vpp;
end
if ((isnan(nlxyz(1))) || (isnan(nlxyz(2))))
    pulu=-1; vkp=0; vnp=0;
end
if ((isinf(nlxyz(1))) || (isinf(nlxyz(2))))
    pulu=-1; vkp=0; vnp=0;
end
pl=[pulu,vnp,vkp];
end
%поиск угла преломления, phi>0
function [ pesdp ] = PoiskTetaShtrNach(dnr,phipad,n1,n2)
phipre=0; pvo=0; a=0; b=pi/2; ep=1e-6; Nit=1e3; k=0; ra=abs(a-b); dephi=0; %находим угол преломления
if (dnr<0)
    phipre=asin(n2/n1); phipre=provUgla(phipre);
end
if (((phipad<phipre) && (dnr<0)) || (dnr>0))
while ((ra>ep) && (k<Nit))
    c=(a+b)/2;
    fa=PraCha10(phipad,a,n1,n2);
    fb=PraCha10(phipad,b,n1,n2);
    fc=PraCha10(phipad,c,n1,n2);
    ro=vybPoKo(fa,fb,fc,a,b,c); 
    a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
    phiprel=provUgla(c);
else
    pvo=1; phiprel=pi/2; dephi=2*pi;
end
pesdp=[phiprel,dephi,pvo,phipre];
end
%поиск угла выхода: v=1 - выход через основание, v=2 - через бок. пов-ть
function [ pesdp ] = PoiskVykhUglaSSKonech(dnr,phi,n1,n2)
phipre=0; pvo=0; a=0; b=pi/2; ep=1e-7; Nit=1e3; k=0; ra=abs(a-b); dephi=0; %находим угол преломления
if (dnr<=0)
    if (n1<=n2)
        phipre=asin(n1/n2); phipre=provUgla(phipre);
    else
        phipre=asin(n2/n1); phipre=provUgla(phipre);
    end
end
if (((phi<phipre) && (dnr<0)) || (dnr>0))
while ((ra>ep) && (k<Nit))
    c=(a+b)/2;
    fa=PraCha20(phi,a,n1,n2);
    fb=PraCha20(phi,b,n1,n2); 
    fc=PraCha20(phi,c,n1,n2);
    ro=vybPoKo(fa,fb,fc,a,b,c); 
    a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
    phiss=provUgla(c); dephi=abs(phi-phiss);
else
    pvo=1; phiss=pi/2; dephi=2*pi;
end
pesdp=[phiss,dephi,pvo,phipre];
end
function ug = provUgla(ugol)
ugo=ugol;
if (ugo>pi/2)
    ugo=pi-ugo;
end
if (ugo<0)
    ugo=abs(ugo);
end
ug=ugo;
end
%------------------пока не трогать--------------
function pd = PopNePop(lam,deltatetass,diaot)
p=0; deko=lam/diaot; %deko=1.24*lam/diaot;
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
    v=[a,b,c];
end
function [ alsr ] = SredGrafRass(vyfr,vyko,xv)
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
%vyfr - выбор фракции: 1 - 60, 2 - 60-100, 3 - 100-150, 4 - 150-200; vyko - выбор массовой доли: 1 - 0,1 %, 2 - 0,2 %, 3 - 0,3 %
switch (vyfr)
    case 1
Tal=privedkEdiPropus(TrkoVer5311()); Tal2=privedkEdiPropus(TrkoVer5312());Tal3=privedkEdiPropus(TrkoVer5313()); 
Tal4=privedkEdiPropus(TrkoVer53101()); Tal5=privedkEdiPropus(TrkoVer53102()); Tal6=privedkEdiPropus(TrkoVer53103()); 
alsre=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv,vyko);
    case 4
Tal=privedkEdiPropus(TrkoVer5681()); Tal2=privedkEdiPropus(TrkoVer5682()); Tal3=privedkEdiPropus(TrkoVer5683()); 
Tal4=privedkEdiPropus(TrkoVer5691()); Tal5=privedkEdiPropus(TrkoVer5692()); Tal6=privedkEdiPropus(TrkoVer5693());
alsre=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv,vyko);
    case 2
Tal=privedkEdiPropus(TrkoVer5421()); Tal2=privedkEdiPropus(TrkoVer5422()); Tal3=privedkEdiPropus(TrkoVer5423()); 
Tal4=privedkEdiPropus(TrkoVer5431()); Tal5=privedkEdiPropus(TrkoVer5432()); Tal6=privedkEdiPropus(TrkoVer5433());
Tal7=privedkEdiPropus(TrkoVer5441()); Tal8=privedkEdiPropus(TrkoVer5442()); Tal9=privedkEdiPropus(TrkoVer5443());
alsre=opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv,vyko);
    case 3
Tal=privedkEdiPropus(TrkoVer5551()); Tal2=privedkEdiPropus(TrkoVer5552()); Tal3=privedkEdiPropus(TrkoVer5553()); 
Tal4=privedkEdiPropus(TrkoVer5561()); Tal5=privedkEdiPropus(TrkoVer5562()); Tal6=privedkEdiPropus(TrkoVer5563());
Tal7=privedkEdiPropus(TrkoVer5571()); Tal8=privedkEdiPropus(TrkoVer5572()); Tal9=privedkEdiPropus(TrkoVer5573());
alsre=opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv,vyko);
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
function [ al ] = opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv,vyko)
nTal=length(Tal);
for k=1:nTal 
    Ta(k)=-log(Tal(k))/xv(1); Ta2(k)=-log(Tal2(k))/xv(2); Ta3(k)=-log(Tal3(k))/xv(3);
    Ta4(k)=-log(Tal4(k))/xv(4); Ta5(k)=-log(Tal5(k))/xv(5); Ta6(k)=-log(Tal6(k))/xv(6); 
end
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
function [ al ] = opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv,vyko)
nTal=length(Tal);
for k=1:nTal
    Ta(k)=-log(Tal(k))/xv(1); Ta2(k)=-log(Tal2(k))/xv(2); Ta3(k)=-log(Tal3(k))/xv(3);
    Ta4(k)=-log(Tal4(k))/xv(4); Ta5(k)=-log(Tal5(k))/xv(5); Ta6(k)=-log(Tal6(k))/xv(6);
    Ta7(k)=-log(Tal7(k))/xv(7); Ta8(k)=-log(Tal8(k))/xv(8); Ta9(k)=-log(Tal9(k))/xv(9); 
end
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
a=1e-5; b=1e5; ep=1e-6; Nit=1e2; k=0; ffv=1e-4; ra=abs(a-b); %ffv = 0,01 % - такая доля частиц не попадает в диапазон фракции
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; fa=erfc(a)-ffv; fb=erfc(b)-ffv; fc=erfc(c)-ffv; 
    ro=vybPoKo(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3); k=k+1;
    ra=abs(a-b);
end
pk=c;
end
%di - для поиска СКО, ide - цилиндр (1) или прямоуг. пар-д (2)
function ras60 = RasFra60(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[250.239,249.740,249.223,250.395,250.336,249.55]; 
mv=[0.283,0.464,0.812,0.22,0.547,0.777]; 
tol=[0.73,0.72,0.72,0.72,0.71,0.7];
rokbr=2.75; rov=0.49; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %в см3
vv(k)=mv(k)/(1e3*rov); %в см3
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %в мкм
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
mini=0; maxi=6e1; mo=(maxi+mini)/2; %МО
si=abs(mini-maxi)/di; %СКО
disp('Фракция менее 60 мкм');
depp=PoisDn(1,nom,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras60=depp; %ras60=alphaRas(tol,length(tol),Tr);
end
function ras60_100 = RasFra60_100(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[250,249.629,249.294,249.706,249.510,249.307,250.328,249.604,249.206]; 
mv=[0.255,0.539,0.809,0.295,0.517,0.756,0.36,0.534,0.843]; 
tol=[0.72,0.71,0.7,0.7,0.73,0.72,0.74,0.7,0.76]; 
rokbr=2.75; rov=0.52; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
mini=6e1; maxi=1e2; mo=(maxi+mini)/2; %МО
si=abs(mini-maxi)/di; %СКО
disp('Фракция 60-100 мкм');
depp=PoisDn(2,nom,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras60_100=depp; %ras60_100=alphaRas(tol,length(tol),Tr);
end
function ras100_150 = RasFra100_150(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[249.913,249.607,249.218,249.929,249.695,249.306,250.405,249.625,249.348];
mv=[0.315,0.473,0.709,0.293,0.528,0.83,0.27,0.493,0.764]; 
tol=[0.74,0.74,0.72,0.72,0.71,0.7,0.78,0.73,0.76]; 
rokbr=2.75; rov=0.53; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
vv(k)=vv(k)*((1e4)^3); %в мкм3
end
mini=1e2; maxi=15e1; mo=(maxi+mini)/2; %МО
si=abs(mini-maxi)/di; %СКО
disp('Фракция 100-150 мкм');
depp=PoisDn(3,nom,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras100_150=depp; %ras100_150=alphaRas(tol,length(tol),Tr);
end
function ras150_200 = RasFra150_200(di,popr,dv,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[250.882,249.590,249.213,250.299,249.441,249.365];
mv=[0.320,0.533,0.849,0.223,0.502,0.797]; 
tol=[0.76,0.72,0.69,0.73,0.73,0.73]; 
rokbr=2.75; rov=0.56; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); vv(k)=mv(k)/(1e3*rov); xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vv(k)=vv(k)*((1e4)^3); vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
end
mini=15e1; maxi=2e2; mo=(maxi+mini)/2; %МО
si=abs(mini-maxi)/di; %СКО
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
function [ vy ] = VykhodIzOsnPP(nlxy,tpl,ak,bk,hk,phiups)
vyh=1;
uvbl=bk-tpl(2); uvbl=atan(abs(uvbl/tpl(1))); %угол, под которым видна левая часть стороны b
uvnk=tpl(2)/tpl(1); uvnk=atan(abs(uvnk)); %угол видимости начала координат
uvav=bk-tpl(2); uvav=abs(uvav/(ak-tpl(1))); uvav=atan(uvav); %угол, под которым видна сторона a верхняя
uvbp=ak-tpl(1); uvbp=abs(tpl(2)/uvbp); uvbp=atan(uvbp); %угол видимости правой стороны b
unl=abs(nlxy(2)/nlxy(1)); unl=atan(unl); %направление распространения луча
mz=opredPutLucha(nlxy,unl,uvbl,uvnk,uvav,uvbp,tpl,ak,bk); pulu=mz(1); panapo=mz(2);
if (pulu>=0)
ugolbeta=pulu/hk; ugolbeta=atan(ugolbeta);
if (ugolbeta>phiups)
    vyh=2;
else vyh=1; %1 - выход через основание, 2 - через бок. пов-ть
end
else vyh=0;
end
vy=[vyh,pulu,ugolbeta,panapo];
end
function [ pl ] = opredPutLucha(nlxy,unl,uvbl,uvnk,uvav,uvbp,tpl,ak,bk)
pulu=-1; padpov=0;
if ((nlxy(1)>=0) && (nlxy(2)>=0))
if ((unl>=uvav) && (unl<=(pi/2)))
    rov=abs(bk-tpl(2));
    pulu=rov/sin(unl);
    padpov=1;
end
if ((unl<uvav) && (unl>=0))
    rov=abs(ak-tpl(1));
    pulu=rov/cos(unl);
    padpov=3;
end
end
if ((nlxy(1)<0) && (nlxy(2)<0))
if ((unl>=uvnk) && (unl<=(pi/2)))
    rov=abs(tpl(2));
    pulu=rov/sin(unl);
    padpov=2;
end
if ((unl<uvnk) && (unl>=0))
    rov=abs(tpl(1));
    pulu=rov/cos(unl);
    padpov=4;
end
end
if ((nlxy(1)<0) && (nlxy(2)>=0))
    if ((unl>=uvbl) && (unl<=(pi/2)))
        rov=abs(bk-tpl(2));
        pulu=rov/sin(unl);
        padpov=1;
    end
if ((unl<uvbl) && (unl>=0))
    rov=abs(tpl(1));
    pulu=rov/cos(unl);
    padpov=4;
end    
end
if ((nlxy(1)>=0) && (nlxy(2)<0))
    if ((unl>=uvbp) && (unl<=(pi/2)))
        rov=abs(tpl(2));
        pulu=rov/sin(unl);
        padpov=2;
    end
    if ((unl<uvbp) && (unl>=0))
    rov=abs(ak-tpl(1));
    pulu=rov/cos(unl);
    padpov=3;
    end
end
if ((isnan(nlxy(1))) || (isnan(nlxy(2))))
    pulu=-1;
end
if ((isinf(nlxy(1))) || (isinf(nlxy(2))))
    pulu=-1;
end
pl=[pulu,padpov];
end
function [ sk ] = PovorotizNovkStar(phi,teta,phiz,nk)
nk=PovorotObratnonaUgolPhi(nk,phi);
nk=PovorotObratnonaUgolTeta(nk,teta);
nk=PovorotObratnonaVokrugOZ(nk,phiz);
sk=nk;
end
function [ koor ] = PovorotObratnonaUgolPhi(koors,phi)
xs=koors(1); ys=koors(2); zs=koors(3);
x=xs; y=ys*cos(phi)+zs*sin(phi); z=-ys*sin(phi)+zs*cos(phi); %первый поворот в плоскости YZ вокруг оси X
koor=[x,y,z];
end
function [ koor ] = PovorotObratnonaUgolTeta(koors,teta)
xs=koors(1); ys=koors(2); zs=koors(3);
x=xs*sin(teta)+zs*cos(teta); y=ys; z=-xs*cos(teta)+zs*sin(teta); %второй поворот в плоскости XZ вокруг оси Y
koor=[x,y,z];
end
function [ nk ] = PovorotObratnonaVokrugOZ(koors,phiz)
xs=koors(1); ys=koors(2); zs=koors(3);
x=xs*cos(phiz)-ys*sin(phiz); y=xs*sin(phiz)+ys*cos(phiz); z=zs; %третий поворот в плоскости XY вокруг оси Z
nk=[x,y,z];
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
function [ upapr ] = opreUgPadPre(ugpalu,veni,dn,n1,n2) %угол (направление) падения луча
switch (veni)
    case 2 %2 - нижнее основание YZ - b, (-1,0,0); 1 - верхнее XZ - a, (0,-1,0)
        phiup=acos(ugpalu(1)); phiup=provUgla(phiup); %работаем в плоскости YZ
        nalu=[ugpalu(2),ugpalu(3)]; pxyz=1; %направление луча на боковой поверхности YZ
    case 1 
        phiup=acos(ugpalu(2)); phiup=provUgla(phiup); %работаем в плоскости XZ
        nalu=[ugpalu(1),ugpalu(3)]; pxyz=2; %направление луча на боковой поверхности XZ
end
    mz=PoiskTetaShtrNach(dn,phiup,n1,n2); pvo=mz(3);
    if (pvo==0)
    phiups=mz(1); phiups=provUgla(phiups); %phiups - угол преломления, phiup - угол падения
    else
        phiups=-1; %полное внутреннее отражение
    end
    upapr=[phiup,phiups,nalu(1),nalu(2),pxyz]; %точка падения на плоскость XY и направление луча
end
function [ vy ] = VykhodIzOsnBokPovPP(ak,bk,hk,narlu,tpl,veni,phiups) %tpl - координаты точки падения в конечных координатах XYZ
vyh=1; %1 - выход через основание, 2 - через бок. пов-ть
if (veni==2) %в нижней плоскости - YZ 
    vyhm=VybVykhOsnBokPov(bk,hk,ak,tpl,narlu,phiups); 
else %в верхней плоскости - XZ 
    vyhm=VybVykhOsnBokPov(ak,hk,bk,tpl,narlu,phiups); 
end
    vy=vyhm;
end
function vy = VybVykhOsnBokPov(ak,hk,bk,tpl,nlxy,phiups)
vyh=1;
uvhl=abs(hk-tpl(2)); uvhl=atan(abs(uvhl/tpl(1))); %угол, под которым видна левая часть стороны h
uvnk=abs(tpl(2)/tpl(1)); uvnk=atan(uvnk); %угол видимости начала координат
uvav=abs(hk-tpl(2)); uvav=uvav/abs(ak-tpl(1)); uvav=atan(uvav); %угол, под которым видна сторона a (b) верхняя
uvhp=abs(ak-tpl(1)); uvhp=abs(tpl(2)/uvhp); uvhp=atan(uvhp); %угол видимости правой стороны h
unl=abs(nlxy(2)/nlxy(1)); unl=atan(unl); %направление распространения луча
mz=opredPutLucha(nlxy,unl,uvhl,uvnk,uvav,uvhp,tpl,ak,hk); 
pulu=mz(1); kpp=mz(2); ugolbeta=0;
if (pulu>0)
ugolbeta=pulu/bk; ugolbeta=atan(ugolbeta); ugolbeta=provUgla(ugolbeta);
if (ugolbeta>phiups)
    vyh=1;
else vyh=2; %1 - выход через основание, 2 - через бок. пов-ть
end
else
    vyh=0;
end
vy=[vyh,ugolbeta,pulu,kpp];
end
function [ dv ] = dlinavolny()
dl=dlvoVer53101(); 
p=length(dl); 
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dv = dl;
end
function [ np ] = PokazPrelMas()
np=[1.00011907361118
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