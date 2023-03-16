%рассеяние только для цилиндров - дополнительная функция
function [ t ] = tmp46()
format long g;
npp=Kramers_ver()-1; npp=npp';
na=3.1; ko=27; de=1e-1;
dv=1.4:de:ko;
dvr=na:de:ko;
q=1; np=0;
for k=1:length(dv)
    if (dv(k)>dvr(q))
        np(q)=npp(k-1);
        q=q+1;
    end
end
fidTr=fopen('Cylindry_Tr-1.txt','w'); fiddpp=fopen('Cylindry_dpp-1.txt','w'); 
fidTrs=fopen('Cylindry_Trs-1.txt','w'); 
fiddnra=fopen('Cylindry_dnra-1.txt','w'); fidrazdT=fopen('Cylindry_razdT-1.txt','w'); 
fidn1=fopen('Cylindry_nKBr-1.txt','w'); fiddv=fopen('Cylindry_dl_vo-1.txt','w'); 
np(q)=npp(length(npp)); 
no=3 %no - номер фракции
nom=2 %nom - номер концентрации
%k=VyvodTrs(dvr,no,nom,fidTrs);
for k=1:length(dv)
    n1(k)=OprPokPreKBr(dv(k));
end
n1=n1'; tr=0; np=np'; npma=max(np); npmi=min(np);
for k=1:length(np)
    trk(k)=nacha(np(k),dvr(k),no,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
    tr(k)=trk(1);
end
fclose(fidTr); fclose(fiddpp); fclose(fiddnra); fclose(fidrazdT); fclose(fidn1); fclose(fidTrs);
t=tr;
end
function [ tm ] = nacha(popr,dv,no,nom,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
format long g; %определяет коэффициент рассеяния для таблетки KBr + вермикулит
%t=2*PoisKorn()
t=2*erfcinv(1e-4)
switch (no) %no - номер фракции
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
function [ dlv ] = dlinyvoln()
format long g;
dl=0; dl=dlvoVer53101(); 
p=length(dl);
for k=1:p  
    dl(k)=1e-2/dl(k);
end
dlv = dl;
end
function [ pk ] = PoisDn(vyfr,vyko,xv,vv,mo,si,mi,ma,dpp,dvr,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
TrRef=SredGrafRass(vyfr,vyko,xv); dlv=(1e6)*dlinyvoln();
q=length(TrRef)
q=length(dlv)
si=si'
numele=length(dlv)-1; na=1; dvko=dvr; n1=OprPokPreKBr(dvr)
xve=opredEffPar(xv,vyfr,vyko); vve=opredEffPar(vv,vyfr,vyko); mo=mo'
dnre=0; q=1; eps=1e-2; TR=1e0; chit=1e1; dpp=dpp'
while (dvr<=dvko)
    w=PoiskNomera(dlv,dvr); dlvo=dlv(w); TRs=TrRef(w-1)+(TrRef(w)-TrRef(w-1))*(dvr-dlv(w-1))/(dlv(w)-dlv(w-1))
    Tra=TRs; Trb=1e0; sc=1; razdT=abs(Tra-Trb); dnc=1; TR=(Tra+Trb)/2;
while ((razdT>eps) && (sc<chit))
    TR=(Tra+Trb)/2; dnc=1;
    dna=npmi; dnb=npma; ep=1e-5; Nit=1e2; k=0; ra=abs(dna-dnb);
while ((ra>ep) && (k<Nit))
    dnc=(dna+dnb)/2;
    na2=n1+dna;
    if (dna==0)
        fa=1;
    else
    fa=opreKoefProp(vve,mo,si,mi,ma,dlvo,n1,na2,1); 
    end
    fa=fa-TR;
    nb2=n1+dnb;
    if (dnb==0)
        fb=1;
    else
    fb=opreKoefProp(vve,mo,si,mi,ma,dlvo,n1,nb2,2); 
    end
    fb=fb-TR;
    nc2=n1+dnc;
    if (dnc==0)
        fc=1;
    else
    fc=opreKoefProp(vve,mo,si,mi,ma,dlvo,n1,nc2,3); 
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
fprintf(fiddv,'%0.15f\n',dvr);
pk=dnre;
end
%Данная функция должна возвращать вероятность попадания при выстрелах, mo - матем. ожид-е, si - СКО, mi, ma - мин. и макс. разм.
function mdn = opreKoefProp(vvk,mo,si,mi,ma,dlvo,n1,n2,nabc)
Nr=1e0; raot=13e3; blis=0; ocp=0; sepoo=0; rosr=0; %raot - размер отверстия, blis - число удачных попаданий
for g=1:Nr %усреднение по реализациям
suv=0; k=1; nk=1e9;
while ((suv<vvk) && (k<nk))
hk=mo+si*randn(); hk=proverka(hk,mi,ma); %определение длины направляющей
tetani=pi*rand(); %направление излучения - угол между плоскостью основания и осью Z (вверх)
phini=pi*(-1/2+rand()); %угол между проекцией оси цилиндра на горизонтальную ось и осью X - коллинеарно направлению излучения (ось X вниз)
nalu=[0    0   -1]; nalusk=PovorotyIzNovKIskh(nalu,tetani,phini); %направление излучения в старых координатах
ocp=ocp+1;
    dk=mo+si*randn(); dk=proverka(dk,mi,ma); rk=dk/2; %определение диаметра цилиндра
    sob=opredPloshOsnBokPov(rk,tetani,phini,hk);
    suv=suv+pi*hk*(rk^2); %проецирование на вертикальную ось, перпендикулярную направлению излучения
    se=sob(1); sb=sob(2); maxx=sob(3); maxy=sob(4); phima=sob(6); phimi=sob(7); pc=[sob(8)        sob(9)      sob(10)];
    sep=se+sb; sepoo=sepoo+sep; %сечение рассеяния
    vys=VysLucha(se,sep,tetani,hk,rk,n1,n2,dlvo,raot,phini,nalusk/DlinaVectora(nalusk),phima,phimi,maxx,maxy,pc,nabc); 
    se=proverkaZnach(vys(1)); se=abs(se); sb=proverkaZnach(vys(2)); sb=abs(sb);
    if ((se<2) && (sb<=1))
    blis=blis+se; %число благополучных исходов, 1 - успешно, 0 - нет
    rosr=rosr+sb; %коэффициент отражения
    else
        continue;
    end
k=k+1;
end
end
w=blis/ocp; rosr=rosr/Nr; plos=pi*(raot^2)/4;
p=sepoo/Nr/plos;
mdn=p*w*(1-rosr)+(1-p); %пропускание за счет рассеяния
end
%выстрел луча по цилиндру - в плоскости угла падения
function [ vy ] = VysLucha(se,spo,teta,hk,rk,n1,n2,dlv,razotv,phi,nisk,phima,phimi,maxx,maxy,pc,iden)
s=spo*rand(); po=0; rop=1; rov=1;
if (s<se)
    pno=PadenOsn(teta,hk,rk,n1,n2,dlv,razotv,phi,nisk,phima,phimi,maxx,maxy,iden); po=pno(1); rop=pno(2); rov=pno(3);
else
    pvbp=PadenBokovPover(teta,hk,rk,n1,n2,dlv,razotv,phi,nisk,phima,phimi,maxx,maxy,pc,iden); po=pvbp(1); rop=pvbp(2); rov=pvbp(3);
end
rov=rop*rov;
vy=[po          rov];
end
function phixy = PrivkOdnoZnach(ug,phinul)
phix1=real(acos(ug(1)));
phix2=real(asin(ug(2)));
if (phix1>(pi/2))
    if (phix2<0)
        phix1=pi-phix1; %x<0, y<0, 3 четверть
        phix1=pi+phix1;
        phix2=abs(phix2);
        phix2=pi+phix2;
    else
        phix2=pi-phix2; %x<0, y<0, 2 четверть
    end
else
    if (phix2<0)
        phix1=-phix1+2*pi; %x>0, y<0, 4 четверть
        phix2=phix2+2*pi;
    end
end
phix1=phix1';
phix2=phix2';
phixy=(phix1+phix2)/2+phinul;
end
function [ kup ] = PoisKTPNOsn(rk,phi,teta,phima,phimi,maxx,maxy,iden)
mazn=PoisKTPNPodProg(rk,phi,teta,phima,phimi,maxx,maxy);
av=mazn(1)/2+maxx/2; bv=mazn(2)/2+maxy/2; avm=mazn(3); av=abs(av)/2+abs(avm)/2; bvm=mazn(4); bv=abs(bv)/2+abs(bvm)/2; phimax=mazn(5); phimix=mazn(6);
kup=PoisKTPNPodProgOsnVtor(rk,phi,teta,phimax,phimix,av,bv,avm,bvm,iden,2*pi*rand());
end
function [ mazn ] = PoisKTPNPodProg(r,phi,teta,phima,phimi,maxx,maxy)
vxyz=[r*cos(phima)     r*sin(phima)    0]; vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); vxyz(3)=0; vxyz=Pereobozn(vxyz); 
vxyz=PlusOdinPovorot(vxyz); k=1; phim(k)=phima; xym(k)=vxyz(2); %vxyz=vxyz';
phima=phima+pi; vxyz=[r*cos(phima)      r*sin(phima)    0]; vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); vxyz(3)=0; vxyz=Pereobozn(vxyz); 
vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phima; xym(k)=vxyz(2); %vxyz=vxyz';
phimi=(phima-pi/2-pi)/2+phimi/2; vxyz=[r*cos(phimi)      r*sin(phimi)    0]; vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); vxyz(3)=0; vxyz=Pereobozn(vxyz); 
vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phimi; xym(k)=vxyz(2); %vxyz=vxyz';
phimi=phimi+pi; vxyz=[r*cos(phimi)      r*sin(phimi)    0]; vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); vxyz(3)=0; vxyz=Pereobozn(vxyz); 
vxyz=PlusOdinPovorot(vxyz); k=k+1; phim(k)=phimi; xym(k)=vxyz(2); %vxyz=vxyz'; xym=xym';
if (xym(1)>xym(2))
        minx=xym(2); phimix=phim(2); maxx=(xym(1)+maxx)/2; phimax=phim(1);
else
        minx=xym(1); phimix=phim(1); maxx=(xym(2)+maxx)/2; phimax=phim(2); 
end
if (xym(3)>xym(4))
        miny=xym(4); maxy=(xym(3)+maxy)/2; phimay=phim(3); phimiy=phim(4);
else
        miny=xym(3); maxy=(xym(4)+maxy)/2; phimay=phim(4); phimiy=phim(3);
end
mazn=[maxx      maxy        minx        miny        phimax      phimix      phimay      phimiy];
end
function [ mazn ] = PoisKTPNPodProgOsnVtorPogProgPerv(r,phimix,phimax,teta,phi,maxx,maxy,minx,miny,ugphi,iden)
%--поиск углов
nuk=0; ra=(ugphi-phimix)/(phimax-phimix); nxy=1e2; hnu=(phimax-phimix)/nxy; xyz=[r*cos(phimix)       r*sin(phimix)        0]; xyz=PovorotyIzIskhKNov(xyz,teta,phi); xyz(3)=0; xymi=xyz; xyz=[r*cos(phimax)       r*sin(phimax)        0]; 
xyz=PovorotyIzIskhKNov(xyz,teta,phi); xyz(3)=0; xyma=xyz; rxy=(xymi+xyma)/2; phis=(phimax+phimix)/2; rxys=[r*cos(phis)       r*sin(phis)        0]; rxys=PovorotyIzIskhKNov(rxys,teta,phi); rxys(3)=0;
if (phimax<phimix)
    phimix=phimix-2*pi; ugphi=ra*(phimax-phimix)+phimix; 
end
    if ((rxys(1)>rxy(1)) && (iden==2))
        phimax=phimix; phimix=phimix-pi; ugphi=ugphi-pi;
    end
if (iden==1) %1 - попали в основание, 2 - в боковую поверхность
    xyz=PoisKTPNPodProgOsnVtorPogProg(phimix,phimax,0,0); phimix=xyz(1); phimax=xyz(2);
    xyz=PoisKTPNPodProgOsnVtorPogProg(phimix,phimax,ugphi,0); phimix=xyz(1); phimax=xyz(2); ugphi=xyz(3); 
    if (phimax>2*pi)
        phimax=phimax-2*pi; phimix=phimix-2*pi; ugphi=ugphi-2*pi;
    end
end
if (((ugphi>phimax) || (ugphi<phimix)) && (iden==2))
            ugphi=ra*(phimax-phimix)+phimix; 
end
if ((ugphi>phimax) || (ugphi<phimix))
            ugphi=ra*(phimax-phimix)+phimix;
end
xyz=[r*cos(phimax)       r*sin(phimax)        0]; xyzmas=xyz'; xyz=[r*cos(phimix)       r*sin(phimix)        0]; xyzmis=xyz'; 
raxys=xyzmas-xyzmis; draxys=DlinaVectora(raxys); psiras=PrivkOdnoZnach([real(raxys(1)/draxys)         real(raxys(2)/draxys)],0); %в старых координатах
xyz=[r*cos(phimax)       r*sin(phimax)        0]; xyz=PovorotyIzIskhKNov(xyz,teta,phi); xyz(3)=0; xyz=Pereobozn(xyz); 
xyzma=xyz'; dxyzma=DlinaVectora(xyzma); psima=PrivkOdnoZnach([real(xyzma(1)/dxyzma)         real(xyzma(2)/dxyzma)],0); 
xyzmi=-xyzma; psimi=psima-pi; raxy=xyzma-xyzmi; draxy=DlinaVectora(raxy); psira=PrivkOdnoZnach([real(raxy(1)/draxy)         real(raxy(2)/draxy)],0); %в новых координатах
maxzx=(abs(maxx)+abs(minx))/2; maxzy=(abs(maxy)+abs(miny))/2; spsis=(psimi+psima)/2; 
naspsis=[cos(spsis)      sin(spsis)      0]; phis=(phimax+phimix)/2;  nasphis=[cos(phis)      sin(phis)      0]; ugpsisphis=acos((naspsis*nasphis')/DlinaVectora(naspsis)/DlinaVectora(nasphis)); 
ugpsi=psimi+(psima-psimi)*rand(); %выбор угла в основании
xyz=PoisKTPNPodProgOsnVtorPogProg(psimi,psima,ugpsi,psira); 
psimi=xyz(1); psima=xyz(2); ugpsi=xyz(3); psira=xyz(4);
hnu=(phimax-phimix)/nxy; hpsi=(psima-psimi)/nxy;
mazn=[phimix        phimax      psimi       psima      psira       psiras      maxzx       maxzy       hnu         hpsi       ugpsisphis        ugpsi       ugphi       nuk];
end
function [ mazn ] = PoisKTPNPodProgOsnVtorPogProgVtor(ugpsi,psira,maxzx,maxzy,ugpsisphis,hpsi,hnu,ugnu,r,phimax,phimix,phi,teta,psimi,psiras,iden)
%----прямое преобразование из старых координат в новые
psik=ugpsi-psira; xyz=[maxzx*cos(psik)      maxzy*sin(psik)       0]; xyz=DopolnitPovorot(xyz,psira);
    if ((ugpsisphis>(pi/2)) || (phimax==0) || (sign(hpsi)==-sign(hnu)))
        xyz=DopolnitPovorot(xyz,pi);
    end
    xs=xyz(1); ys=xyz(2); %новые координаты
    xyz=[r*cos(ugnu)       r*sin(ugnu)        0]; 
    if (sign(hpsi)==-sign(hnu))
        xyz=DopolnitPovorot(xyz,pi);
    end
    xv=xyz(1); yv=xyz(2); %старые координаты
    xyz=PovorotyIzIskhKNov(xyz,teta,phi); xyz(3)=0; xyz=Pereobozn(xyz); 
    xn=xyz(1); yn=xyz(2); %новые координаты
%-----обратное преобразование из новых координат в старые
    if (sign(hpsi)==-sign(hnu))
        xyz=DopolnitPovorot(xyz,pi);
    end
xyz=DopolnitPovorot(xyz,-psira); 
nuk=PrivkOdnoZnach([real(xyz(1)/maxzx)         real(xyz(2)/maxzy)],psira+(phimix-psimi)); nuk=real(nuk);
xyz=[-max([maxzx  maxzy])*sin(nuk)      max([maxzx  maxzy])*cos(nuk)       0];
        if ((ugpsisphis>pi/2) || (phimax==0))
            xyz=DopolnitPovorot(xyz,pi);
        end
        xw=xyz(1); yw=xyz(2); %старые координаты
        nuk=PrivkOdnoZnach([(xw/r)         (yw/r)],0); nuk=real(nuk);
        mz=PoisKTPNPodProgOsnVtorPogProg(phimix,phimax,nuk,0); nuk=mz(3); k=1; Nit=1e2; koe=1;
        if (iden==1) %попали в основание
            koe=1;
        else koe=2;
        end
while ((nuk>phimax) && (k<Nit))
            nuk=nuk-koe*pi; k=k+1;
end
mazn=[xs        ys      xv      yv      xn      yn      xw      yw      nuk         ugnu        ugpsi]; %xs, ys (1,2 - НП), xn, yn (5,6 - ПВК) - новые, xv, yv (3,4), xw, yw (7,8 - НП) - старые
end
function [ maug ] = PoisKTPNPodProgOsnVtorPogProg(psimi,psima,ugpsi,ugol)
Nit=1e2; k=1;
    while (((ugpsi>psima) || (ugpsi<psimi)) && (k<Nit))
        if (psimi>ugpsi)
        psimi=psimi-pi;
        psima=psima-pi;
        end
        if (psimi>psima)
        phite=psimi;
        psimi=psima;
        psima=phite;
        ugol=ugol+pi;
        end
    if (psima<ugpsi)
        psima=psima+pi;
        psimi=psimi+pi;
    end
    if (psimi>psima)
        phite=psimi;
        psimi=psima;
        psima=phite;
        ugol=ugol+pi;
    end
    k=k+1;
    end
    if (ugpsi>(2*pi))
        ugpsi=ugpsi-2*pi;
        psimi=psimi-2*pi;
        psima=psima-2*pi;
    end
    k=1;
    while ((ugol>(2*pi)) && (k<Nit))
        ugol=ugol-2*pi;
        k=k+1;
    end
    maug = [psimi       psima       ugpsi       ugol];
end
function [ vesk ] = PoisKTPNPodProgOsnVtor(r,phi,teta,phimax,phimix,maxx,maxy,minx,miny,iden,ugphi)
mz=PoisKTPNPodProgOsnVtorPogProgPerv(r,phimix,phimax,teta,phi,maxx,maxy,minx,miny,ugphi+phimix,1);
phimix=mz(1); phimax=mz(2); psimi=mz(3); psima=mz(4); psiras=mz(6); maxzx=mz(7); maxzy=mz(8); hnu=mz(9); hpsi=mz(10); ugpsisphis=mz(11); ugpsi=mz(12); ugphi=mz(13); psira=mz(5); 
mz=0; mz=PoisKTPNPodProgOsnVtorPogProgVtor(ugpsi,psira,maxzx,maxzy,ugpsisphis,hpsi,hnu,ugphi,r,phimax,phimix,phi,teta,psimi,psiras,1);
xk=mz(1); yk=mz(2); xw=mz(7); yw=mz(8); nuk=mz(9); ugphi=mz(10);
dolya=rand(); xw=dolya*xw; yw=dolya*yw; xk=dolya*xk; yk=dolya*yk;
if (iden==4)
    cem=1e2; mxw=0:xw/cem:xw; myw=0:yw/cem:yw; %старые координаты
    mxk=0:xk/cem:xk; myk=0:yk/cem:yk; %новые координаты
        nnu=cem; hnun=2*pi/nnu; nu=(phimax-pi):hnun:(phimax+pi); iota=(psima-pi):hnun:(psima+pi); %nu - в старых, iota - в новых
for k=1:length(nu)   
    ugpsi=iota(k); ugnu=nu(k);
    mz=0; mz=PoisKTPNPodProgOsnVtorPogProgVtor(ugpsi,psira,maxzx,maxzy,ugpsisphis,hpsi,hnu,ugnu,r,phimax,phimix,phi,teta,psimi,psiras,1);
    xs(k)=mz(1); ys(k)=mz(2); xv(k)=mz(3); yv(k)=mz(4); xn(k)=mz(5); yn(k)=mz(6); xw(k)=mz(7); yw(k)=mz(8); mz=0;
end
        %pl=plot(xn,yn,':b',xs,ys,':k'); %новые координаты
        %pl=plot(xv,yv,':b',xw,yw,':k'); %старые координаты
        subplot(1,2,1); pl=plot(xw,yw,':b',mxw,myw,':k'); %старые координаты
        set(pl,'LineWidth',3); hold on; grid on; xlabel({'x'}); ylabel({'y'}); title({'y(x)'}); hold off;
        subplot(1,2,2); pl=plot(xs,ys,':b',mxk,myk,':k'); %новые координаты
        set(pl,'LineWidth',3); hold on; grid on; xlabel({'x'}); ylabel({'y'}); title({'y(x)'}); hold off;
end
%xw=xv; yw=yv;
        if ((ugphi>=phimix) && (ugphi<=phimax))
vesk=[xw        yw      0       phimax      phimix      nuk         ugphi];
        else
vesk=[-xw         -yw       0       (phimax+pi)      (phimix+pi)        nuk         ugphi];
        end
end
function [ vesk ] = PoisKTPNPodProgBokVtor(r,phi,teta,phimax,phimix,maxx,maxy,minx,miny,rpc,iden)
mz=PoisKTPNPodProgOsnVtorPogProgPerv(r,phimix,phimax,teta,phi,maxx,maxy,minx,miny,phimix+(phimax-phimix)*rand(),2);
phimix=mz(1); phimax=mz(2); psimi=mz(3); psima=mz(4); psira=mz(5); psiras=mz(6); maxzx=mz(7); maxzy=mz(8); hnu=mz(9); hpsi=mz(10); ugpsisphis=mz(11); ugpsi=mz(12); ugphi=mz(13); mz=0;
mz=PoisKTPNPodProgOsnVtorPogProgVtor(ugpsi,psira,maxzx,maxzy,ugpsisphis,hpsi,hnu,ugphi,r,phimax,phimix,phi,teta,psimi,psiras,2);
xk=mz(1); yk=mz(2); xw=mz(7); yw=mz(8); nuk=mz(9); ugphi=mz(10);
dolya=rand(); rpcp=dolya*rpc; rpcp=PovorotyIzNovKIskh(rpcp,teta,phi); %поиск вектора-переноса
        if ((ugphi>=phimix) && (ugphi<=phimax))
vesk=[xw+rpcp(1)        yw+rpcp(2)      rpcp(3)       phimax      phimix        nuk         ugphi];
        else
vesk=[-xw+rpcp(1)         -yw+rpcp(2)       rpcp(3)       (phimax+pi)      (phimix+pi)        nuk         ugphi];
        end
end
function [ ktpl ] = PoisKTPNBokPodProg(r,phi,teta,phima,phimi,maxx,maxy,rpc,iden)
    mazn=PoisKTPNPodProg(r,phi,teta,phima,phimi,maxx,maxy); maxx=mazn(1); maxy=mazn(2); minx=mazn(3); miny=mazn(4); phimax=phima/2+mazn(5)/2; phimix=mazn(6);
    ktpl=PoisKTPNPodProgBokVtor(r,phi,teta,phimax,phimix,maxx,maxy,minx,miny,rpc,iden);
end
function [ kpl ] = PoisKTPNBok(r,phi,teta,phima,phimi,maxx,maxy,pc,iden)
    kpl=PoisKTPNBokPodProg(r,phi,teta,phima,phimi,maxx,maxy,pc,iden);
end
%падает и выходит через основание или боковую стенку
function pr10 = PraCha10(phi,phis,n1,n2)
pr10=n1*sin(phi)-n2*sin(phis);
end
function pr20 = PraCha20(phis,phiss,n1,n2)
pr20=n2*sin(phis)-n1*sin(phiss);
end
%показатель преломления бромида калия
function opn = OprPokPreKBr(la)
lam=[1      2       10      20      25]; 
pp=[1.487   1.48        1.52        1.48        1.453]; 
ko=polyfit(lam,pp,length(pp)-1); 
opn=polyval(ko,la);
end
%поиск угла преломления, phi>0
function [ pesdp ] = PoiskPhiShtrNach(phi,n1,n2)
dnr=n2-n1; phipre=0; pvo=0; a=0; b=pi/2; ep=1e-7; Nit=1e2; k=0; ra=abs(a-b); dephi=0; %находим угол преломления
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
    fa=PraCha10(phi,a,n1,n2);
    fb=PraCha10(phi,b,n1,n2);
    fc=PraCha10(phi,c,n1,n2);
    ro=vybPoKo(fa,fb,fc,a,b,c); 
    a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
    phis=provUgla(c);
else
    pvo=1; phis=pi/2; dephi=2*pi;
end
pesdp=[phis     dephi       pvo         phipre];
end
function [ pesdp ] = PoiskPhiShtrVykh(phi,n1,n2)
dnr=n1-n2; phipre=0; pvo=0; a=0; b=pi/2; ep=1e-7; Nit=1e2; k=0; ra=abs(a-b); dephi=0; %находим угол преломления
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
    phiss=provUgla(c);
else
    pvo=1; phiss=pi/2; dephi=2*pi;
end
pesdp=[phiss    dephi       pvo         phipre];
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
function [ mpv ] = PlusOdinPovorot(axy)
x=axy(1); y=axy(2); z=axy(3);
phiksi=x/y; phiksi=atan(phiksi); 
%debo=DlinaVectora([x      y]); phiksi=PrivkOdnoZnach([x/debo   y/debo],0);
ksi=x*cos(phiksi)-y*sin(phiksi); eta=x*sin(phiksi)+y*cos(phiksi); dzeta=z;
mpv=[ksi    eta     dzeta];
end
function [ mpv ] = Pereobozn(axy)
x=axy(1); y=axy(2); z=axy(3);
ksi=y; eta=-x; dzeta=z;
mpv=[ksi    eta         dzeta];
end
function [ mpv ] = VraschenKoorVykh(axy,ugol)
x=axy(1); y=axy(2); z=axy(3);
xs=-x*cos(2*ugol)-y*sin(2*ugol); ys=x*sin(2*ugol)-y*cos(2*ugol); zs=z;
mpv=[xs    ys   zs];
end
function [ mpv ] = VraschenKoorVykhObrat(axy,ugol)
xs=axy(1); ys=axy(2); zs=axy(3);
x=-xs*cos(2*ugol)+ys*sin(2*ugol); y=-xs*sin(2*ugol)-ys*cos(2*ugol); z=zs;
mpv=[x          y           z];
end
function [ mpv ] = DopolnitVraschen(axy,ugvr)
x=axy(1); y=axy(2); z=axy(3);
xs=x*cos(ugvr)-y*sin(ugvr); ys=x*sin(ugvr)+y*cos(ugvr); zs=z;
mpv=[xs    ys   zs];
end
function [ mpv ] = DopolnitVraschenObrat(axy,ugvr)
xs=axy(1); ys=axy(2); zs=axy(3);
x=xs*cos(ugvr)+ys*sin(ugvr); y=-xs*sin(ugvr)+ys*cos(ugvr); z=zs;
mpv=[x          y         z];
end
function [ mpv ] = DopolnitPovorot(axy,phi)
x=axy(1); y=axy(2); z=axy(3);
ksi=x*cos(phi)-y*sin(phi); eta=x*sin(phi)+y*cos(phi); dzeta=z;
mpv=[ksi    eta         dzeta];
end
function [ v ] = VykhodIzCylOsn(nalu,tpisk,teta,phi,h,r,phiups)
    dlo=PoiskDlPutiLuchaOsn(r,nalu,tpisk); dlo=abs(dlo); vyh=1;
    beta=atan(dlo/h); beta=provUgla(beta); phiups=provUgla(phiups);
    if (phiups>beta)
        vyh=2; %выход из боковой поверхности
    else vyh=1; %выход из основания
    end
    %vyh=OprTochVykhKritSluch(phi,teta,nalu,tpisk,h,r,vyh); %в плоскости основания
v=[vyh      beta        dlo];
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
    v=[a        b       c];
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
ro=0; a=1e-5; b=1e5; ep=1e-6; Nit=1e2; k=0; ffv=1e-4; ra=abs(a-b); %ffv = 0,01 % - такая доля частиц не попадает в диапазон фракции
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; 
    fa=erfc(a)-ffv; 
    fb=erfc(b)-ffv; 
    fc=erfc(c)-ffv; 
    ro=vybPoKo(fa,fb,fc,a,b,c); 
    a=ro(1); b=ro(2); c=ro(3); k=k+1;
    ra=abs(a-b);
end
pk=c;
end
function [ ras60 ] = RasFra60(di,popr,dv,noko,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
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
mini=0; maxi=6e1; mo=(maxi+mini)/2; %МО
si=abs((mini-maxi)/di); %СКО
disp('Фракция менее 60 мкм');
depp=PoisDn(1,noko,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras60=depp; 
end%di - для поиска СКО
function [ ras60_100 ] = RasFra60_100(di,popr,dv,noko,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
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
mini=6e1; maxi=1e2; mo=(maxi+mini)/2; %МО
si=abs((mini-maxi)/di); %СКО
disp('Фракция 60-100 мкм');
depp=PoisDn(2,noko,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras60_100=depp; 
end
function [ ras100_150 ] = RasFra100_150(di,popr,dv,noko,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
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
mini=1e2; maxi=15e1; mo=(maxi+mini)/2; %МО
si=abs((mini-maxi)/di); %СКО
disp('Фракция 100-150 мкм');
depp=PoisDn(3,noko,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras100_150=depp;
end
function [ ras150_200 ] = RasFra150_200(di,popr,dv,noko,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[250.882 249.590 249.213 250.299 249.441 249.365];
mv=[0.320 0.533 0.849 0.223 0.502 0.797]; 
tol=[0.76 0.72 0.69 0.73 0.73 0.73]; 
rokbr=2.75; rov=0.56; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); vv(k)=mv(k)/(1e3*rov); xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vv(k)=vv(k)*((1e4)^3); vkbr(k)=vkbr(k)*((1e4)^3); %в мкм3
end
mini=15e1; maxi=2e2; mo=(maxi+mini)/2; %МО
si=abs((mini-maxi)/di); %СКО
disp('Фракция 150-200 мкм');
depp=PoisDn(4,noko,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras150_200=depp; 
end
function p = proverka(hk,mi,ma)
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
%phin - угол направления падения, xp, yp - координаты падения на круге, ищем в плоскости основания
function dlot = PoiskDlPutiLuchaOsn(rok,nalu,tpisk)
ko=nalu(2)/nalu(1); tpisk(3)=0; dlo=0;
if ((isinf(ko)) || (isnan(ko)))
else
xtope(1)=0; ytope(1)=0; xtope(2)=0; ytope(2)=0;
yp=tpisk(2); xp=tpisk(1); sc=yp-ko*xp;
b=(2*ko*sc); a=(abs(ko)^2+1); c=(abs(sc)^2-abs(rok)^2); diskr=(abs(b)^2)-(4*a*c);
if (diskr<0)
    rvp(1)=0; rvp(2)=0; rvv(1)=0; rvv(2)=0; rvp(3)=0; rvv(3)=0;
else
    diskr=sqrt(diskr);
    xtope(1)=(-b-diskr)/2/a; 
    ytope(1)=ko*xtope(1)+sc; 
    rvp=[xtope(1)       ytope(1)        0];
    xtope(2)=(-b+diskr)/2/a; 
    ytope(2)=ko*xtope(2)+sc; 
    rvv=[xtope(2)      ytope(2)        0];
end
rvp=rvp'; rvv=rvv'; tpisk=tpisk';
a=nalu*rvp;
b=nalu*rvv;
if (a>0)
    dxyz=rvp-tpisk; 
    dlo=DlinaVectora(dxyz);
else
if (b>0)
    dxyz=rvv-tpisk;
    dlo=DlinaVectora(dxyz);
end
end
end
dlot=abs(dlo);
end
function [ nisk ] = PovorotyIzNovKIskh(vnk,teta,phi)
ksi=vnk(1); eta=vnk(2); dzeta=vnk(3);
xs=ksi; ys=eta*cos(phi)-dzeta*sin(phi); zs=eta*sin(phi)+dzeta*cos(phi);
x=xs*sin(teta)+zs*cos(teta); y=ys; z=-xs*cos(teta)+zs*sin(teta);
nisk=[x     y   z];
end
function [ vnk ] = PovorotyIzIskhKNov(vik,teta,phi)
x=vik(1); y=vik(2); z=vik(3);
xs=x*sin(teta)-z*cos(teta); ys=y; zs=x*cos(teta)+z*sin(teta);
ksi=xs; eta=ys*cos(phi)+zs*sin(phi); dzeta=-ys*sin(phi)+zs*cos(phi);
vnk=[ksi     eta   dzeta];
end
function dlve = DlinaVectora(vect)
s=0;
for k=1:length(vect)
    s=s+(vect(k))^2;
end
dlve=sqrt(s);
end
function ugol = ProvUglaNu(nu,phimax,phimix)
if (nu<phimix)
        nu=nu+pi;
end
    if (nu>phimax)
        nu=nu-pi;
    end
    ugol = nu;
end
function [ mz ] = PoisTochPeresech(nalus,rk,iden,ktpos)
phiyx=atan(nalus(2)/nalus(1)); ktpos(3)=0; nalus(3)=0; ko=tan(phiyx); dxy=0; nalus=nalus/DlinaVectora(nalus);
if ((isinf(ko)) || (isnan(ko)))
else
xtope=0; ytope=0; %ищем точки пересечения
yp=ktpos(2); xp=ktpos(1); sc=yp-ko*xp;
b=(2*ko*sc); a=(abs(ko)^2+1); c=(abs(sc)^2-abs(rk)^2); diskr=(abs(b)^2)-(4*a*c);
if (diskr<0)
    rvp=[0        0       0]; rvv=rvp; dxy=0;
else
    ep=1e-8; diskr=sqrt(diskr); xtope=(-b-diskr)/2/a; ytope=ko*xtope+sc; rvp=[xtope       ytope        0];
    xtope=(-b+diskr)/2/a; ytope=ko*xtope+sc; rvv=[xtope      ytope        0];
    if (DlinaVectora(rvp-ktpos)>ep)
    rt=rvp; rvp=rvv; rvv=rt;
    end
end
end
dvp=rvv-ktpos; dxy=DlinaVectora(dvp);
mz=[dxy       rvv(1)      rvv(2)      nalus(1)        nalus(2)      phiyx];
end
%падение луча в боковую поверхность
function [ pvbp ] = PadenBokovPover(teta,hk,rk,n1,n2,dlv,razotv,phi,nisk,phima,phimi,maxx,maxy,pc,iden) 
    mz=PoisKTPNBok(rk,phi,teta,phima,phimi,maxx,maxy,pc,iden); ktpsk=[mz(1)         mz(2)       mz(3)]; 
    phimax=mz(4); phimix=mz(5); nu=mz(6); ugphi=mz(7); nu=ProvUglaNu(nu,phimax,phimix); 
    %disp('пад. на бок. пов.: s>se'); 
    ktpos=[ktpsk(1)       ktpsk(2)         0]; %проверка пройдена - длина вектора rk
    phitp=-atan(ktpos(2)/ktpos(1)); phitp=phitp+(1-sign(ktpos(1)))*pi/2; phitp=phitp+(1-sign(phitp))*pi; phitp=pi-PrivkOdnoZnach(ktpos/rk,0)/2+phitp/2; %угол вращения
    nisk=DopolnitVraschen(nisk,phitp); ktpos=DopolnitVraschen(ktpos,phitp); ktpsk=DopolnitVraschen(ktpsk,phitp); %дополнительные вращения систем координат - точка на оси X
po=2; rop=1; rov=1;
if (nisk(1)<0) %nisk(1)=-nisk(1); nisk(2)=-nisk(2); %меняем на противоположное направление излучения, т.к. оно должно идти против оси X, сохраняем правую тройку
    nu=acos((nisk*ktpos')/DlinaVectora(nisk)/DlinaVectora(ktpos)); nu=provUgla(nu); nu=provUgla(atan(DlinaVectora([nisk(2)          nisk(3)])/abs(nisk(1))))/2+nu/2; %угол падения
        mz=0; mz=PoiskPhiShtrNach(nu,n1,n2); nus=provUgla(mz(1)); ppvovh=mz(3); %угол преломления в плоскости падения
        if (ppvovh==0)
        rop=PoiskReflPerv(nu,nus,n1,n2);
        phiyz=provUgla(atan(nisk(2)/nisk(3))); nalus=[sign(nisk(1))*cos(nus)            sign(nisk(2))*sin(nus)*sin(phiyz)        sign(nisk(3))*sin(nus)*cos(phiyz)]; %проверка пройдена: y-,z-составляющие nalus, nisk одинаковы 
        if (nalus(3)<0)
            roz=abs(hk)-abs(ktpsk(3)); roz=abs(roz); %поиск расстояния до основания
        else
            roz=abs(ktpsk(3));
        end
        nisk=n1*nisk; nalus=n2*nalus;
        niskn=nisk/DlinaVectora(nisk); nalusn=nalus/DlinaVectora(nalus);
        mz=0; mz=PoisTochPeresech(nalus,rk,iden,ktpos); dxy=mz(1); ktvl(1)=mz(2); ktvl(2)=mz(3); nlxy(1)=mz(4); nlxy(2)=mz(5); nlxy(3)=0; phiyx=mz(6);
        beta=provUgla(atan(roz/dxy)); gamma=DlinaVectora([nalusn(1)         nalusn(2)]); gamma=provUgla(atan(abs(nalusn(3))/gamma));
        if (beta>gamma) 
            %disp('вых. ч. бок. пов.');
            kztv=dxy*tan(gamma)*sign(nalusn(3))/DlinaVectora(nlxy); %kztv=kztv*DlinaVectora(nalus)/DlinaVectora([nalus(1)       nalus(2)]);
        ktvl(3)=ktpsk(3)+kztv; %излучение выходит через боковую поверхность, проверка пройдена, точка выхода лежит на боковой поверхности
        naluss=(VraschenKoorVykh(nalus,phiyx)+DopolnitVraschen(nalus,pi-2*phiyx))/2;
        ktvl=(VraschenKoorVykh(ktvl,phiyx)+DopolnitVraschen(ktvl,pi-2*phiyx))/2;
        provt=DlinaVectora([ktvl(1)         ktvl(2)])-rk; %проверка пройдена - точка выхода лежит на окружности  основания
        phiyz=provUgla(atan(naluss(2)/naluss(3)));
        nusv=provUgla(acos(naluss(1)/DlinaVectora(naluss))); %угол падения в точке выхода
        mz=0; mz=PoiskPhiShtrVykh(nusv,n1,n2); nuss=provUgla(mz(1)); ppvovy=mz(3); %угол преломления в плоскости падения
        if (ppvovy==0) 
        rov=PoiskReflVtor(nusv,nuss,n1,n2);
        phiyz=provUgla(atan(nalus(2)/nalus(3)));
        naluvy=n1*[sign(naluss(1))*cos(nuss),sign(naluss(2))*sin(nuss)*sin(phiyz),sign(naluss(3))*sin(nuss)*cos(phiyz)]; %тангенсальные составляющие сохраняются
        naluvyn=naluvy/DlinaVectora(naluvy);
        naluvyn=VraschenKoorVykhObrat(naluvyn,phiyx);
        naluvyskn=DopolnitVraschenObrat(naluvyn,phitp);
        denunus=acos(niskn*naluvyskn'/DlinaVectora(niskn)/DlinaVectora(naluvyskn));
        po=PopNePop(dlv,denunus,razotv); %попадание в диапазон малых углов - выход через бок. пов-ть
        else
            po=0; rop=1; rov=1;
        end
        else
            %disp('вых. ч. осн.');
            roxy=roz/tan(gamma);
            nalusn=nalus/DlinaVectora(nalus); nlxy=nalusn; nlxy(3)=0; nlxy=nlxy/DlinaVectora(nlxy);
            kztv=roxy*tan(gamma)*sign(nlxy(3))/DlinaVectora(nlxy); %provt=DlinaVectora(nlxy)-1;  %проверка пройдена - длина вектора равна одному
            ktvy(1)=ktpsk(1)+roxy*nalus(1); ktvy(2)=ktpsk(2)+roxy*nalus(2);          
            ktvy(3)=ktpsk(3)+kztv; %излучение выходит через основание %provt=DlinaVectora([ktvy(1)         ktvy(2)])-rk; %проверка пройдена, точка выхода внутри круга
            nusv=provUgla(acos(nalus(3)/DlinaVectora(nalus))); %угол падения внутри цилиндра на основание из боковой поверхности
            nusv=nusv/2+provUgla(atan(DlinaVectora([nalus(1),nalus(2)])/abs(nalus(3))))/2;
            mz=0; mz=PoiskPhiShtrVykh(nusv,n1,n2); nuss=provUgla(mz(1)); ppvovyh=mz(3);
            if (ppvovyh==0)
            rov=PoiskReflVtor(nusv,nuss,n1,n2); %коэффициент отражения
            phixy=provUgla(atan(nalus(1)/nalus(2)));
            naluvy=[sign(nalus(1))*sin(nuss)*sin(phixy)            sign(nalus(2))*sin(nuss)*cos(phixy)        sign(nalus(3))*cos(nuss)]; %тангенсальные составляющие сохраняются
            naluvyn=naluvy/DlinaVectora(naluvy);
            naluvysk=DopolnitVraschenObrat(naluvyn,phitp); %naluv=[sign(nalus(1))*sin(phixy)*sin(nuss)             sign(nalus(2))*cos(phixy)*sin(nuss)           sign(nalus(3))*cos(nuss)]; %naluvn=naluv/DlinaVectora(naluv);
            denunus=acos((niskn*naluvysk')/DlinaVectora(niskn)/DlinaVectora(naluvysk));
            po=PopNePop(dlv,denunus,razotv); %попадание в диапазон малых углов - выход через бок. пов-ть
            else
                po=0; rop=1; rov=1;
            end
        end
        else
            po=0; rop=1; rov=1;
        end
end
pvbp = [po      rop         rov];
end
%падение луча в основание - выход  через  боковую поверхность
function [ vcbp ] = VykhCherezBokPovOsn(nisk,dlo,rk,n1,n2,phiup,phiups,naluskoc,ktpsk,dlv,razotv,iden)
%disp('вых. ч. бок. пов.'); %луч достигает быстрее боковой поверхности
        niskn=nisk/DlinaVectora(nisk); 
        phixy=provUgla(atan(naluskoc(1)/naluskoc(2))); %направление луча в основании цилиндра
        nalu=n2*[sign(nisk(1))*sin(phiups)*sin(phixy)          sign(nisk(2))*sin(phiups)*cos(phixy)          sign(nisk(3))*cos(phiups)];
        nalun=nalu/DlinaVectora(nalu); provt=DlinaVectora(nalun)-1; %проверка пройдена - длина вектора равна одному
        phiups=provUgla(phiups); 
        zk=sign(nisk(3))*abs(dlo)/tan(phiups); 
        nlxy=nalu; nlxy(3)=0;  nlxy=nlxy/DlinaVectora(nlxy);
        ktvy(1)=ktpsk(1)+nlxy(1)*abs(dlo); 
        ktvy(2)=ktpsk(2)+nlxy(2)*abs(dlo); 
        %zk=zk*DlinaVectora([nlxy(1)         nlxy(2)])/DlinaVectora([nalun(1)        nalun(2)]);
        ktvy(3)=ktpsk(3)+zk;
        ugpo=PrivkOdnoZnach([ktvy(1)/rk             ktvy(2)/rk],0); %угол поворота, нацеливаемся на точку выхода
        nalus=DopolnitVraschen(nalu,ugpo);
        ktvys=DopolnitVraschen(ktvy,ugpo);
        if (iden==4)
        nnu=1e2; hnu=2*pi/nnu; nu=0:hnu:2*pi; hxy=dlo/nnu;
for k=1:length(nu)
            x0(k)=rk*cos(nu(k)); x1(k)=ktpsk(1)+k*hxy*naluskoc(1);
            y0(k)=rk*sin(nu(k)); y1(k)=ktpsk(2)+k*hxy*naluskoc(2);
end
        pl=plot(x0,y0,':b',x1,y1,':k'); 
        set(pl,'LineWidth',3); 
        hold on; grid on; 
        xlabel({'x'}); ylabel({'y'}); title({'y(x)'});
        end
        ugpa=provUgla(acos((nalus*ktvys')/DlinaVectora(nalus)/DlinaVectora(ktvys))); %угол падения внутри цилиндра на боковую поверхность
        phiyz=provUgla(atan(nalus(2)/nalus(3)));
    mz=0; mz=PoiskPhiShtrVykh(ugpa,n1,n2); phiupss=provUgla(mz(1)); ppvovy=mz(3); %угол преломления на выходе из цилиндра, луч выходит через боковую поверхность
    if (ppvovy==0)
    naluss=n1*[sign(nalus(1))*cos(phiupss)          sign(nalus(2))*sin(phiupss)*sin(phiyz)          sign(nalus(3))*sin(phiupss)*cos(phiyz)]; 
    nalusk=DopolnitVraschenObrat(naluss,ugpo); 
    naluskn=nalusk/DlinaVectora(nalusk); provt=DlinaVectora(naluskn)-1; %проверка пройдена - длина вектора равна одному
    deususs=acos((naluskn*niskn')/DlinaVectora(niskn)/DlinaVectora(naluskn));
    po=PopNePop(dlv,deususs,razotv); %попадание в диапазон малых углов
        rop=PoiskReflPerv(phiup,phiups,n1,n2);
        rov=PoiskReflVtor(ugpa,phiupss,n1,n2);
    else
        po=0; rop=1; rov=1;
    end
        vcbp = [po      rop         rov];
end
%луч падает в основание цилиндра
function [ pno ] = PadenOsn(teta,hk,rk,n1,n2,dlv,razotv,phi,nisk,phima,phimi,maxx,maxy,iden)
    po=2; rop=1; rov=1; %начальные значения 
    %disp('пад. на осн.: s<se');
    mz=PoisKTPNOsn(rk,phi,teta,phima,phimi,maxx,maxy,iden); 
    ktpsk=[mz(1)        mz(2)       mz(3)]; %координаты точки падения луча на основание
    phimax=mz(4); phimix=mz(5); nu=mz(6); ugphi=mz(7); nu=ProvUglaNu(nu,phimax,phimix); provt=DlinaVectora(ktpsk)-rk; %отрицательна - проверка пройдена - точка входа лежит внутри окружности  основания
    if ((provt<=0) && (nisk(3)<0)) %направление луча должно быть против оси Z, если это не так, то это ошибка
    naluskoc=nisk; naluskoc(3)=0; naluskoc=naluskoc/DlinaVectora(naluskoc); %направление излучения на основании цилиндра
    phiup=provUgla(acos(nisk(3)/DlinaVectora(nisk))); phiups=DlinaVectora([nisk(1)             nisk(2)])/abs(nisk(3)); phiup=phiup/2+provUgla(atan(phiups))/2; %угол падения
    mz=0; mz=PoiskPhiShtrNach(phiup,n1,n2); phiups=provUgla(mz(1)); ppvovh=mz(3); %phiups - угол преломления, phiup - угол падения - проверен
    if (ppvovh==0)
    mz=VykhodIzCylOsn(naluskoc,ktpsk,teta,phi,hk,rk,phiups); v=mz(1); dlo=abs(mz(3)); %ищем, откуда будет выходить луч, в плоскости самого основания
    if (dlo>0)
    if (v==2) %выход через боковую поверхность
        mz=0; mz=VykhCherezBokPovOsn(nisk,dlo,rk,n1,n2,phiup,phiups,naluskoc,ktpsk,dlv,razotv,iden); po=mz(1); rop=mz(2); rov=mz(3);
    else 
        %disp('вых. ч. осн.');
        mz=0; mz=PoiskPhiShtrNach(phiup,n1,n2); phiups=provUgla(mz(1)); ppvos=mz(3); 
        mz=0; mz=PoiskPhiShtrVykh(phiups,n1,n2); phiupss=provUgla(mz(1)); ppvoss=mz(3);
        if ((ppvos==0) && (ppvoss==0))
            po=1;  %через основание - никакого ухода луча нет
            rop=PoiskReflPerv(phiup,phiups,n1,n2); 
            rov=PoiskReflVtor(phiups,phiupss,n1,n2); 
        else 
            rop=1; rov=1; po=0;
        end
    end
    else
        po=1; rop=0; rov=0; %%если точка лежит за окружностью или на ней, то излучение проходит свободно
    end
    else 
        po=0; rop=1; rov=1;
    end
    end %повторяем процедуру, т.к. это явно ошибка
    pno=[po         rop         rov];
end
function t = proverkaZnach(d)
tt=d;
if (isnan(d))
    tt=2;
end
if (isinf(d))
    tt=2;
end
t=tt;
end
function [ pobp ] = opredPloshOsnBokPov(r,teta,phi,h)
pc=[0  0   -h];
pc=PovorotyIzIskhKNov(pc,teta,phi);
rpc=DlinaVectora(pc); %перенос центра
phima=PoiskUglaMax(phi,teta,r);
vxyz=[r*cos(phima)  r*sin(phima)   0]; 
vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); %находим полуоси эллипса - проекции основания на вертикальной плокости
vxyz(3)=0; 
ma=DlinaVectora(vxyz); %большая полуось
phimi=phima-pi/2;
vxyz=[r*cos(phimi)  r*sin(phimi)   0];
vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); 
vxyz(3)=0;
mi=DlinaVectora(vxyz); %малая полуось
so=pi*abs(ma*mi); 
sbp=abs(rpc*ma);
phimi=phima+pi;
pobp=[so    sbp         ma       mi          rpc         phima       phimi          pc(1)       pc(2)       pc(3)];
end
function rop = PoiskReflPerv(phieta,phietas,n1,n2)
rpa=(n2*cos(phieta)-n1*cos(phietas))/(n2*cos(phieta)+n1*cos(phietas)); %коэффициент отражения параллельный плоскости падения
rpe=(n1*cos(phieta)-n2*cos(phietas))/(n1*cos(phieta)+n2*cos(phietas)); %коэффициент отражения перпендикулярный к плоскости падения
rop=(abs(rpa^2)+abs(rpe^2))/2;
end
function rov = PoiskReflVtor(ugalp,phietass,n1,n2) %ugalp, phietass - угол падения, угол преломления
rpe=(n2*cos(ugalp)-n1*cos(phietass))/(n2*cos(ugalp)+n1*cos(phietass)); %коэффициент отражения перпендикулярный к плоскости падения 
rpa=(n1*cos(ugalp)-n2*cos(phietass))/(n1*cos(ugalp)+n2*cos(phietass)); %коэффициент отражения параллельный плоскости падения
rov=(abs(rpe^2)+abs(rpa^2))/2;
end
function [ vy ] = VybABC(fa,fb,fc,a,b,c)
if (fc>fa) 
        if (fc>fb)
            if (fa>fb)
                b=c;
            else
                a=c;
            end
        else
        a=c;    
        end
else
        b=c;
end
vy=[a   b   c];
end
function pk = PoiskUglaMax(phi,teta,ak)
bk=0; nu=0:1e-2:2*pi; n=length(nu);
for k=1:n
    nuk=nu(k); vxyz=[ak*cos(nuk)        ak*sin(nuk)         0];
    vxyz=PovorotyIzIskhKNov(vxyz,teta,phi); 
    vxyz=Pereobozn(vxyz);
    x0(k)=vxyz(1); y0(k)=vxyz(2); vxyz(3)=0; 
    bk(k)=DlinaVectora(vxyz);
end
ma=max(bk);
diap=PoiskDiap(nu,bk,ma,n);
a=diap(1); b=diap(2); c=(a+b)/2;
ep=1e-6; Nit=1e2; k=0; ra=abs(a-b); 
while ((ra>ep) && (k<Nit))
    c=(a+b)/2;
    xya=PovorotyIzIskhKNov([ak*cos(a)   ak*sin(a)    0],teta,phi); xya(3)=0; fa=DlinaVectora(xya);
    xyb=PovorotyIzIskhKNov([ak*cos(b)   ak*sin(b)    0],teta,phi); xyb(3)=0; fb=DlinaVectora(xyb);
    xyc=PovorotyIzIskhKNov([ak*cos(c)   ak*sin(c)    0],teta,phi); xyc(3)=0; fc=DlinaVectora(xyc);
    ro=VybABC(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3);
ra=abs(a-b); k=k+1;
end
pk=c;
end
function [ diap] = PoiskDiap(phi,ar,m,n)
f=1; l=1; p=1;
for k=2:(n-1)
    if (ar(k)==m)
        if (f>0)
        l=k-1;
        p=k+1;
        f=0;
        end
    end
end
ug(1)=phi(l); ug(2)=phi(p);
a=min(ug); b=max(ug);
diap=[a     b];
end
function pn = PoiskNomera(dlvo,dvr)
n=length(dlvo); f=1; q=1;
for k=1:n
    if ((dlvo(k)>dvr) && (f>0))
        f=0; q=k;
    end
end
pn=q;
end
function vykh = OprTochVykhKritSluch(phi,teta,nalu,tpisk,h,r,vyhod) %в плоскости основания
phin=atan(nalu(2)/nalu(1)); phitp=atan(tpisk(2)/tpisk(1)); ep=1e-10;
vyh=vyhod; ro=0;
if ((phi==0) && (teta==0))
vyh=1;
end
if ((phin==0) || (phin==pi) || (phin==2*pi) || (phin==-pi) || (phin==-2*pi) || (phin==pi/2) || (phin==3*pi/2) || (phin==-pi/2) || (phin==-3*pi/2))
    if ((phitp==0) || (phitp==pi) || (phitp==2*pi) || (phitp==-pi) || (phitp==-2*pi) || (phitp==pi/2) || (phitp==3*pi/2) || (phitp==-pi/2) || (phitp==-3*pi/2))
        if (abs(nalu*tpisk')>ep)
            if (abs(tpisk(2))>ep)
                if ((nalu*tpisk')>0)
            dy=r-tpisk(2); dx=tpisk(1); ro=sqrt(dx^2+dy^2);
                else
                dy=r-tpisk(2); dx=tpisk(1); ro=sqrt(dx^2+dy^2);
                end
            else
                if ((nalu*tpisk')>0)
                dx=r-tpisk(1); dy=tpisk(2); ro=sqrt(dx^2+dy^2);
                else
                dx=r-tpisk(1); dy=tpisk(2); ro=sqrt(dx^2+dy^2);
                end
            end
        else
            x=tpisk(1); ro=sqrt(r^2-x^2);
        end
    end
end
if (ro<h)
    vyh=2; %выход из боковой поверхности
    else vyh=1; %выход из основания
end
vykh=vyh;
end
function t = VyvodTrs(dvr,vyfr,vyko,fidTrs)
        dlv=(1e6)*dlinyvoln();
        n=length(dvr);
        if ((vyko==3) && (vyfr==4))
        mkbr=[250.882 249.590 249.213 250.299 249.441 249.365];
        mv=[0.320 0.533 0.849 0.223 0.502 0.797]; 
        tol=[0.76 0.72 0.69 0.73 0.73 0.73]; 
        rokbr=2.75; rov=0.56; n=length(mv);
        for k=1:n
        vkbr(k)=mkbr(k)/(1e3*rokbr); 
        vv(k)=mv(k)/(1e3*rov); 
        xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
        end
        TrRef=SredGrafRass(vyfr,vyko,xv);
        for k=1:length(dvr)
        w=PoiskNomera(dlv,dvr(k));
        TRs(k)=TrRef(w-1)+(TrRef(w)-TrRef(w-1))*(dvr(k)-dlv(w-1))/(dlv(w)-dlv(w-1));
        fprintf(fidTrs,'%0.15f\n',TRs(k));
        end
        end
t=length(TRs);
end