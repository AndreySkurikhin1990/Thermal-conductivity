%��������� ������ ��� ��� ��� ����
function [ t ] = tmp77()
str10='Pryamougolnye_Parallelepipedy_'; 
str1000='itom-1000.txt'; str620='itom-620.txt'; str860='itom-860.txt'; str440='itom-440.txt'; 
str3=str440;
str2='Tr_';  str1=str10; str1=strcat(str1,str2); str1=strcat(str1,str3); fidTr=fopen(str1,'w'); 
str2='dpp_'; str1=str10; str1=strcat(str1,str2); str1=strcat(str1,str3); fiddpp=fopen(str1,'w'); 
str2='Trs_'; str1=str10; str1=strcat(str1,str2); str1=strcat(str1,str3); fidTrs=fopen(str1,'w'); 
str2='dnra_'; str1=str10; str1=strcat(str1,str2); str1=strcat(str1,str3); fiddnra=fopen(str1,'w'); 
str2='razdT_'; str1=str10; str1=strcat(str1,str2); str1=strcat(str1,str3); fidrazdT=fopen(str1,'w'); 
str2='nKBr_'; str1=str10; str1=strcat(str1,str2); str1=strcat(str1,str3); fidn1=fopen(str1,'w'); 
str2='dl_vo_'; str1=str10; str1=strcat(str1,str2); str1=strcat(str1,str3); fiddv=fopen(str1,'w'); 
n440=(1.543-1)*0.6+0.4*(1.56-1)+1; %n620=(1.543-1)*0.5+0.5*(1.56-1)+1; %n860=(1.543-1)*0.25+(1-0.25)*(1.56-1)+1; %n1000=(1.543-1)*0.2+(1-0.2)*(1.56-1)+1;
npp=Kramers_n_Itom_uk()-n440;
na=1.4; ko=27; de=1e-1; ep=de/2;
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
npma=max(np); npmi=min(np); 
for k=1:length(np)
trk=nacha(np(k),dvr(k),npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
tr(k)=trk(1);
end
fclose(fidTr); fclose(fiddpp); fclose(fiddnra); fclose(fidrazdT); fclose(fidTrs); fclose(fidn1); fclose(fiddv);
t=0;
end
function tm = nacha(popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
format long g; %���������� ����������� ��������� ��� �������� KBr + ����������
t=2*erfcinv(1e-4);
vyfr=1; %no - ����� ������� - ������ ����� 1
vyko=1; %nom - ����� ������������
vypl=3; %����� �������� ���������: 1 - 620, 2 - 860, 3 - 1000 
ti=RasFraItom(t,popr,dv,vyfr,vyko,vypl,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
tm=ti;
end
function pk = PoisDn(vyfr,vyko,vypl,xv,vv,mo,si,mi,ma,dpp,dvr,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
TrRef=SredGrafRass(vypl); dlv=dlinavolny(); dvko=dvr;
xve=opredEffPar(xv,vyfr,vyko)
ide=1; n1=OprPokPreKBr(dvr); TR=0; dnre=0; vve=opredEffPar(vv,vyfr,vyko);
dnre=0; q=1; eps=1e-3; dpp=dpp'
while (dvr<=dvko)
    w=PoiskNomera(dlv,dvr); dlvo=dlv(w); chit=1e1; TRs=TrRef(w-1)+(TrRef(w)-TrRef(w-1))*(dvr-dlv(w-1))/(dlv(w)-dlv(w-1));
    TRs=exp(-(1e-6)*TRs*xve)
    Tra=TRs; Trb=1e0; sc=0; razdT=abs(Tra-Trb); TR=(Tra+Trb)/2; dnc=1e2;
while ((razdT>eps) && (sc<chit))
    TR=(Tra+Trb)/2;
    dnb=npma; dna=npmi; dnb=1e3; 
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
    ep=1e-7; Nit=1e2; k=0; ra=abs(dna-dnb); dnc=(dna+dnb)/2; ddnc=dnc; kk=0; nikk=1e2;
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
%������ ������� ������ ���������� ����������� ��������� ��� ���������
%mo - �����. ����-�, si - ���, mi, ma - ���. � ����. ����., no - ����� (���. ��� ����. ���-�)
function mdn = opreKoefProp(vvk,mo,si,mi,ma,dlvo,dn,ide)
n1=OprPokPreKBr(dlvo); n2=n1+dn; rosr=0;
Nr=1e0; raot=13e3; blis=0; ocp=0; sepoo=0; %raot - ������ ��������� � ��������, blis - ����� ������� ���������
for g=1:Nr %���������� �� �����������
suv=0; k=0; nk=1e12;
while ((suv<vvk) && (k<=nk))
hk=mo+si*randn(); hk=prov(hk,mi,ma); %����������� ����� ������������
ak=mo+si*randn(); ak=prov(ak,mi,ma); bk=mo+si*randn(); bk=prov(bk,mi,ma); %����������� ���� ��������� ��� ��
tetani=pi*rand(); %����������� ��������� - ���� ����� ���������� ��������� � ���� Z (�����)
phini=pi*rand(); %���� ����� ��������� ��� �������� �� �������������� ��� � ���� X - ����������� ����������� ���������
upvoz=pi*rand(); %���� �������� ������ ��� Z
suv=suv+hk*ak*bk; %����� ��
    vat=[ak,0,0]; vat=PovorotIzStarkNov(vat,upvoz,tetani,phini); dva=DlinaVectora(vat); %���� ��������� �� ������������ ���������
    vbt=[0,bk,0]; vbt=PovorotIzStarkNov(vbt,upvoz,tetani,phini); dvb=DlinaVectora(vbt); %������ ��������� �� ������������ ���������
    vabt=vat+vbt; dvab=DlinaVectora(vabt); %������� ��������� �� ������������ ���������
    vab=vat-vbt; dvabm=DlinaVectora(vab); %������� ��������� �� ������������ ���������
    vht=[0,0,-hk]; vht=PovorotIzStarkNov(vht,upvoz,tetani,phini); dvh=DlinaVectora(vht); %������ - �������� ���������� ������� ����������� �� ������������ ���������, ��� �������� �������� a=(a,0,0), b=(0,b,0), h=(0,0,-h), XZ - � ����� ����������� ������������ ���������    
    vaht=vat+vht; vbht=vbt+vht; vabht=vat+vbt+vht; uab=acos(abs(vbt*vat'))/dva/dvb; %������ ���� ����� ���������� �������� a � b
    uah=acos(abs(vat*vht'))/dva/dvh; %������ ���� ����� ���������� a � h
    ubh=acos(abs(vbt*vht'))/dvb/dvh; %������ ���� ����� ���������� b � h
    svo=dvb*dva*sin(uab); svb=dva*dvh*sin(uah)+dvb*dvh*sin(ubh); svb=svb+svo; %�������� ��������� � ������� �����������
    mz=VysLuchaPP(svo,svb,tetani,hk,abs(ak),abs(bk),n1,n2,dn,dlvo,raot,upvoz,phini,dva,dvb,dvh,uab,uah,ubh,ide,vat,vbt,vabt,vht); 
    vys=mz(1); vys=proverka(vys); blis=blis+vys; %1 - �������, 0 - ���
    ko=mz(2); ko=proverka(ko); rosr=rosr+ko; sepoo=sepoo+hk*(ak*cos(upvoz)+bk*sin(upvoz)); ocp=ocp+1;
k=k+1;
end
end
w=blis/ocp;
rosr=rosr/ocp;
p=sepoo/Nr/(pi*(raot^2)/4);
mdn=p*w*(1-rosr)+(1-p); %����������� �� ���� ���������
end
function rop = PoiskReflPerv(phi,phis,n1,n2)
rpa=(n2*cos(phi)-n1*cos(phis))/(n2*cos(phi)+n1*cos(phis)); %����������� ��������� ������������ ��������� �������
rpe=(n1*cos(phi)-n2*cos(phis))/(n1*cos(phi)+n2*cos(phis)); %����������� ��������� ���������������� � ��������� �������
rop=(abs(rpa^2)+abs(rpe^2))/2;
end
function rov = PoiskReflVtor(phips,phiss,n1,n2)
rpe=(n2*cos(phips)-n1*cos(phiss))/(n2*cos(phips)+n1*cos(phiss)); %����������� ��������� ���������������� � ��������� ������� 
rpa=(n1*cos(phips)-n2*cos(phiss))/(n1*cos(phips)+n2*cos(phiss)); %����������� ��������� ������������ ��������� �������
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
if (s<svo) %������ � ���������
    phiup=acos(abs(ugpalu(3))); phiup=provUgla(abs(phiup)); %���� (�����������) ������� ���� �� ���������
    mz=PoiskTetaShtrNach(dn,phiup,n1,n2); phiups=mz(1); pvon=mz(3); %phiups - ���� �����������, phiup - ���� �������
    if ((phiups>=0) && (pvon==0))
        phiups=provUgla(abs(phiups));
        rop=PoiskReflPerv(phiup,phiups,n1,n2);
        nalu=[ugpalu(1),ugpalu(2)]; %����������� ���� � ��������� XY (���������)
        topa=[ak*rand(),bk*rand(),0]; %����� ������� �� ��������� XY
        mz=VykhodIzOsnPP(nalu,topa,ak,bk,hk,phiups); v=mz(1); vpp=mz(4); plo=mz(2); ugbe=mz(3); %����, ������ ����� �������� ���, � ��������� ������ ���������
    if (v==2) %1 - ����� ����� ���������, 2 - ����� ���. ���-��
        phiups=pi/2-abs(phiups); %����� ���� �������
        mz=PoiskVykhUglaSSKonech(dn,phiups,n1,n2); pvok=mz(3); 
    if (pvok==0)
    phiupss=pi/2-mz(1); %������� ����� ������� �����������
    rov=PoiskReflVtor(phiups,mz(1),n1,n2);
    po=PopNePop(dlv,abs(phiup-phiupss),razotv); %��������� � �������� ����� �����
    else %��� ������ ���������� ��������� ���� ������ ���� ��� ��������� �� ������� �����������
        nlxyz=n1*[sign(ugpalu(1))*cos(phipaluxy)*sin(phiup),sign(ugpalu(2))*sin(phipaluxy)*sin(phiup),sign(ugpalu(3))*cos(phiup)]; %����������� ���� ��� �������
        nlxyzs=n2*[sign(nlxyz(1))*cos(phipaluxy)*sin(phiups),sign(nlxyz(2))*sin(phipaluxy)*sin(phiups),sign(nlxyz(3))*cos(phiups)]; %����������� ���� ��� �����������
        nnrlo=[nlxyzs(1),nlxyzs(2)]; nnrlo=nnrlo/DlinaVectora(nnrlo); topaxy=[sign(nlxyzs(1))*abs(nnrlo(1)*plo)+topa(1),sign(nlxyzs(2))*abs(nnrlo(2)*plo)+topa(2)];
        switch (vpp)
            case 3
            nlxyzs(1)=-nlxyzs(1); veni=2; %��������� �� �������, YZ
            case 4
            nlxyzs(1)=-nlxyzs(1); veni=2; %��������� �� ������, YZ
            case 1 
            nlxyzs(2)=-nlxyzs(2); veni=1; %��������� �� ��������, XZ
            case 2 
            nlxyzs(2)=-nlxyzs(2); veni=1; %��������� �� �������, XZ
        end
        topaz=DlinaVectora(topaxy)/tan(ugbe); tpl=[topaxy(1),topaxy(2),sign(nlxyz(3))*abs(topaz)];
        mz=VykhodIzOsnPPPVO(ak,bk,hk,nlxyzs,tpl,veni,phiups,vpp,3,1,0); %����, ������ ����� �������� ���
        vkp=mz(5); ugobet=mz(2);
        if ((((vkp==1) || (vkp==2)) && ((vpp==3) || (vpp==4))) || (((vkp==3) || (vkp==4)) && ((vpp==1) || (vpp==2))))
        phiups=pi/2-abs(phiups); %����� ���� ������� ����� ������� ����������� ���������
        end
        mz=PoiskVykhUglaSSKonech(dn,phiups,n1,n2); pvokk=mz(3); phiupss=mz(1);
        if ((pvokk==0) && (ugobet>=0))
        rov=PoiskReflVtor(phiups,phiupss,n1,n2);
        po=PopNePop(dlv,abs(phiup-phiupss),razotv); %��������� � �������� ����� �����
        else
            po=0; rop=0; rov=0;
        end
    end
    else
        mz=PoiskVykhUglaSSKonech(dn,abs(phiups),n1,n2); pvok=mz(3);
        if (pvok==0)
            phiupss=mz(1); phiupss=provUgla(abs(phiupss));
            po=1; rov=PoiskReflVtor(phiups,phiupss,n1,n2); %����� ��������� - �������� ����� ���� ���
        else po=0; rov=0; rop=0;
        end
    end
    else rop=0; rov=0; po=0;
    end
    if (ide==9) %���������� �������
        cem=1e2; mxa=0:(vat(1)/cem):vat(1); mya=0:(vat(2)/cem):vat(2); mxh=0:(vht(1)/cem):vht(1); myh=0:(vht(2)/cem):vht(2); %����� ���������� %ska=0:(ak/cem):ak; skh=0:(hk/cem):hk; skh=-skh; skb=0:(bk/cem):bk; %������ ����������
        pl=plot(mxa+vbt(1),mya+vbt(2),':r',mxh+vbt(1),myh+vbt(2),':g',mxh+vabt(1),myh+vabt(2),':b',mxa+vbt(1)+vht(1),mya+vbt(2)+vht(2),':m'); %����� ���������� %pl=plot(ska,skh,':b',ska,skh,':k'); %������ ����������
        set(pl,'LineWidth',3); hold on; grid on; xlabel({'x'}); ylabel({'y'}); title({'y(x)'}); hold off;
    end
else %��������� � ������� �����������
    vyav=vh*va*sin(uah); vybv=vh*vb*sin(ubh); vypo=vyav+vybv;
    vxy=vypo*rand();
    if (vxy<=vyav) %1 - ������� ���������, XZ
        tpx=ak*rand(); 
        if (phiz>pi/2) 
        tpy=bk; 
        else tpy=0;
        end
        tpz=-hk*rand(); veni=1;
    else %2 - ������ ���������, YZ
        tpx=ak; tpy=bk*rand(); tpz=-hk*rand(); veni=2;
    end %���������� ����� ������� � ������ �����������
       tpl=[tpx,tpy,tpz]; upp=opreUgPadPre(ugpalu,veni,dn,n1,n2); 
       phiup=upp(1); phiups=upp(2); narlu=[upp(3),upp(4)]; %���� ���� ������� � �����������
       if (phiups>=0)
           rop=PoiskReflPerv(phiup,phiups,n1,n2);
           mz=VykhodIzOsnBokPovPP(ak,bk,hk,narlu,tpl,veni,phiups); v=mz(1); vpp=mz(4); plbp=mz(3); ugbe=mz(2); %����, ������ ����� �������� ���, � ��������� ������ ���������
           phiups=provUgla(abs(phiups));
    switch (v)
        case 2 %������� ����������� ��������
                mz=PoiskVykhUglaSSKonech(dn,phiups,n1,n2); pvok=mz(3);
        if (pvok==0)
                phiupss=abs(mz(1)); phiupss=provUgla(phiupss);
                po=1; rov=PoiskReflVtor(phiups,phiupss,n1,n2); %����� ��������� - �������� ����� ���� ���
        else po=0; rov=0;
        end
        case 1 %����� ������� ��� ������ ��������� ��� ������� ����������� �����
                mz=PoiskVykhUglaSSKonech(dn,pi/2-phiups,n1,n2); pvok=mz(3);
                if (pvok==0)
                phiupss=pi/2-abs(mz(1)); phiupss=provUgla(phiupss);
                rov=PoiskReflVtor(pi/2-phiups,pi/2-phiupss,n1,n2);
                po=PopNePop(dlv,abs(phiup-phiupss),razotv); %��������� � �������� ����� �����
                else %��� ������ ���������� ��������� ���� ������ ���� ��� ��������� �� ������� ����������� � �� ���������
                    narlu=ugpalu;
                    if (veni==1) %XZ: ak*hk
                    phixyz=abs(narlu(3)/narlu(1)); phixyz=atan(phixyz); phixyz=provUgla(phixyz);
                    nlxyz=n1*[sign(ugpalu(1))*cos(phixyz)*sin(phiup),sign(ugpalu(2))*cos(phiup),sign(ugpalu(3))*sin(phixyz)*sin(phiup)]; %����������� ���� ��� �������
                    nlxyzs=n2*[sign(nlxyz(1))*cos(phixyz)*sin(phiups),sign(nlxyz(2))*cos(phiups),sign(nlxyz(3))*sin(phixyz)*sin(phiups)]; %����������� ���� ��� �����������
                    nnrlo=[nlxyzs(1),nlxyzs(3)]; nnrlo=nnrlo/DlinaVectora(nnrlo); 
                    topaxz=[sign(nlxyzs(1))*abs(nnrlo(1)*plbp)+tpl(1),sign(nlxyzs(3))*abs(nnrlo(2)*plbp)+tpl(3)]; %����������� ���� �� ������� �����������
                    switch (vpp)
            case 3
            nlxyzs(1)=-nlxyzs(1); venin=2; %��������� �� �������, YZ
            case 4
            nlxyzs(1)=-nlxyzs(1); venin=2; %��������� �� ������, YZ
            case 1 
            nlxyzs(3)=-nlxyzs(3); venin=1; %��������� �� ��������, XY
            case 2 
            nlxyzs(3)=-nlxyzs(3); venin=1; %��������� �� �������, XY
                    end
                    topay=DlinaVectora(topaxz)/tan(ugbe); tpl=[topaxz(1),sign(nlxyzs(2))*abs(topay),topaxz(2)];
                    mz=VykhodIzOsnPPPVO(ak,hk,bk,nlxyzs,tpl,venin,phiups,vpp,2,0,1); %����, ������ ����� �������� ���
                    elseif (veni==2) %YZ: bk*hk
                        phixyz=abs(narlu(3)/narlu(2)); phixyz=atan(phixyz); phixyz=provUgla(phixyz);
                        nlxyz=n1*[sign(ugpalu(1))*cos(phiup),sign(ugpalu(2))*cos(phixyz)*sin(phiup),sign(ugpalu(3))*sin(phixyz)*sin(phiup)]; %����������� ���� ��� �������
                        nlxyzs=n2*[sign(nlxyz(1))*cos(phiups),sign(nlxyz(2))*cos(phixyz)*sin(phiups),sign(nlxyz(3))*sin(phixyz)*sin(phiups)]; %����������� ���� ��� �����������
                        nnrlo=[nlxyzs(2),nlxyzs(3)]; nnrlo=nnrlo/DlinaVectora(nnrlo); topayz=[sign(nlxyzs(2))*abs(nnrlo(1)*plbp)+tpl(2),sign(nlxyzs(3))*abs(nnrlo(2)*plbp)+tpl(3)];
                        switch (vpp)
            case 3
            nlxyzs(2)=-nlxyzs(2); venin=2; %��������� �� �������, XZ
            case 4
            nlxyzs(2)=-nlxyzs(2); venin=2; %��������� �� ������, XZ
            case 1 
            nlxyzs(3)=-nlxyzs(3); venin=1; %��������� �� ��������, XY
            case 2 
            nlxyzs(3)=-nlxyzs(3); venin=1; %��������� �� �������, XY
                        end
                        topax=DlinaVectora(topayz)/tan(ugbe); tpl=[topax,topayz(1),topayz(2)];
                        mz=VykhodIzOsnPPPVO(bk,hk,ak,nlxyzs,tpl,venin,phiups,vpp,1,0,2); %����, ������ ����� �������� ���
                    end
                    vkp=mz(5); ugobet=mz(2);
        if ((((vkp==1) || (vkp==2)) && ((vpp==3) || (vpp==4))) || (((vkp==3) || (vkp==4)) && ((vpp==1) || (vpp==2))))
        phiups=pi/2-abs(phiups); %����� ���� ������� ����� ������� ����������� ���������
        end
        mz=PoiskVykhUglaSSKonech(dn,phiups,n1,n2); pvokk=mz(3); phiupss=mz(1);
        if ((pvokk==0) && (ugobet>=0))
        rov=PoiskReflVtor(phiups,phiupss,n1,n2);
        po=PopNePop(dlv,abs(phiup-phiupss),razotv); %��������� � �������� ����� �����
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
end %����� ������� ���� �� ������� �����������
vy=[po,rop*rov];
end
function pr10 = PraCha10(phi,phis,n1,n2)
pr10=n1*sin(phi)-n2*sin(phis);
end
function pr20 = PraCha20(phis,phiss,n1,n2)
pr20=n2*sin(phis)-n1*sin(phiss);
end
%���������� ����������� ������� �����
function opn = OprPokPreKBr(la)
lam=[1,2,10,20,25]; pp=[1.487,1.48,1.52,1.48,1.453]; 
ko=polyfit(lam,pp,length(pp)-1); opn=polyval(ko,la);
end
function [ vy ] = VykhodIzOsnPPPVO(ak,bk,hk,narlu,tpl,veni,phiups,vpp,pk,os,vn) %����, ������ ����� �������� ���
vyhm=[0,0,0,0,0];
if (os==1) %��������� - 1, ������� ����������� - 0
    vyhm=VybVykhOsnPVO(ak,bk,hk,tpl,narlu,phiups,veni,vpp,pk); %������� ������� �� �������� b
elseif (os==0)
if (vn==1) %ak*hk 
    vyhm=VybVykhOsnPVO(ak,hk,bk,tpl,narlu,phiups,veni,vpp,pk); %������� ������� �� �������� b
elseif (vn==2) %bk*hk
    vyhm=VybVykhOsnPVO(bk,hk,ak,tpl,narlu,phiups,veni,vpp,pk); %������� ������� �� �������� a
end
end
    vy=vyhm; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
end
function vy = VybVykhOsnPVO(ak,bk,hk,tpl,nlxyz,phiups,veni,vpp,peko)
vyh=1;
uvav=abs(bk-tpl(2)); uvav=atan(abs(ak/uvav)); %����, ��� ������� ����� ������� ����� ������� a
uvnk=abs(ak/tpl(2)); uvnk=atan(uvnk); %���� ��������� ������ ���������, 2
uvnk1=abs(bk/tpl(1)); uvnk1=atan(uvnk1); %���� ��������� ������ ���������, 1
uvav1=abs(bk/(ak-tpl(1))); uvav1=atan(uvav1);
unl=abs(nlxyz(2)/nlxyz(1)); unl=atan(unl); %����������� ��������������� ����
mz=opredPutLuchaPVO(nlxyz,unl,uvav,uvnk,uvav1,uvnk1,tpl,ak,bk,veni,vpp); 
pulu=mz(1); vnp=mz(2); vkp=mz(3); ugolbeta=-1;
if (pulu>0)
ugolbeta=pulu/abs(hk-tpl(peko)); ugolbeta=atan(ugolbeta); ugolbeta=provUgla(ugolbeta);
if (ugolbeta>phiups)
    vyh=1;
else vyh=2; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
end
else
    vyh=0;
end
vy=[vyh,ugolbeta,pulu,vnp,vkp];
end
function [ pl ] = opredPutLuchaPVO(nlxyz,unl,uvav,uvnk,uvav1,uvnk1,tpl,ak,bk,veni,vpp)
pulu=-1; vnp=0; vkp=0;
if ((nlxyz(1)>=0) && (nlxyz(2)>=0))
    if ((vpp==4) && (veni==2)) %������� �������
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
    if ((veni==2) && (vpp==3)) %������� �������
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
    if ((veni==2) && (vpp==3)) %������� �������
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
    if ((veni==2) && (vpp==4)) %������� �������
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
%����� ���� �����������, phi>0
function [ pesdp ] = PoiskTetaShtrNach(dnr,phipad,n1,n2)
phipre=0; pvo=0; a=0; b=pi/2; ep=1e-6; Nit=1e3; k=0; ra=abs(a-b); dephi=0; %������� ���� �����������
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
%����� ���� ������: v=1 - ����� ����� ���������, v=2 - ����� ���. ���-��
function [ pesdp ] = PoiskVykhUglaSSKonech(dnr,phi,n1,n2)
phipre=0; pvo=0; a=0; b=pi/2; ep=1e-7; Nit=1e3; k=0; ra=abs(a-b); dephi=0; %������� ���� �����������
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
%------------------���� �� �������--------------
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
function [ alsr ] = SredGrafRass(vypl)
alsre=0;
if (vypl == 1)
fileID = fopen('Koefficient_pogloscheniya_itom620.txt','r'); 
alsre=fscanf(fileID,'%f'); 
fclose(fileID);
elseif (vypl == 2)
fileID = fopen('Koefficient_pogloscheniya_itom860.txt','r'); 
alsre=fscanf(fileID,'%f'); 
fclose(fileID);
elseif (vypl == 3)
fileID = fopen('Koefficient_pogloscheniya_itom1000.txt','r'); 
alsre=fscanf(fileID,'%f'); 
fclose(fileID);
end
alsr=alsre;
end
%di - ��� ������ ���, ide - ������� (1) ��� �������. ���-� (2)
function ras60 = RasFraItom(di,popr,dv,vyfr,vyko,vypl,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv)
mkbr=[230.078,231.006]; mkbr=mkbr(vyko); %� ��
mv=[0.95,1.19]; mv=mv(vyko); %� ��
tol=[0.64,0.64]; tol=tol(vyko); % � ��
switch (vypl)
    case (1)
        rov=0.62;
    case (2)
        rov=0.86;
    case (3)
        rov=1;
end
rokbr=2.75;
vkbr=mkbr/(1e3*rokbr); %� ��3
vv=mv/(1e3*rov); %� ��3
xv=(vv/(vv+vkbr))*tol*1e3; %� ���
%vkbr=vkbr*((1e4)^3); %� ���3
%vv=vv*((1e4)^3); %� ���3
if (vyfr==1)
mini=0; maxi=63; 
elseif (vyfr==2)
mini=63; maxi=1e2; 
elseif (vyfr==3)
mini=1e2; maxi=15e1; 
else
mini=15e1; maxi=2e2; 
end
mo=(maxi+mini)/2; %��
si=abs(mini-maxi)/di; %���
depp=PoisDn(vyfr,vyko,vypl,xv,vv,mo,si,mini,maxi,popr,dv,npmi,npma,fidTr,fiddpp,fidTrs,fiddnra,fidrazdT,fidn1,fiddv);
ras60=depp;
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
no=vyko; no=vyfr; no=1;
oep=mapa(no);
end
function [ vy ] = VykhodIzOsnPP(nlxy,tpl,ak,bk,hk,phiups)
vyh=1;
uvbl=bk-tpl(2); uvbl=atan(abs(uvbl/tpl(1))); %����, ��� ������� ����� ����� ����� ������� b
uvnk=tpl(2)/tpl(1); uvnk=atan(abs(uvnk)); %���� ��������� ������ ���������
uvav=bk-tpl(2); uvav=abs(uvav/(ak-tpl(1))); uvav=atan(uvav); %����, ��� ������� ����� ������� a �������
uvbp=ak-tpl(1); uvbp=abs(tpl(2)/uvbp); uvbp=atan(uvbp); %���� ��������� ������ ������� b
unl=abs(nlxy(2)/nlxy(1)); unl=atan(unl); %����������� ��������������� ����
mz=opredPutLucha(nlxy,unl,uvbl,uvnk,uvav,uvbp,tpl,ak,bk); pulu=mz(1); panapo=mz(2);
if (pulu>=0)
ugolbeta=pulu/hk; ugolbeta=atan(ugolbeta);
if (ugolbeta>phiups)
    vyh=2;
else vyh=1; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
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
x=xs; y=ys*cos(phi)+zs*sin(phi); z=-ys*sin(phi)+zs*cos(phi); %������ ������� � ��������� YZ ������ ��� X
koor=[x,y,z];
end
function [ koor ] = PovorotObratnonaUgolTeta(koors,teta)
xs=koors(1); ys=koors(2); zs=koors(3);
x=xs*sin(teta)+zs*cos(teta); y=ys; z=-xs*cos(teta)+zs*sin(teta); %������ ������� � ��������� XZ ������ ��� Y
koor=[x,y,z];
end
function [ nk ] = PovorotObratnonaVokrugOZ(koors,phiz)
xs=koors(1); ys=koors(2); zs=koors(3);
x=xs*cos(phiz)-ys*sin(phiz); y=xs*sin(phiz)+ys*cos(phiz); z=zs; %������ ������� � ��������� XY ������ ��� Z
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
function [ upapr ] = opreUgPadPre(ugpalu,veni,dn,n1,n2) %���� (�����������) ������� ����
switch (veni)
    case 2 %2 - ������ ��������� YZ - b, (-1,0,0); 1 - ������� XZ - a, (0,-1,0)
        phiup=acos(ugpalu(1)); phiup=provUgla(phiup); %�������� � ��������� YZ
        nalu=[ugpalu(2),ugpalu(3)]; pxyz=1; %����������� ���� �� ������� ����������� YZ
    case 1 
        phiup=acos(ugpalu(2)); phiup=provUgla(phiup); %�������� � ��������� XZ
        nalu=[ugpalu(1),ugpalu(3)]; pxyz=2; %����������� ���� �� ������� ����������� XZ
end
    mz=PoiskTetaShtrNach(dn,phiup,n1,n2); pvo=mz(3);
    if (pvo==0)
    phiups=mz(1); phiups=provUgla(phiups); %phiups - ���� �����������, phiup - ���� �������
    else
        phiups=-1; %������ ���������� ���������
    end
    upapr=[phiup,phiups,nalu(1),nalu(2),pxyz]; %����� ������� �� ��������� XY � ����������� ����
end
function [ vy ] = VykhodIzOsnBokPovPP(ak,bk,hk,narlu,tpl,veni,phiups) %tpl - ���������� ����� ������� � �������� ����������� XYZ
vyh=1; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
if (veni==2) %� ������ ��������� - YZ 
    vyhm=VybVykhOsnBokPov(bk,hk,ak,tpl,narlu,phiups); 
else %� ������� ��������� - XZ 
    vyhm=VybVykhOsnBokPov(ak,hk,bk,tpl,narlu,phiups); 
end
    vy=vyhm;
end
function vy = VybVykhOsnBokPov(ak,hk,bk,tpl,nlxy,phiups)
vyh=1;
uvhl=abs(hk-tpl(2)); uvhl=atan(abs(uvhl/tpl(1))); %����, ��� ������� ����� ����� ����� ������� h
uvnk=abs(tpl(2)/tpl(1)); uvnk=atan(uvnk); %���� ��������� ������ ���������
uvav=abs(hk-tpl(2)); uvav=uvav/abs(ak-tpl(1)); uvav=atan(uvav); %����, ��� ������� ����� ������� a (b) �������
uvhp=abs(ak-tpl(1)); uvhp=abs(tpl(2)/uvhp); uvhp=atan(uvhp); %���� ��������� ������ ������� h
unl=abs(nlxy(2)/nlxy(1)); unl=atan(unl); %����������� ��������������� ����
mz=opredPutLucha(nlxy,unl,uvhl,uvnk,uvav,uvhp,tpl,ak,hk); 
pulu=mz(1); kpp=mz(2); ugolbeta=0;
if (pulu>0)
ugolbeta=pulu/bk; ugolbeta=atan(ugolbeta); ugolbeta=provUgla(ugolbeta);
if (ugolbeta>phiups)
    vyh=1;
else vyh=2; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
end
else
    vyh=0;
end
vy=[vyh,ugolbeta,pulu,kpp];
end
function [ dv ] = dlinavolny()
dv=(1e6)*RasshDiapDlinVoln();
end