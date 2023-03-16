%����������� �� ��� ��
function [ t ] = tmp39()
t=nacha();
end

function [ tm ] = nacha()
format long g; %���������� ����������� ��������� ��� �������� KBr+����������
t=2*PoisKorn();
nom=2; no=4; %no - ����� �������
switch (no) %nom = 1 - ��� ���������, 2 - ��� ������������� ����������������, 3 - ��� �����
    case 1
ti=RasFra60(t,nom);
    case 2
ti=RasFra60_100(t,nom);
    case 3
ti=RasFra100_150(t,nom);
    case 4
ti=RasFra150_200(t,nom);
end
tm=ti;
end
function [ pk ] = PoisDn(vyfr,vyko,xv,vv,mo,si,mi,ma,ide)
TrRef=SredGrafRass(vyfr,vyko,xv); dlv=1e6*dlivoln(); numele=1e0; na=1e3;
xve=opredEffPar(xv,vyfr,vyko); vve=opredEffPar(vv,vyfr,vyko); dnre=0; q=1;
for w=na:na+numele
    dlvo=dlv(w); TR=TrRef(w); nkbr=OprPokPreKBr(dlvo);
    dna=-1; dnb=1e3; ep=1e-6; Nit=1e3; k=0; ra=abs(dna-dnb);
while ((ra>ep) && (k<Nit))
    dnc=(dna+dnb)/2;
    ppa=nkbr+dna; Ra=abs(ppa^2-nkbr^2)/(ppa^2+nkbr^2); 
    fa=opreKoefProp(vve,mo,si,mi,ma,dlvo,ide,dna); fa=fa*((1-Ra)^2)-TR;
    ppb=nkbr+dnb; Rb=abs(ppb^2-nkbr^2)/(ppb^2+nkbr^2); 
    fb=opreKoefProp(vve,mo,si,mi,ma,dlvo,ide,dnb); fb=fb*((1-Rb)^2)-TR;
    ppc=nkbr+dnc; Rc=abs(ppc^2-nkbr^2)/(ppc^2+nkbr^2); 
    fc=opreKoefProp(vve,mo,si,mi,ma,dlvo,ide,dnc); fc=fc*((1-Rc)^2)-TR;
    ro=vybPoKo(fa,fb,fc,dna,dnb,dnc); dna=ro(1); dnb=ro(2); dnc=ro(3); 
    k=k+1;
    ra=abs(dna-dnb);
end
dnre(q)=dnc; q=q+1;
end
pk=dnre;
end
%������ ������� ������ ���������� ����������� ��������� ��� ���������
%mo - �����. ����-�, si - ���, mi, ma - ���. � ����. ����., no - ����� (���. ��� ����. ���-�)
function mdn = opreKoefProp(vvk,mo,si,mi,ma,dlvo,no,dn)
Nr=5e2; raot=13e3; blis=0; ocp=0; sepoo=0; %raot - ������ ���������, blis - ����� ������� ���������
for g=1:Nr %���������� �� �����������
suv=0; k=1; nk=1e9;
while ((suv<vvk) && (k<nk))
hk=mo+si*randn(); hk=prov(hk,mi,ma); %����������� ����� ������������
tetani=pi*rand(); %����������� ��������� - ���� ����� ���������� ��������� � ���� Z (�����)
phini=-pi/2+pi*rand(); %���� ����� ��������� ��� �������� �� �������������� ��� � ���� X - ����������� ����������� ���������
phiv=abs(acos(cos(phini)*sin(tetani))); %���� ������� - ����� ������������ ��������� � ���� ��������
ocp=ocp+1;
if (tetani<pi/2)
    nuk=pi/2-tetani; %hkp=hk*cos(nuk); ������������� �� �������������� ���
    else
    nuk=pi-tetani; nuk=pi/2-nuk; %������ ���, ����� nu ���� ������ pi/2
end
    uzvp=acos(sin(tetani)*cos(phini)); %���� ��� Z � ������������ ����������
    if (uzvp>pi/2)
        uzvp=pi-uzvp;
    end
    switch (no)
case 1 %�������� - �������� ���������� �������� ������ ��� Z
    dk=mo+si*randn(); dk=prov(dk,mi,ma); rk=dk/2; %����������� �������� ��������
    suv=suv+pi*hk*(rk^2); rkp=rk*abs(cos(tetani)); 
    rkpv=rk*cos(uzvp); hkpv=hk*sin(uzvp); %������������� �� ������������ ���, ���������������� ����������� ���������
    se=pi*rkpv*rk; sp=dk*hkpv+se/2; 
    sepo=se+sp; sepoo=sepoo+sepo; %������� ���������
    vys=VysLucha(se,sepo,tetani,nuk,hk,rk,rkp,dn,dlvo,raot,phini,phiv,hkpv); blis=blis+vys; %1 - �������, 0 - ���
case 2 %������������� ���������������
    ak=mo+si*randn(); ak=prov(ak,mi,ma); bk=mo+si*randn(); bk=prov(bk,mi,ma); %����������� ���� ��������� ��� ��
    suv=suv+hk*ak*bk; %����� ��
    upvoz=pi*(-1+2*rand());
    vat=[ak*cos(upvoz)*sin(tetani) -ak*cos(upvoz)*cos(tetani)*sin(phini)-ak*sin(upvoz)*cos(phini) -ak*cos(upvoz)*cos(tetani)*cos(phini)-ak*sin(phini)*sin(upvoz)];
    vbt=[bk*sin(upvoz)*sin(tetani) bk*cos(upvoz)*cos(phini)-bk*sin(upvoz)*cos(tetani)*sin(phini) bk*sin(upvoz)*cos(tetani)*cos(phini)+bk*cos(upvoz)*sin(phini)];
    vabt=[(ak*cos(upvoz)+bk*sin(upvoz))*sin(tetani) (-ak*sin(upvoz)+bk*cos(upvoz))*cos(phini)-(ak*cos(upvoz)+bk*cos(upvoz))*sin(phini)*cos(tetani)  (-ak*sin(upvoz)+bk*cos(upvoz))*sin(phini)+(ak*cos(upvoz)+bk*cos(upvoz))*cos(phini)*cos(tetani)];
    vht=[hk*cos(tetani) hk*sin(tetani)*sin(phini) -hk*sin(tetani)*cos(phini)]; %��� �������� �������� a=(a,0,0), b=(0,b,0), h=(0,0,-h), XZ - � ����� ����������� ������������ ���������
    vaht=[ak*cos(upvoz)*sin(tetani)+hk*cos(tetani)  (-ak*cos(upvoz)*cos(tetani)+hk*sin(tetani))*sin(phini)-ak*sin(upvoz)*cos(phini) (ak*cos(upvoz)*cos(tetani)-hk*sin(tetani))*cos(phini)-ak*sin(upvoz)*sin(phini)];
    vbht=[bk*sin(upvoz)*sin(tetani)+hk*cos(tetani)  (-bk*sin(upvoz)*cos(tetani)+hk*sin(tetani))*sin(phini)+bk*cos(upvoz)*cos(phini)    (bk*sin(upvoz)*cos(tetani)-hk*sin(tetani))*cos(phini)+bk*cos(upvoz)*sin(phini)];
    vabht=[(ak*cos(upvoz)+bk*sin(upvoz))*sin(tetani)+hk*cos(tetani) (hk*sin(tetani)-(ak*cos(upvoz)+bk*sin(upvoz))*cos(tetani))*sin(phini)+(-ak*sin(upvoz)+bk*cos(upvoz))*cos(phini) ((ak*cos(upvoz)+bk*sin(upvoz))*cos(tetani)-hk*sin(tetani))*cos(phini)+(-ak*sin(upvoz)+bk*sin(upvoz))*sin(phini)];
    vab=[vat(1)-vbt(1) vat(2)-vbt(2) vat(3)-vbt(3)]; dvab=sqrt((vab(1))^2+(vab(2))^2+(vab(3))^2); %��������� �� ������������ ���������
    vaab=[vat(1)-vabt(1) vat(2)-vabt(2) vat(3)-vabt(3)]; dvaab=sqrt((vaab(1))^2+(vaab(2))^2+(vaab(3))^2); %���� ��������� �� ������������ ���������
    vbab=[vbt(1)-vabt(1) vbt(2)-vabt(2) vbt(3)-vabt(3)]; dvbab=sqrt((vbab(1))^2+(vbab(2))^2+(vbab(3))^2); %������ ��������� �� ������������ ���������
    vbbh=[vbt(1)-vbht(1) vbt(2)-vbht(2) vbt(3)-vbht(3)]; dvbbh=sqrt((vbbh(1))^2+(vbbh(2))^2+(vbbh(3))^2); %������ - ��� �������� ������� ������� ����������� �� ������������ ���������
    vaah=[vat(1)-vaht(1) vat(2)-vaht(2) vat(3)-vaht(3)]; dvaah=sqrt((vaah(1))^2+(vaah(2))^2+(vaah(3))^2); %������ - ��� �������� ������ ������� ����������� �� ������������ ���������
    vababh=[vabt(1)-vabht(1) vabt(2)-vabht(2) vabt(3)-vabht(3)]; dvababh=sqrt((vababh(1))^2+(vababh(2))^2+(vababh(3))^2); %������� ������ - ��� �������� ������ ������� ����������� �� ������������ ���������
    uab=acos(abs(vaab*vbab')/dvaab/dvbab); %���� ����� ���������� �������� a � b
    uah=acos(abs(vaab*vaah')/dvaab/dvaah); %���� ����� ���������� a � h
    ubh=acos(abs(vbab*vbbh')/dvbab/dvbbh); %���� ����� ���������� b � h
    svo=dvaab*dvbab*sin(uab); svb=dvaah*dvaab*sin(uah)+dvbab*dvbbh*sin(ubh); svb=svb+svo; %�������� ��������� � ������� �����������
    vys=VysLuchaPP(svo,svb,tetani,hk,ak,bk,dn,dlvo,raot,upvoz,phini,dvaab,dvbab,dvbbh,dvaah,dvababh,uab,uah,ubh); blis=blis+vys; %1 - �������, 0 - ���
        case 3 %����
            ds=mo+si*randn(); ds=prov(ds,mi,ma); rs=ds/2; %������� ����
            phitp=pi*rand(); 
            if (phitp>pi/2)
            phip=pi-phitp; phip=pi/2-phip; 
            else
                phip=pi/2-phitp; 
            end %����� ������� - ������� ����
            phis=PoiskTetaShtrNach(phip,dn); phiss=PoiskTetaShtrNach(phis,dn); 
            vys=PopNePop(dlvo,abs(phip-phiss),raot); suv=suv+4*pi*(rs^3)/3; blis=blis+vys; sepoo=sepoo+pi*(rs^2);
    end
k=k+1;
end
end
w=blis/ocp;
p=sepoo/Nr/(pi*(raot^2)/4);
mdn=p*w+(1-p); %����������� �� ���� ���������
end
function vy = VysLuchaPP(svo,svb,teta,hk,ak,bk,dn,dlv,razotv,phiz,phi,va,vb,dvbbh,dvaah,dvababh,uab,uah,ubh)
s=svb*rand(); po=0;
if (s<svo)
    bv=vb*sin(uab); ztp=bv*rand(); xtp=va*rand(); phia=atan(ztp/xtp); ytp=0;
    if (phia>uab) 
        xtp=xtp+va;
    end %����� ������� ���� �� ���������
    kotp=[xtp ytp ztp]; mxyz=Povoroty(phi,teta,phiz,kotp); xtp3=mxyz(1); ytp3=mxyz(2); ztp3=mxyz(3); ak=abs(ak); bk=abs(bk);
    if ((abs(xtp3)<=ak) && (abs(ytp3)<=bk))
    ugpalu=[-sin(phi)*cos(teta)*cos(phiz)-cos(phi)*sin(phiz)    -sin(phi)*cos(teta)*sin(phiz)+cos(phi)*cos(phiz)    -sin(phi)*sin(teta)]; %���� (�����������) ������� ����
    phiup=acos(abs(ugpalu(3)));
    phiups=PoiskTetaShtrNach(phiup,dn);  %phiups - ���� �����������, phiup - ���� �������
    nalu=[ugpalu(1) ugpalu(2)]; topa=[abs(xtp3) abs(ytp3)]; %����� ������� �� ��������� XY � ����������� ����
    v=VykhodIzOsnPP(nalu,topa,ak,bk,hk); %����, ������ ����� �������� ���, � ��������� ������ ���������
    if (v==2)
    phiupss=PoiskVykhUglaSSKonech(dn,abs(phiups));  %������� ����� ������� �����������
    po=PopNePop(dlv,abs(phiup-phiupss),razotv); %��������� � �������� ����� �����
    else po=1;  %����� ��������� - �������� ����� ���� ���
    end
    end
else
    hv=(dvbbh+dvaah+dvababh)/6+hk/2; va=(va+ak)/2; vyav=va*sin(uah); vb=(vb+bk)/2; vybv=vb*sin(ubh); vypo=vyav+vybv;
    vyz=vypo*rand(); vyx=hv*rand(); uxz=atan(vyz/vyx); vyy=0;
    if (vyz<vybv)
        veni=1; %1 - ������ ���������, YZ
        if (uxz>ubh)
            vyx=vyx+hv;
        end
    else 
        veni=2; %2 - ������� ���������, XZ
        if (uxz>uah)
            vyx=vyx+hv;
        end
    end %����� ������� ���� �� ���������
    kotp=[vyx vyy vyz]; mxyz=Povoroty(phi,teta,phiz,kotp); xtp3=mxyz(1); ytp3=mxyz(2); ztp3=mxyz(3); ak=abs(ak); bk=abs(bk); %kotp - ���������� ����� ������� � �������� �����������, mxyz - � ���������
   if ((abs(xtp3)<=ak) && (abs(ytp3)<=bk))
       upp=opreUgPadPre(phi,teta,phiz,veni,dn); phiup=upp(1); phiups=upp(2); narlu=[upp(3) upp(4)]; %���� ���� ������� � �����������
       v=VykhodIzOsnBokPovPP(ak,bk,hk,narlu,mxyz,veni); %����, ������ ����� �������� ���, � ��������� ������ ���������
    switch (v)
        case 2 %������� ����������� ��������
            phiupss=phiup; dephi=abs(phiup-phiupss); po=1;  %����� ��������� - �������� ����� ���� ���
        case 1 %����� ������� � ������ ���������
                phiups=acos(narlu(2));
                phiupss=PoiskTetaShtrNach(phiups,dn);
                switch (veni)
                    case 1 %������ ��������� ��������� YZ, ���������� 2, 3
                        naralun=[cos(phiup) narlu(1)    narlu(2)]; dvnaralun=sqrt((naralun(1))^2+(naralun(2))^2+(naralun(3))^2);
                        naraluk=[cos(phiups)    narlu(1)    cos(phiupss)]; dvnaraluk=sqrt((naraluk(1))^2+(naraluk(2))^2+(naraluk(3))^2);
                    case 2 %������� ��������� ��������� XZ, ���������� 1, 3
                        naralun=[narlu(1)   cos(phiup)  narlu(2)]; dvnaralun=sqrt((naralun(1))^2+(naralun(2))^2+(naralun(3))^2);
                        naraluk=[narlu(1)   cos(phiups)    cos(phiupss)]; dvnaraluk=sqrt((naraluk(1))^2+(naraluk(2))^2+(naraluk(3))^2);
                end
                dephi=acos((naralun*naraluk')/dvnaraluk/dvnaralun); po=PopNePop(dlv,dephi,razotv); %��������� � �������� ����� ����� 
        case 3 %������� ����������� �����
                switch (veni)
                    case 1 %������ ��������� ��������� YZ, ���������� 2, 3
                        phiupsd=acos(narlu(1)); phiupss=PoiskTetaShtrNach(phiupsd,dn); %����� ����� ������ ��������� - XZ
                        naralun=[cos(phiup) narlu(1)    narlu(2)]; dvnaralun=sqrt((naralun(1))^2+(naralun(2))^2+(naralun(3))^2);
                        naraluk=[cos(phiups)	cos(phiupss)	narlu(2)]; dvnaraluk=sqrt((naraluk(1))^2+(naraluk(2))^2+(naraluk(3))^2);
                    case 2 %������� ��������� ��������� XZ, ���������� 1, 3
                        phiupsd=acos(narlu(1)); phiupss=PoiskTetaShtrNach(phiupsd,dn); %����� ����� ������ ��������� - YZ
                        naralun=[narlu(1)   cos(phiup)  narlu(2)]; dvnaralun=sqrt((naralun(1))^2+(naralun(2))^2+(naralun(3))^2);
                        naraluk=[cos(phiupss)   cos(phiups)    narlu(2)]; dvnaraluk=sqrt((naraluk(1))^2+(naraluk(2))^2+(naraluk(3))^2);
                end %������� ����� ������� �����������
                dephi=acos((naralun*naraluk')/dvnaraluk/dvnaralun); po=PopNePop(dlv,dephi,razotv); %��������� � �������� ����� ����� 
    end
   end
end %����� ������� ���� �� ������� �����������
vy=po;
end
%������� ���� �� �������� - � ��������� ���� �������
function vy = VysLucha(se,spo,teta,nu,hk,rk,rkp,dn,dlv,razotv,phip,phiup,hkpv)
s=spo*rand(); po=0;
if ((teta==pi/2) && (phiup==0))
    po=1;
end
switch (abs(phiup))
    case pi/2
    po=2;
    case 0
    po=1;
    otherwise
if (s<se)
    p=1; zl=0; rlx=rkp*(-1+2*rand()); phitp=poisphitpn(rlx/rkp); %����� �������������� �� �������������� ���
    rly=rlx*tan(phitp); ksi=rlx*cos(phip)+rly*sin(phip); eta=-rlx*sin(phip)+rly*cos(phip); %rlx, rly - ����� �������, ksi, eta - ������� �� phip, phitp - ����� �������
    ksi=ksi/cos(teta); phitp=poisphitp(eta,ksi); %phitp, ksi, eta ���������� ����� ������� ���� �� ���������
    rl=prover(ksi,eta,phitp); %������ �� ���������, X - ������������� ���
    phin=-atan(tan(phip)*cos(teta)); %����������� ��������� �� ��������� ��������
    phiups=PoiskTetaShtrNach(phiup,dn);  %phiups - ���� �����������, phiup - ���� �������
    v=VykhodIzCylOsn(dn,phitp,phin,phiup,rl,rk,hk,phiups); %����, ������ ����� �������� ���, � ��������� ������ ���������
    if (v>1)
    phiupss=PoiskVykhUglaSSKonech(dn,abs(phiups));  %������� ����� ������� �����������
    po=PopNePop(dlv,abs(phiup-phiupss),razotv); %��������� � �������� ����� �����
    else po=1;  %����� ��������� - �������� ����� ���� ���
    end
else p=2; rl=rk; zl=hk*rand(); %������ �� ���. ���-�� (� ����������)
    if (teta>pi/2)
    phitp=pi*(-1+2*rand())/2; %����� ������� ���� �� ������� ����� �� ��������� � OZ
    else
    phitp=pi*(1+2*rand())/2; %����� ������� ���� �� ������ �����, OZ �����
    end
    if (((phitp==0) || (phitp==pi)) && ((phip==0) || (phip==pi/2) || (phip==-pi/2)))
    nus=PoiskTetaShtrNach(abs(nu),dn);
    else
        zl=rk*cos(phitp); yl=rk*sin(phitp); xl=hkpv*rand(); xl=xl/sin(nu); %������� ����� ������� ���� � ����� �����������
        xnl=1; ynl=0; znl=0; %������� ����������� ����
        xn=xl*cos(phip)-yl*sin(phip); yn=xl*sin(phip)+yl*cos(phip); zn=zl; 
        xls=xnl*cos(phip)-ynl*sin(phip); yls=xnl*sin(phip)+ynl*cos(phip); zls=znl;
        xnn=xn*cos(teta)+zn*sin(teta); znn=-xn*sin(teta)+zn*cos(teta); ynn=yn;
        xlss=xls*cos(teta)-yls*sin(teta); ylss=xls*sin(teta)+yls*cos(teta); zlss=zls;
        rlni=sqrt(xlss^2+ylss^2+zlss^2); rltp=sqrt(xnn^2+ynn^2+znn^2); phipn=acos((xnn*xlss+ynn*ylss+znn*zlss)/rlni/rltp); %���� ������� � ���������, ������������ ���������
        nus=PoiskTetaShtrNach(abs(phipn),dn); %���� ����������� � ��������� �������, ����������� ���������
        phipz=acos(zlss/rlni); %���� ����� ������������ ��������� � ���������� � ����������� ���������, ������������ ��� Z
    end
    v=VykhodIzCylBokPov(dn,xl,2*rk,hk,nus,phipn,phipz,phitp); %1 - ����� ����� ���������, 2 - ����� ���. ���-��
    if (v>1)
        po=PopNePop(dlv,2*abs(phipn-nus),razotv); %��������� � �������� ����� ����� - ����� ����� ���. ���-��
    else
        nuss=PoiskVykhUglaSSKonech(dn,abs(nus));  %������� ����� ���������
        po=PopNePop(dlv,abs(phipn-nuss),razotv); %��������� � �������� ����� �����
    end
end
end
vy=po;
end
%������ � ������� ����� ��������� ��� ������� ������
function pr10 = PraCha10(phi,phis)
pr10=sin(phi)/sin(phis)-1; %dn/n1
end
%��������� ������ �� ���. ������ � ������� �� �����-�
function pr20 = PraCha20(phis,phiss)
pr20=cos(phiss)/cos(phis)-1;
end
%���������� ����������� ������� �����
function opn = OprPokPreKBr(la)
lam=[1 2 10 20 25]; pp=[1.487 1.48 1.52 1.48 1.453]; ko=polyfit(lam,pp,length(pp)-1); opn=polyval(ko,la);
end
%����� ���� �����������, phi>0
function pk = PoiskTetaShtrNach(phi,dnr)
a=0; b=pi/2; ep=1e-6; Nit=1e2; k=0; ra=abs(a-b);
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; fa=PraCha10(phi,a)-dnr; fb=PraCha10(phi,b)-dnr; fc=PraCha10(phi,c)-dnr;
    ro=vybPoKo(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
pk=c;
end
%����� ���� ������: v=1 - ����� ����� ���������, v=2 - ����� ���. ���-��
function pk = PoiskVykhUglaSSKonech(dnr,phis)
a=0; b=pi/2; ep=1e-6; Nit=1e2; k=0; ra=abs(a-b); 
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; 
fa=PraCha20(phis,a)-dnr; fb=PraCha20(phis,b)-dnr; fc=PraCha20(phis,c)-dnr;
ro=vybPoKo(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3); k=k+1; ra=abs(a-b);
end
pk=c;
end
%------------------���� �� �������--------------
%------------------���� �� �������--------------
%------------------���� �� �������--------------
function pr = prover(rx,ry,phi)
if ((phi==0) || (phi==pi))
    ry=0;
end
if ((phi==pi/2) || (phi==3*pi/2))
    rx=0;
end
pr=sqrt(rx^2+ry^2);
end
function pd = PopNePop(lam,deltatetass,raot)
p=0; deko=lam/raot;
if (deltatetass>deko)
    p=0;
else p=1;
end
pd=p;
end
function v = VykhodIzCylBokPov(dn,z,dk,hk,hi,phip,phipz,phitp)
vyh=1;
if (((phitp==pi/2) || (phitp==-pi/2)) && (phip==0) && (phipz==0))
    vyh=1;
end
if (((phitp==pi/2) || (phitp==3*pi/2)) && (phip==0) && (phipz==0))
    vyh=1;
end
if (phip==0)
beta=atan(abs((hk-z)/dk));
if (abs(phipz)>beta)
    vyh=1;
else vyh=2;
end
if ((dn<0) && (abs(phip)>beta))
        vyh=1; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
        else vyh=2;
end
end
if (hi==pi/2)
    vyh=2;
end
dlpuos=dk*abs(cos(hi));
dlpuob=abs(hk-z); 
if (dlpuob<dlpuos)
    vyh=1;
else vyh=2;
end
v=vyh;
end
function v = VykhodIzCylOsn(dn,phitp,phin,teta,rltp,rk,hk,hi)
vyh=1;
if (phin==0)
    if (phitp==pi)
beta=atan(abs((rltp-rk)/hk));
    end
    if (phitp==0)
beta=atan(abs((rltp+rk)/hk));
    end
if (abs(hi)>beta)
    vyh=2;
else vyh=1;
end
    if ((dn<0) && (abs(teta)>beta))
        vyh=2; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
        else vyh=1;
    end
else
    dlo=PoiskDlPutiLuchaOsn(phin,rk,rltp*cos(phitp),rltp*sin(phitp));
    if (dlo<hk)
        vyh=2;
    else vyh=1;
    end
end
v=vyh;
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
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
%vyfr - ����� �������: 1 - 60, 2 - 60-100, 3 - 100-150, 4 - 150-200; vyko - ����� �������� ����: 1 - 0,1 %, 2 - 0,2 %, 3 - 0,3 %
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
ro=0; a=1e-5; b=1e5; ep=1e-6; Nit=1e2; k=0; ffv=1e-4; ra=abs(a-b); %ffv = 0,01 % - ����� ���� ������ �� �������� � �������� �������
while ((ra>ep) && (k<Nit))
    c=(a+b)/2; fa=erfc(a)-ffv; fb=erfc(b)-ffv; fc=erfc(c)-ffv; 
    ro=vybPoKo(fa,fb,fc,a,b,c); a=ro(1); b=ro(2); c=ro(3); k=k+1;
    ra=abs(a-b);
end
pk=c;
end
%������������ ������ - ��������� ������� �� ������� ������
function pr1 = PraCha1(p,dp)
prc=dp^2+dp-(1/6)*(dp^3)-(1/24)*(dp^4)+(1/120)*(dp^5)+(1/720)*(dp^6);
prz=1-cot(p)*dp-(1/2)*(dp^2)+(1/6)*cot(p)*(dp^3)+(1/24)*(dp^4)-(1/120)*cot(p)*(dp^5);
pr1=prc/prz; %dn/n1
end
function pr2 = PraCha2(p,dp,dp3)
pr=cot(dp-dp3)-(1/2)*(dp3^2)-(cot(p)^2)*dp*dp3+(1/2)*(dp^2)+(cot(p)*dp)^2; %1 , 2
pr=pr+(1/6)*cot(p)*dp3^3-(1/2)*cot(p)*dp*(dp3^2)-(1/2)*cot(p)*dp3*(dp^2)-(cot(p)^3)/2*dp3*(dp^2)+(5/6)*cot(p)*(dp^3)+(cot(p)^3)*(dp^3); %3
pr=pr+(dp3^4)/24+(1/6)*(cot(p)^2)*dp*(dp3^3)-(1/4)*(dp3^2)*(dp^2)-(1/2)*(cot(p)^2)*(dp3^2)*(dp^2)-(5/6)*(cot(p)^2)*dp3*(dp^3)-(cot(p)^4)*dp3*(dp^3)+(7/6)*(cot(p)^2)*(dp^4)+(5/24)*(dp^4)+(cot(p)^4)*(dp^4); %4
pr=pr-(cot(p)/120)*(dp3^5)+(1/24)*cot(p)*dp*(dp3^4)+(1/12)*cot(p)*(dp^2)*(dp3^3)+(1/6)*(cot(p)^3)*(dp3^3)*(dp^2)-(5/12)*cot(p)*(dp3^2)*(dp^3)-(1/2)*(cot(p)^3)*(dp3^2)*(dp^3)-(5/24)*cot(p)*dp3*(dp^4);
pr=pr-(1/3)*(cot(p)^3)*dp3*(dp^4)+(3/2)*(cot(p)^3)*dp3*(dp^4)+(cot(p)^5)*dp3*(dp^4)+(61/120)*cot(p)*(dp^5)+(3/2)*(cot(p)^3)*(dp^5)+(cot(p)^5)*(dp^5); %5
pr2=pr; %dn/n1
end
function pr4 = PraCha4(dn,n,p,dp)
prn=dn*sin(p)-n*cos(p)*dp-((1/2)*n*sin(p)*(dp^2)+dn*dp*cos(p))+((1/6)*n*cos(p)*(dp^3)-(1/2)*dn*sin(p)*(dp^2));
prn=prn+((1/24)*n*sin(p)*(dp^4)+(1/6)*dn*cos(p)*(dp^3))+((1/24)*dn*sin(p)*(dp^4)-(1/120)*n*cos(p)*(dp^5));
pr4=prn;
end
function pr3 = PraCha3(p,dp,n,dn,de)
prn=dn*cos(p)+n*sin(p)*dp-((1/2)*n*cos(p)*(dp^2)-dn*dp*sin(p))-((1/2)*dn*cos(p)*(dp^2)+(1/6)*n*sin(p)*(dp^3));
prn=prn+((1/24)*n*cos(p)*(dp^4)-(1/6)*dn*sin(p)*(dp^3))+((1/24)*dn*cos(p)*(dp^4)+(1/120)*n*sin(p)*(dp^5));
prk=-n*sin(p)*de-(1/2)*n*cos(p)*(de^2)+(1/6)*n*sin(p)*(de^3)+(1/24)*n*cos(p)*(de^4)-(1/120)*n*sin(p)*(de^5);
pr3=prk-prn;
end
%di - ��� ������ ���, ide - ������� (1) ��� �������. ���-� (2)
function [ ras60 ] = RasFra60(di,ide)
mkbr=[250.239 249.740 249.223 250.395 250.336 249.55]; 
mv=[0.283 0.464 0.812 0.22 0.547 0.777]; 
tol=[0.73 0.72 0.72 0.72 0.71 0.7];
rokbr=2.75; rov=0.49; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %� ��3
vv(k)=mv(k)/(1e3*rov); %� ��3
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %� ���
vkbr(k)=vkbr(k)*((1e4)^3); %� ���3
vv(k)=vv(k)*((1e4)^3); %� ���3
end
mini=0; maxi=6e1; mo=(maxi-mini)/2; %��
si=(mini-maxi)/di; %���
disp('������� ����� 60 ���');
depp=PoisDn(1,1,xv,vv,mo,si,mini,maxi,ide);
ras60=depp; %ras60=alphaRas(tol,length(tol),Tr);
end
function [ ras60_100 ] = RasFra60_100(di,ide)
mkbr=[250 249.629 249.294 249.706 249.510 249.307 250.328 249.604 249.206]; 
mv=[0.255 0.539 0.809 0.295 0.517 0.756 0.36 0.534 0.843]; 
tol=[0.72 0.71 0.7 0.7 0.73 0.72 0.74 0.7 0.76]; 
rokbr=2.75; rov=0.52; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %� ���3
vv(k)=vv(k)*((1e4)^3); %� ���3
end
mini=6e1; maxi=1e2; mo=(maxi-mini)/2; %��
si=(mini-maxi)/di; %���
disp('������� 60-100 ���');
depp=PoisDn(2,1,xv,vv,mo,si,mini,maxi,ide);
ras60_100=depp; %ras60_100=alphaRas(tol,length(tol),Tr);
end
function [ ras100_150 ] = RasFra100_150(di,ide)
mkbr=[249.913 249.607 249.218 249.929 249.695 249.306 250.405 249.625 249.348];
mv=[0.315 0.473 0.709 0.293 0.528 0.83 0.27 0.493 0.764]; 
tol=[0.74 0.74 0.72 0.72 0.71 0.7 0.78 0.73 0.76]; 
rokbr=2.75; rov=0.53; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vkbr(k)=vkbr(k)*((1e4)^3); %� ���3
vv(k)=vv(k)*((1e4)^3); %� ���3
end
mini=1e2; maxi=15e1; mo=(maxi-mini)/2; %��
si=(mini-maxi)/di; %���
disp('������� 100-150 ���');
depp=PoisDn(3,1,xv,vv,mo,si,mini,maxi,ide);
ras100_150=depp; %ras100_150=alphaRas(tol,length(tol),Tr);
end
function [ ras150_200 ] = RasFra150_200(di,ide)
mkbr=[250.882 249.590 249.213 250.299 249.441 249.365];
mv=[0.320 0.533 0.849 0.223 0.502 0.797]; 
tol=[0.76 0.72 0.69 0.73 0.73 0.73]; 
rokbr=2.75; rov=0.56; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); vv(k)=mv(k)/(1e3*rov); xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
vv(k)=vv(k)*((1e4)^3); vkbr(k)=vkbr(k)*((1e4)^3); %� ���3
end
mini=15e1; maxi=2e2; mo=(maxi-mini)/2; %��
si=(mini-maxi)/di; %���
disp('������� 150-200 ���');
depp=PoisDn(4,1,xv,vv,mo,si,mini,maxi,ide);
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
function [ al ] = alphaRas(tol,n,Tr)
alp=0;
for k=1:n
    alp(k)=-log(Tr(k))/tol(k);
end
al=alp;
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
%phin - ���� ����������� �������, xp, yp - ���������� ������� �� �����, ���� � ��������� ���������
function dlot = PoiskDlPutiLuchaOsn(phin,a,xp,yp)
dlo=0; ep=1e-6; Nit=1e2; ko=tan(phin); sc=yp-ko*xp; zn=1;
for j=1:2
zn=-zn; k=0; xa=min([0 zn*a]); xb=max([0 zn*a]); ra=abs(xa-xb); 
while ((ra>ep) && (k<Nit))
    xc=(xa+xb)/2; 
    fa=zn*sqrt(1-xa^2)-(ko*xa+sc); 
    fb=zn*sqrt(1-xb^2)-(ko*xb+sc); 
    fc=zn*sqrt(1-xc^2)-(ko*xc+sc); 
    ro=vybPoKo(fa,fb,fc,xa,xb,xc); 
    xa=ro(1); xb=ro(2); xc=ro(3); k=k+1;
    ra=abs(xa-xb);
end
xpe(j)=xc; ype(j)=ko*xc+sc;
end
v1=[xpe(1)-xp ype(1)-yp]; v2=[xpe(2)-xp ype(2)-yp];
v3=[a*cos(phin) a*sin(phin)]; 
if (((v1*v3')>0) && ((v2*v3')<0))
dlo=sqrt(v1(1)^2+v1(2)^2);
end
if (((v1*v3')<0) && ((v2*v3')>0))
dlo=sqrt(v2(1)^2+v2(2)^2);
end
dlot=dlo;
end
function phitp=poisphitp(eta,ksi)
    if (((eta<0) && (ksi<0)) || ((eta>0) && (ksi<0)))
    phitp=pi+atan(eta/ksi); 
    end
    if (((eta>0) && (ksi>0)) || ((eta<0) && (ksi>0)))
    phitp=atan(eta/ksi); 
    end
end
function phitp=poisphitpn(rlx)
    if (rlx>0) 
        phitp=acos(rlx)*(-1+2*rand()); 
    else
        phitp=pi+acos(abs(rlx))*(-1+2*rand()); 
    end
end
function vy = VykhodIzOsnPP(nlxy,tpl,ak,bk,hk)
vyh=1;
r0=[ak bk]; Rmr0=[tpl(1)-r0(1) tpl(2)-r0(2)]; nva=[-ak 0]; nvb=[0 -bk]; %������������ �������
roav=abs(Rmr0(1)*nva(2)-Rmr0(2)*nva(1))/ak; robp=abs(Rmr0(1)*nvb(2)-Rmr0(2)*nvb(1))/bk;
r0=[ak 0]; Rmr0=[tpl(1)-r0(1) tpl(2)-r0(2)]; roan=abs(Rmr0(1)*nva(2)-Rmr0(2)*nva(1))/ak; 
r0=[0 bk]; Rmr0=[tpl(1)-r0(1) tpl(2)-r0(2)]; robl=abs(Rmr0(1)*nvb(2)-Rmr0(2)*nvb(1))/bk; %���������� �� ������ ��������������
untpb=atan((-bk+tpl(2))/tpl(1))+pi; untpab=atan((bk-tpl(2))/(ak-tpl(1))); unl=atan(nlxy(2)/nlxy(1));
if (nlxy(2)>0)
if (unl<0)
    unl=unl+pi;
end
if ((unl>untpab) && (unl<untpb))
pulu=roav/sin(unl);
end
if ((unl>=untpb) && (unl<=pi))
    unl=pi-unl;
    pulu=abs(robl/cos(unl));
end
if ((unl>=0) && (unl<=untpab))
pulu=abs(robp/cos(unl));
end
end
if (nlxy(2)<0)
    untpb=atan((tpl(2))/tpl(1))-pi;
    untpa=atan(-tpl(2)/(ak-tpl(1)));
if (unl>0)
    unl=unl-pi;
end
if ((unl<untpa) && (unl>untpb))
pulu=abs(roan/sin(unl));
end
if ((unl<=untpb) && (unl>=-pi))
    unl=pi+unl;
    pulu=abs(robl/cos(unl));
end
if ((unl<=0) && (unl>=untpa))
pulu=abs(robp/cos(unl));
end
end
if ((unl==-pi/2) || (unl==pi/2))
    if (nlxy(2)<0)
    pulu=abs(tpl(1));
    else pulu=abs(tpl(1)-ak);
    end
end
if ((unl==pi) || (unl==0))
    if (nlxy(1)<0)
    pulu=abs(tpl(2));
    else pulu=abs(tpl(2)-bk);
    end
end
if (pulu<hk)
    vyh=1;
else vyh=2; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
end
vy=vyh;
end
function [ povo ] = Povoroty(phi,teta,phiz,kotp)
xtp=kotp(1); ytp=kotp(2); ztp=kotp(3);
    xtp1=xtp; ytp1=ytp*cos(phi)+ztp*sin(phi); ztp1=-ytp*sin(phi)+ztp*cos(phi); %������ ������� � ��������� YZ ������ ��� X
    xtp2=xtp1*sin(teta)+ztp1*cos(teta); ytp2=ytp1; ztp2=ztp1*sin(teta)-xtp1*cos(teta); %������ ������� � ��������� XZ ������ ��� Y
    xtp3=xtp2*cos(phiz)-ytp2*sin(phiz); ytp3=xtp2*sin(phiz)+ytp2*cos(phiz); ztp3=ztp2; %������ ������� � ��������� XY ������ ��� Z
    povo=[xtp3  ytp3    ztp3];
end
function [ upapr ] = opreUgPadPre(phi,teta,phiz,veni,dn)
    ugpalu=[-sin(phi)*cos(teta)*cos(phiz)-cos(phi)*sin(phiz)    -sin(phi)*cos(teta)*sin(phiz)+cos(phi)*cos(phiz)    -sin(phi)*sin(teta)]; %���� (�����������) ������� ����
    if (veni==1) %1 - ������ ��������� - b, (-1,0,0); 2 - ������� - a, (0,-1,0)
        phiup=acos(abs(-ugpalu(1))); %�������� � ��������� YZ
        nalu=[ugpalu(2) ugpalu(3)]; pxyz=1; %����������� ���� �� ������� ����������� YZ
    else phiup=acos(abs(-ugpalu(2))); %�������� � ��������� XZ
        nalu=[ugpalu(1) ugpalu(3)]; pxyz=2; %����������� ���� �� ������� ����������� XZ
    end
    phiups=PoiskTetaShtrNach(phiup,dn);  %phiups - ���� �����������, phiup - ���� �������
    upapr=[phiup    phiups  nalu(1) nalu(2) pxyz]; %����� ������� �� ��������� XY � ����������� ����
end
function vy = VykhodIzOsnBokPovPP(ak,bk,hk,narlu,tpl,veni) %tpl - ���������� ����� ������� � �������� ����������� XYZ
vyh=1; %1 - ����� ����� ���������, 2 - ����� ���. ���-��
if (veni==1) %� ������ ��������� - YZ 
    vyh=VybVykhOsnBokPov(bk,hk,ak,tpl,narlu);
else %� ������� ��������� - XZ 
    vyh=VybVykhOsnBokPov(ak,hk,bk,tpl,narlu);
end
    vy=vyh;
end
function vy = VybVykhOsnBokPov(ak,hk,bk,tpl,narlu)
r0=[ak -hk]; Rmr0=[tpl(1)-r0(1) tpl(2)-r0(2)]; nva=[-ak 0]; nvh=[0  hk]; %������������ �������
roan=abs((Rmr0(1)*nva(2)-Rmr0(2)*nva(1))/ak); rohp=abs((Rmr0(1)*nvh(2)-Rmr0(2)*nvh(1))/hk); %���������� �� ������ ��������������: ������ � ������
r0=[ak 0]; Rmr0=[tpl(1)-r0(1) tpl(2)-r0(2)]; roav=abs((Rmr0(1)*nva(2)-Rmr0(2)*nva(1))/ak); %�������
r0=[0 -hk]; Rmr0=[tpl(1)-r0(1) tpl(2)-r0(2)]; rohl=abs((Rmr0(1)*nvh(2)-Rmr0(2)*nvh(1))/hk); %�����
untpa=(ak-tpl(1))/(-tpl(2)); untpa=atan(untpa); unl=narlu(2)/narlu(1); unl=atan(unl);
if (narlu(2)>0)
    untp=tpl(2)/tpl(1); untp=atan(untp); 
if (untp<0)
    untp=untp+pi;
end
if (unl<0)
    unl=unl+pi;
end
if ((unl>untpa) && (unl<untp))
pulu=roav/sin(unl); osbo=1; %osbp: 1 - ������ ���������, 2 - ����� ���������, 3 - ������� ����������� �����, 4 - ������, 5 - ��������
end
if ((unl>=untp) && (unl<=pi))
    unl=pi-unl; pulu=abs(rohl/cos(unl)); osbo=3;
end
if ((unl>=0) && (unl<=untpa))
pulu=abs(rohp/cos(unl)); osbo=4;
end
end
if (narlu(2)<0)
    untph=atan((-hk-tpl(2))/(-tpl(1)));
    if (untph>0)
    untph=untph-pi;
    end
    untpah=atan((-hk-tpl(2))/(ak-tpl(1)));
if (unl>0)
    unl=unl-pi;
end
if ((unl<untpah) && (unl>untph))
pulu=abs(roan/sin(unl)); osbo=2;
end
if ((unl<=untph) && (unl>=-pi))
    pulu=abs(rohl/cos(unl+pi)); osbo=3;
end
if ((unl<=0) && (unl>=untpah))
pulu=abs(rohp/cos(unl)); osbo=4;
end
end
if ((unl==-pi/2) || (unl==pi/2))
    if (narlu(2)>0)
    pulu=abs(tpl(2));
    else pulu=abs(-hk-tpl(2));
    end
end
if ((unl==pi) || (unl==0))
    if (narlu(1)<0)
    pulu=abs(tpl(1));
    else pulu=abs(ak-tpl(1));
    end
end
if (pulu>bk) %pulu - ���� ���� �� ���������
    vyh=2;
else
    vyh=1;
    switch (osbo)
        case {3,4}
    vyh=3;         
        case {1,2}
    vyh=1;
        otherwise
    vyh=2;
    end %1 - ����� ����� ���������, 2 - ����� ���. ���-�� �������� 3 - ����� ������� ���-�� �����
end    
vy=vyh;
 end