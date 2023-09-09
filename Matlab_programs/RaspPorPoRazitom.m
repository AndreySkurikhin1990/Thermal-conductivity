function [ t ] = RaspPorPoRazitom(vvi,vyb)
rpn=1e0; dp=1e0; 
srp=0; l=0; h=1e-6; p=h; qg=2e2;
prgr0=PravGranPoromitom(); 
legr0=LevGranPoromitom(prgr0);
rpr0=NapMasRaspitom(vvi);
n=length(rpr0);
raspr=rasrasporporaz0(prgr0, legr0, rpr0, rpn, 0, n, dp);
qg=length(raspr);
for k=1:qg
    srra(k)=(p+l)/2e0; 
    srp=srp+srra(k)*raspr(k); 
    legr(k)=l; 
    prgr(k)=p; 
    p=p+h; 
    l=l+h;
end
srp=(srp+SreRazPoritom(raspr, srra))/2e0;
switch (vyb)
    case (0)
        t=legr0;
    case (1)
        t=prgr0;
    case (2)
        t=raspr;
    case (3)
        t=srra;
    case (4)
        t=srp;
    otherwise
        t=[0];
end
end

%»“ŒÃ-440
function [ t ] = rasPorpoRazmitom440()
rapo=[0.57, 0.3, 0.6, 0.87, 1.35, 2.07, 3.72, 3.81, 5.38 ...
    7.6, 9.67, 10.87, 34.68, 10.78, 6.25, 1.47, 0.0];
s=sum(rapo); 
rapo=rapo*1e2*1e-2/s;
rapo=fliplr(rapo);
t=rapo;
end
%»“ŒÃ-620
function [ t ] = rasPorpoRazmitom620()
rapo=[0.26, 0.15, 0.26, 0.22, 0.81, 0.81, 1.88 ...
    3.95, 5.54, 7.35, 7.09, 9.01, 34.9, 13.59, 7.5, 2.14, 4.54];
s=sum(rapo); 
rapo=rapo*1e2*1e-2/s;
rapo=fliplr(rapo);
t=rapo;
end
%»“ŒÃ-860
function [ t ] = rasPorpoRazmitom860()
rapo=[0.4, 0.09, 0.44, 0.22, 0.66, 1.02, 1.33 ...
    2.66, 4.07, 10.71, 12.17, 11.29, 35.06, 11.24, 7.13, 1.51, 0.0];
s=sum(rapo); 
rapo=rapo*1e2*1e-2/s;
rapo=fliplr(rapo);
t=rapo;
end
%»“ŒÃ-1000
function [ t ] = rasPorpoRazmitom1000()
rapo=[0.23, 0.19, 0.04, 0.61, 0.23, 1.03, 0.8 ...
    2.47, 5.66, 10.87, 14.18, 12.5, 32.61, 11.59, 5.25, 1.75, 0.0];
s=sum(rapo); 
rapo=rapo*1e2*1e-2/s;
rapo=fliplr(rapo);
t=rapo;
end
function [ t ] = PravGranPoromitom()
rapo=[15e1, 95, 85, 75, 65, 55, 45, 35, 25, 15, 7.5, 4, 2, 0.75, 0.3, 0.055, 0.0055]; 
rapo=(1e-6)*rapo; 
rapo=fliplr(rapo);
t=rapo;
end
function [ t ] = LevGranPoromitom(pravgr, levgr)
n=length(pravgr);
levgr(1)=0; 
for k=2:n
    levgr(k)=pravgr(k-1); 
end
t=levgr;
end

function t = SreRazPoritom(rpr, srra)
srpo=0; 
n=length(srra);
for k=1:n
    if (rpr(k)>0) 
        srpo=srpo+rpr(k)*srra(k); 
    end
end
t=srpo/sum(rpr);
end

function [ r ] = NapMasRaspitom(vvi)
if (vvi==0) 
    r=rasPorpoRazmitom440(); %»“ŒÃ-440
elseif (vvi==1) 
    r=rasPorpoRazmitom620(); %»“ŒÃ-620
elseif (vvi==2) 
    r=rasPorpoRazmitom860(); %»“ŒÃ-860
elseif (vvi==3) 
    r=rasPorpoRazmitom1000(); %»“ŒÃ-1000
end
end

function [ t ] = rasrasporporaz0(prgr0, legr0, raspr0, rpn, vyb, n, de)
e=1e-1; ht=1e0; e1=ht+e; koef=1e6;
prgr00=prgr0*koef;
legr00(1)=0; 
for k=2:n
    legr00(k)=prgr00(k-1);
end
s=0; 
for k=1:n
    if (prgr00(k)<e1) 
        p=k; s=s+raspr0(k); %‡ÁÏÂ ÔÓ ‰Ó 1 ÏÍÏ
    end
end
k=1; rprm(k)=s; q=1; 
for k=p:n 
    r=prgr00(k)-legr00(k);
    if ((prgr00(k)>e1) && (r<e1))
		rprm(q)=raspr0(k);
		q=q+1;
    end
    if (r>e1)
        m=raspr0(k)/r;
		t=r; j=0; 
        while (t>e) 
            j=j+1; 
            t=t-ht;
        end
            jk=j;
        for j=1:jk
		rprm(q)=m;
		q=q+1; 
        end
    end
end
prgr01=rpn; n=q-1; qg=n; 
for k=1:qg
    prgrm(k)=prgr01;
if (k>1) 
    legrm(k)=prgrm(k-1);
else
    legrm(k)=0;
end
    ras01(k)=rprm(k);
    prgr01=prgr01+de;
end
%prgrm=prgrm
%ras01=ras01
if (vyb==0) 
    t=ras01;
elseif (vyb==1) 
    t=prgrm; 
elseif (vyb==2) 
    t=legrm;
end
end