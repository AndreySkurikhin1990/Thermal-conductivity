function [ oktptks ] = opredKTPTverKarkitom(vmitom,tem)
te0=273.15; srk=RazPorItom(vmitom); srp=srk(1); por=srk(3);
y0=30e-3; %qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); 
lentem=length(tem);
koef=poisMasKoefItom(vmitom);
for k=1:lentem
ts=tem(k); laef=koef(2)*ts+koef(1); laefm(k)=laef;
end
x0=srp/(por^(1/3))/2;
cvp=1+17; cve=cvp+12; cved=cve+1; na = 1;
srk=0;
for k=1:lentem
    srk(k)=0;
end
dkosc=oprmassialmgitom(vmitom,0);
dkosct=oprmassialmgitom(vmitom,1);
stchk=oprmassialmgitom(vmitom,2);
wsio=stchk(1); walo=stchk(2); wmgo=stchk(3);
kusci=oprStCheritom(vmitom,tem,3,0);
tkusci=oprStCheritom(vmitom,tem,4,0);
stchk=0; stchk=oprStCheritom(vmitom,tem,1,1)
cvp=1+17+4; cve=cvp+12+2; cved=cve+4; f=(cve-cvp)/2; na=1;
s=0; p=0; ep=1e-4;
for k=1:lentem
    lavo(k)=opredTeploprovVozd(tem(k));
    for j=1:cved
        sr(j,k)=0;
    end
end
tena=3e2; dtos=1e2; 
srk=DulnevZernkvi(por, tem, srk, na, laefm, lavo, wsio, walo, wmgo, x0, lentem, tena, dtos, length(kusci), kusci, tkusci, stchk,x0);
if (sum(srk)>ep)
srk=srk'
for k=1:length(srk)
sr(na,k)=srk(k); 
end
srk=0;
end
na=na+1; 
for k=na:cved
if (k<cvp) 
        srk=opredTvChaKoeTepSrInSpSha(k, por, tem, te0, laefm, x0)
elseif (k<cve) 
        srra=oprEdinVelkvi(5, vmitom); 
        prgr=oprEdinVelkvi(7, vmitom); 
        legr=oprEdinVelkvi(6, vmitom); 
        raspr=oprEdinVelkvi(3, vmitom);
        tm=oprEdinVelkvi(1, vmitom); srp=tm(1); %dmi=length(srra);
        srk=DulnevKoefTepVermN(k-cvp, por, tem, te0, laefm, srp, 0, f, 1, raspr, legr, prgr, srra)
elseif (k<=cved) 
        srk=VasFrayObsh(laefm, lavo, por, lentem);
        if (sum(srk)>ep)
            srk=srk';
        for j=1:length(srk)
            sr(k,j)=srk(j); 
        end
        end
end
srk=0;
end
srzn=0;
if (sum(sr)>ep)
q=1; p=0; 
for j=1:cved
	f=1; 
    for k=1:lentem
if (sr(j,k)<ep)
    f=0; 
    break; 
end
    end
if (f>0)
    no(q)=j; q=q+1; p=p+1;
end
end
for j=1:lentem
    s=0; 
    for k=1:p
        s=s+sr(no(k),j); 
    end
        if (abs(p)>0) 
            srzn(j)=s/p; 
        else
            srzn(j)=0;
        end
end
oktptks=srzn;
else oktptks=0;
end
end

function [ ktpi ] = poisMasKoefItom(no)
te200 = 2e2 + tei0; 
te380 = 38e1 + tei0;
switch (no)
    case (0)
ktpit(1) = 9e0*1e-2; ktpit(2) = 12e0*1e-2;
kitom440(2) = (ktpit(2) - ktpit(1)) / (te380 - te200); 
kitom440(1) = ktpit(1) - kitom440(2) * te200;
kti = kitom440;
    case (1)
		%ktpit(1) = 12.0*1e-2; ktpit(2) = 139.0*1e-3; %из Диссертации
		ktpit(1) = 18.0*1e-2; ktpit(2) = 19.0*1e-2; %Данные 2017 года
		kitom620(2) = (ktpit(2) - ktpit(1)) / (te380 - te200); 
        kitom620(1) = ktpit(2) - kitom620(2) * te380;
		kti = kitom620;
    case (2)
		%ktpit(1) = 183.0*1e-3; ktpit(2) = 194.0*1e-3; %из Диссертации
		ktpit(1) = 26.0*1e-2; ktpit(2) = 37.0*1e-2; %Данные 2017 года
		kitom860(2) = (ktpit(2) - ktpit(1)) / (te380 - te200); 
        kitom860(1) = ktpit(2) - kitom860(2) * te380;
		kti = kitom860;
    case (3)
		%ktpit(1) = 23.0*1e-2; ktpit(2) = 25.0*1e-2; %из Диссертации
		ktpit(1) = 42.0*1e-2; ktpit(2) = 52.0*1e-2; %Данные 2017 года
		kitom1000(2) = (ktpit(2) - ktpit(1)) / (te380 - te200); 
        kitom1000(1) = ktpit(2) - kitom1000(2) * te380;
		kti = kitom1000;
    otherwise
        disp('Net takoy marki ITOM!');
end
	ktpi=kti;
end

function [ om ] = oprmassialmgitom(vybitom, vvm)
ko=1e-2; tnd = 6e2; dtd = 2e2; dkoscil = 6; cemi = 11;
pori440 = (8e1+82e0)*1e-2/2e0;
pori620 = (75e0+78e0)*1e-2/2e0;
pori860 = (65e0+68e0)*1e-2/2;
pori1000 = (62e0+65e0)*1e-2/2e0;
dkoscit(1) = tnd; 
for k = 2:dkoscil
    dkoscit(k) = dkoscit(k - 1) + dtd;
end
switch (vybitom)
    case (0)
		saloi = 26e0*ko; smgoi = 22e0*ko; ssioi = 52e0*ko; poritom = pori440; %для трехкомпонентной смеси
		wal = 23e0*ko; wmg = 19e0*ko; wsi = 49e0*ko; %для многокомпонентной смеси
		dkoscim(1) = 6.07; dkoscim(2) = 5.36; dkoscim(3) = 6.19;
		dkoscim(4) = 13.48; dkoscim(5) = 19.93; dkoscim(6) = 27.69; %ИТОМ-440
    case (1)
		saloi = 29e0*ko; smgoi = 16e0*ko; ssioi = 55e0*ko;
		wal = 26e0*ko; wmg = 15e0*ko; wsi = 5e1*ko; poritom = pori620;
		dkoscim(1) = 6.41; dkoscim(2) = 5.53; dkoscim(3) = 6.32;
		dkoscim(4) = 13.56; dkoscim(5) = 19.86; dkoscim(6) = 27.66; %ИТОМ-620
    case (2)
		saloi = 33e0*ko; smgoi = 11e0*ko; ssioi = 56e0*ko;
		wal = 3e1*ko; wmg = 1e1*ko; wsi = 52e0*ko; poritom = pori860;
		dkoscim(1) = 7.28; dkoscim(1) = 6.1; dkoscim(2) = 6.83;
		dkoscim(3) = 14.02; dkoscim(4) = 19.86; dkoscim(5) = 27.22; %ИТОМ-860
    case (3)
		saloi = 35e0*ko; smgoi = 9e0*ko; ssioi = 56e0*ko; 
		wal = 33e0*ko; wmg = 7e0*ko; wsi = 53e0*ko; poritom = pori1000;
		dkoscim(1) = 7.45; dkoscim(2) = 6.24; dkoscim(3) = 6.95;
		dkoscim(4) = 14.15; dkoscim(5) = 19.89; dkoscim(6) = 27.09; %ИТОМ-1000
    otherwise
        disp('Net takoy marki ITOM!');
end
        for k = 1:dkoscil
		tm = dkoscim(k)*ko; 
        dkoscim(k) = 1e0 - tm;
        end
        switch (vvm)
            case (0)
                om=dkoscim;
            case (1)
                om=dkoscit;
            case (2)
                om=[wal,wmg,wsi,saloi,smgoi,ssioi,poritom];
        end
end

function [ rt ] = RazPorItom(no) %расчет размера пор %выбор номера ИТОМ
po=[(80+82)/2,(75+78)/2,(65+68)/2,(62+65)/2]*1e-2; %пористость
rapo=[150,95,85,75,65,55,45,35,25,15,7.5,4,2,0.75,0.3,0.055,0.0055]*1e-6; %в метрах
por=po(no+1);
switch no
    case 0
    rpr = [0.57,0.3,0.6,0.87,1.35,2.07,3.72,3.81,5.38,7.6,9.67,10.87,34.68,10.78,6.25,1.47,0]*1e-2; %распределение частиц по размерам
    case 1
    rpr = [0.26,0.15,0.26,0.22,0.81,0.81,1.88,3.95,5.54,7.35,7.09,9.01,34.9,13.59,7.5,2.14,4.54]*1e-2;
    case 2
    rpr = [0.4,0.09,0.44,0.22,0.66,1.02,1.33,2.66,4.07,10.71,12.17,11.29,35.06,11.24,7.13,1.51,0]*1e-2;
    case 3
    rpr = [0.23,0.19,0.04,0.61,0.23,1.03,0.8,2.47,5.66,10.87,14.18,12.5,32.61,11.59,5.25,1.75,0]*1e-2;
end
srp=0;
for k=1:length(rpr)
    srp = srp + rpr(k) * rapo(k);
end
srts = (1 - por) * srp / por; %средний размер твердого скелета
rt=[srp,srts,po];
end