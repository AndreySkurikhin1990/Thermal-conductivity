function [ mti ] = ModelirTeploperItom(ttss,ttsv,tem)
nom=1; %выбор Г и Х темп., ПТП для моделирования т/о в ИТОМ
no=1; te0=273.15; %выбор номера ИТОМ
h=30e-3; Thna=32e1-te0; Thde=1e-1;
temi=0; temi=1e2:1e2:5e2; n=length(temi); %моделирование процесса теплообмена в ИТОМ
qob=[2723,4235,5796]; te=[633,773,907]-te0;
te200=[0.09,0.12,0.183,0.23];
te380=[0.12,0.139,0.194,0.25];
tei=[200,380]+te0;
kitom440=polyfit(tei,[te200(1),te380(1)],1);
kitom620=polyfit(tei,[te200(2),te380(2)],1);
kitom860=polyfit(tei,[te200(3),te380(3)],1);
kitom1000=polyfit(tei,[te200(4),te380(4)],1);
switch no
    case 1
        kti=kitom440;
    case 2
        kti=kitom620;
    case 3
        kti=kitom860;
    case 4
        kti=kitom1000;
end
kti=kti';
ep=1e-6; nit=1e6; d=1e-2; koeq=polyfit(te,qob,2);
for k=1:n
del=1; j=1; ts=temi(k); laef(k)=polyval(kti,ts); qo(k)=polyval(koeq,ts); Thna=Thna+(k-1)*Thde;
while ((del>ep) && (j<nit))
Th(k)=Thna+(j-1)*d; Tg(k)=Th(k)+abs(qo(k)*h/laef(k)); del=(2*ts-Tg(k)-Th(k)); j=j+1;
end
tes(k)=(Th(k)+Tg(k))/2;
end
thitom=Th(nom);
ptpitom=qo(nom);
tgitom=Tg(nom);
laef=laef'; tesr = tes(nom);
posh=21.8e-2; pov = 66.35e-2; %пористость шамота и вермикулита
hpv = srRaPorVerm(); htv = (1 - pov) * hpv / pov;
hpsh = SreRazPorSha(); htsh = (1 - posh) * hpsh / posh;
sv=[60,50,25,20]*1e-2; %содержание вермикулита
psv = sv(no); pssh = 1 - psv; csi = 1e2; csv = floor(psv * csi); cssh = floor(pssh * csv); 
hob = (hpv + htv) * csv + (hpsh + htsh) * cssh; 
cpsi = floor(h / hob) - 2; 
nsshl = floor(cssh / 2); nsshr = csi - nsshl; tk = tgitom; titoms = 0; %titomss=0; 
x=0; x(1) = 0; dono=0; dh=1e-4; koo = 0:dh:(h-dh); dtdx=abs(thitom-tgitom)/h; ts = tk; %распределение температуры линейно
for k=1:length(koo)
    rastemitom(k)=tgitom-dtdx*koo(k);
end
for k=1:cpsi
for p=1:csi
    l = p; p = dono + p;
    if ((l<=nsshl) || (l>=nsshr))
        ttsshk=opredKTPTKTochSha(ttss, tem, tk, length(tem));
        tk = tk - ptpitom * htsh / ttsshk;
        titoms(p) = tk;  %расчет температуры;
        ktpvo = opredKTPVozd(tk);
        tk = tk - ptpitom * hpsh / ktpvo;
        %titomss(p) = tk;
        hob = (hpsh + htsh);
    else
        ttsvk=opredKTPTKTochSha(ttsv, tem, tk, length(tem));
        tk = tk - ptpitom * htv / ttsvk;
        titoms(p) = tk;
        ktpvo = opredKTPVozd(tk);
        tk = tk - ptpitom * hpv / ktpvo;
        %titomss(p) = tk;
        hob = (htv + hpv);
    end
    if (p > 1)
            x(p) = x(p-1) + hob;
    end
    p = l;
end
dono=dono + p;
end
ts =  (ts + tk) / 2;
kxt = polyfit(x,titoms,1);
ektp = ptpitom*h/(polyval(kxt,0) - polyval(kxt,h));
ektp = polyval(kti,tesr);
pl=plot(x,titoms,'-b',koo,rastemitom,'-k');
set(pl,'LineWidth',3); hold on; grid on; 
xlabel({'Координата, м'}); 
ylabel({'Температура, К'}); 
title({'График зависимости температуры от координаты'});
legend('Результат моделирования','Результат аппроксимации','Location','Best');
mti=[0];
end

function ktpo = opredKTPTKTochSha(ktptks, te, temp, n)
f = 1; p = 0; ep=1e-4; ktp = 0;
if ((temp>te(1)) && (temp<te(n)))
for k = 1:n
if ((te(k) >= temp) && (f>0))
        p = k; f = 0; break;
end
end
elseif (temp>=te(n))
    p=n; f=0;
elseif (temp<=te(1))
    p=2; f=0;
end
if ((f==0) && (p>1))
    x2=te(p);
    x1=te(p-1);
	dt = x2 - x1;
if (abs(dt) > ep)
    y2=ktptks(p);
    y1=ktptks(p - 1);
    b=y1;
			ko = (y2 - y1) / dt;
            if (p==n)
                b=y2;
            end
			ktp = b + ko*(temp - x1);
else 
    ktp=0;
end
end
ktpo=ktp;
end

function opr = opredKTPVozd(T)
te0=273.15;
te = arrTempAir()+te0;
lamb = koefTeploprovAir();
n = length(te);
lam=opredKTPTKTochSha(lamb,te,T,n);
opr = lam;
end