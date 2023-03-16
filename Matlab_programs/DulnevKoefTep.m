function [ dul ] = DulnevKoefTep()
format long g;
Te=arrTem_VVF2();
%lae = [0.180183842490065   0.244572803850782   0.127756995443625   0.174128767881838   0.104329016643457   0.214722590630864];
%tee = [620.15  890.65  627.0875    882.4   620.525     905.3375];
te0 = 273.15;
Te = Te + te0;
laef=arrKTP_VVF2();
kola=polyfit(Te,laef,length(Te)-1);
T=1e2+1:1e-1:9e2-1; T=T+te0;
for k=1:length(T)
    lae(k)=polyval(kola,T(k));
end
%pl=plot(T,lae,'-b');
%set(pl,'LineWidth',2); hold on; grid on;
%k=0; po = 0.5575; %8-4 mm
k=1; por = 0.6635; %2-0.7 mm
l = srRaPor(k);
lb = l / (por^(1 / 3));
n=length(T);
lamadi = 0; lamizo = 0; 
lamodo=0; lamcyladi=0; 
lamcylizo = 0; lamcub=0; 
lamcyl=0; lamsr=0;
for k=1:n
lamadi(k) = DulKoefTep1Adi(T(k), por, lae(k));
lamizo(k) = DulKoefTep1Izoterm(T(k), por, lae(k));
lamodo(k) = DulKoefTep1Odoev(T(k), por, lae(k));
lamcyladi(k) = DulKoefTep1AdiCyl(T(k), por, lae(k), lb);
lamcylizo(k) = DulKoefTep1IzotermCyl(T(k), por, lae(k), lb);
lamcub(k)=(lamadi(k)+lamizo(k))/2;
lamcyl(k)=(lamcyladi(k)+lamcylizo(k))/2;
lamsr(k)=(lamcub(k)+lamcyl(k)+lamodo(k))/3;
end
lamadi=lamadi'
lamizo=lamizo'
lamodo=lamodo'
lamcyladi=lamcyladi'
lamcylizo=lamcylizo'
lamcub=lamcub'
lamcyl=lamcyl'
T=T-te0;
%pl=plot(T,lamadi,'-b',T,lamizo,'-k',T,lamodo,'-c',T,lamcyladi,'-g',T,lamcylizo,'-r');
%%pl=plot(T,lamsr,'-k');
pl=plot(T,lamcub,'-b',T,lamodo,'-k',T,lamcyl,'-g');
set(pl,'LineWidth',2); hold on; grid on;
xlabel({'Температура, °C'}); 
ylabel({'КТП твердой фазы, Вт/(м*К)'}); 
%title({'График зависимости КТП ТФ от T - Адиаб. разб-е (цил. поры)'});
title({'График зависимости КТП ТФ от температуры'});
%legend('Адиабатическое разбиение (куб.)','Изотермическое разбиение (куб.)','По Оделевскому', 'Адиабатическое разбиение (цил.)','Изотермическое разбиение (цил.)');
legend('Куб','По Оделевскому', 'Цилиндр','position','best');
dul = lamsr';
end

function SRP = srRaPor(n)
rp = 0;
legr = LevGranPorom();
prgr = PravGranPorom();
sr = (legr + prgr) / 2;
d = 0;
if (n == 0)
rp = rasPorpoRazm84();
elseif (n == 1)
rp = rasPorpoRazm207();
end
le = length(sr);
for k=1:le
    d = d + rp(k) * sr(k);
end
SRP = d;
end

function du1 = DulKoefTep1Adi(T, por, lame)
lama = 1e4;
lamb = 1e-6;
ep = 1e-7;
%koe = (lae(2) - lae(1)) / (tee(2) - tee(1));
%lame = lae(1) + koe * (T - tee(1));
lam2 = opredTeploprovVozd(T);
k=0;
while (abs(lama-lamb) > ep)
    lamc = (lama + lamb) / 2;
    nua = lam2 / lama;
    nub = lam2 / lamb;
    nuc = lam2 / lamc;
    fa = (nua - (nua - 1) * (por^(1 / 3)) * (1 - por^(2 / 3))) / (nua - (por^(1/3)) * (nua - 1)) - lame / lama;
    fb = (nub - (nub - 1) * (por^(1 / 3)) * (1 - por^(2 / 3))) / (nub - (por^(1/3)) * (nub - 1)) - lame / lamb;
    fc = (nuc - (nuc - 1) * (por^(1 / 3)) * (1 - por^(2 / 3))) / (nuc - (por^(1/3)) * (nuc - 1)) - lame / lamc;
    if (fc*fb > 0) 
        if (fa*fc < 0) lamb=lamc; end; end;
    if (fc*fa > 0) 
        if (fb*fc < 0) lama=lamc; end; end;
    k=k+1;
    if (k>1e2) 
        break; end;
end
du1 = lamc;
end

function du2 = DulKoefTep1Izoterm(T, por, lame)
lama = 1e4;
lamb = 1e-6;
ep = 1e-10;
%koe = (lae(2) - lae(1)) / (tee(2) - tee(1));
%lame = lae(1) + koe * (T - tee(1));
lam2 = opredTeploprovVozd(T);
k=0;
while (abs(lama - lamb) > ep)
    lamc = (lama + lamb) / 2;
    nua = lam2 / lama;
    nub = lam2 / lamb;
    nuc = lam2 / lamc;
    fa = (1 + (nua - 1) * (por^(2 / 3)) ) / (1 + (por^(2 / 3)) * (nua - 1)) * (1 - por^(1 / 3)) - lame / lama;
    fb = (1 + (nub - 1) * (por^(2 / 3)) ) / (1 + (por^(2 / 3)) * (nub - 1)) * (1 - por^(1 / 3)) - lame / lamb;
    fc = (1 + (nuc - 1) * (por^(2 / 3)) ) / (1 + (por^(2 / 3)) * (nuc - 1)) * (1 - por^(1 / 3)) - lame / lamc;
    if (fc*fb > 0) 
        if (fa*fc < 0) lamb=lamc; end; end;
    if (fc*fa > 0) 
        if (fb*fc < 0) lama=lamc; end; end;
    k=k+1;
    if (k>1e2) 
        break; end;
end
du2 = lamc;
end

function du31 = DulKoefTep1OdolevMatr(lam2, m2, lame)
lama = 1e4; lamb = 1e-9; ep = 1e-7; k=0; nit=1e2; ra=1e4;
while ((ra > ep) && (k < nit))
    lamc = (lama + lamb) / 2;
    nua = 1 - lama / lam2; %lam2 - КТП воздуха
    nub = 1 - lamb / lam2;
    nuc = 1 - lamc / lam2;
    fa = 1 - (1 - m2) / (1 / nua - m2 / 3) - lame / lam2;
    fb = 1 - (1 - m2) / (1 / nub - m2 / 3) - lame / lam2;
    fc = 1 - (1 - m2) / (1 / nuc - m2 / 3)  - lame / lam2;
    if ((fc*fb > 0) && (fa*fc < 0)) 
        lamb=lamc; 
    end
    if ((fc*fa > 0) && (fb*fc < 0)) 
        lama=lamc; 
    end
    k=k+1;
    ra=abs(lama - lamb);
end
du31 = lamc;
end
%нет регулярной структуры, статистическая смесь, частицы распределены хаотически
function du32 = DulKoefTep1OdolevStatSm(lam2, m2, lame)
lama = 1e4; lamb = 1e-9; ep = 1e-7; k=0; nit=1e2; ra=1e4; %lam2 - КТП воздуха, ищем КТП твердого скелета
while ((ra > ep) && (k < nit))
    lamc = (lama + lamb) / 2;
    v1 = m2; v2 = 1 - m2;
    nua = ((3*v1 - 1) * lama + (3*v2 - 1) * lam2) / 4;
    nub = ((3*v1 - 1) * lamb+ (3*v2 - 1) * lam2) / 4;
    nuc = ((3*v1 - 1) * lamc + (3*v2 - 1) * lam2) / 4;
    fa = nua + sqrt (nua^2 + lama * lam2 / 2) - lame;
    fb = nub + sqrt (nub^2 + lamb * lam2 / 2) - lame;
    fc = nuc + sqrt (nuc^2 + lamc * lam2 / 2) - lame;
    if ((fc*fb > 0) && (fa*fc < 0)) 
        lamb=lamc; 
    end
    if ((fc*fa > 0) && (fb*fc < 0)) 
        lama=lamc; 
    end
    k=k+1;
    ra=abs(lama - lamb);
end
du32 = lamc;
end

function du4 = DulKoefTep1AdiCyl(T, m2, lame, d)
fra = 0.9;
lama = 1e4;
lamb = 1e-6;
ep = 1e-7;
%koe = (lae(2) - lae(1)) / (tee(2) - tee(1));
%lame = lae(1) + koe * (T - tee(1));
lam2 = opredTeploprovVozd(T);
k=0; h = d;
while (abs(lama - lamb) > ep)
    lamc = (lama + lamb) / 2;
    h1 = fra * d;
    d1 = d * sqrt(4 * m2 / pi / fra);
    F1 = ( pi / 4 ) * d1^2; 
    R1a = (h - h1) / 2 / F1 / lama;
    R1b = (h - h1) / 2 / F1 / lamb;
    R1c = (h - h1) / 2 / F1 / lamc;
    R2 = (h - h1) / 2 / F1 / lam2;
    F12 = ( pi / 4 ) * (d^2 - d1^2); 
    R3a = h / F12 / lama;
    R3b = h / F12 / lamb;
    R3c = h / F12 / lamc;
    F = (pi / 4) * d^2; 
    R = h / lame / F;
    fa = R3a * (2 * R1a + R2) / (2 * R1a + R2 + R3a) - R;
    fb = R3b * (2 * R1b + R2) / (2 * R1b + R2 + R3b) - R;
    fc = R3c * (2 * R1c + R2) / (2 * R1c + R2 + R3c) - R;
    if (fc * fb > 0) 
        if (fa * fc < 0) 
            lamb=lamc; 
        end
    end
    if (fc * fa > 0) 
        if (fb * fc < 0) 
            lama=lamc; 
        end
    end
    k=k+1;
    if (k > 1e2) 
        break; end;
end
du4 = lamc;
end

function du5 = DulKoefTep1IzotermCyl(T, m2, lame, d)
fra = 0.9;
lama = 1e5;
lamb = 1e-6;
ep = 1e-7;
%koe = (lae(2) - lae(1)) / (tee(2) - tee(1));
%lame = lae(1) + koe * (T - tee(1));
lam2 = opredTeploprovVozd(T);
k=0; h = d;
while (abs(lama - lamb) > ep)
    lamc = (lama + lamb) / 2;
    h1 = fra * d;
    d1 = d * sqrt(4 * m2 / pi / fra);
    F = ( pi / 4 ) * d^2; 
    R1a = (h - h1) / 2 / F / lama;
    R1b = (h - h1) / 2 / F / lamb;
    R1c = (h - h1) / 2 / F / lamc;
    F1 = ( pi / 4 ) * d1^2; 
    R2 = h1 / F1 / lam2;
    F12 = ( pi / 4 ) * (d^2 - d1^2); 
    R3a = h1 / F12 / lama;
    R3b = h1 / F12 / lamb;
    R3c = h1 / F12 / lamc;
    F = ( pi / 4 ) * d^2; 
    R = h / lame / F;
    fa = 2 * R1a + R3a * R2 / (R2 + R3a) - R;
    fb = 2 * R1b + R3b * R2 / (R2 + R3b) - R;
    fc = 2 * R1c + R3c * R2 / (R2 + R3c) - R;
    if (fc * fb > 0) 
        if (fa * fc < 0) 
            lamb=lamc; end; end;
    if (fc * fa > 0) 
        if (fb * fc < 0) 
            lama=lamc; end; end;
    k=k+1;
    if (k > 1e2) 
        break; end;
end
du5 = lamc;
end