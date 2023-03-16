function dul = DulnevKoefTep207Vzaim()
format long g;
Te=arrTem_VVF2();
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
%lae = [0.180183842490065 0.244572803850782 0.127756995443625 0.174128767881838 0.104329016643457 0.214722590630864];
%tee = [620.15 890.65 627.0875 882,4 620.525 905.3375];
%k = 0; po = 0.5575; %8-4 mm
k = 1; por = 0.6635; %2-0.7 mm
%T=1e2+1e0:1e-1:9e2-1e0;
te0 = 273.15;
T = T + te0;
%koe = (lae(2) - lae(1)) / (tee(2) - tee(1));
dul = PoiskLam1(k, por, T, lae);
end

function lam = PoiskLam1(n, m2, T, ktp)
legr = LevGranPorom();
prgr = PravGranPorom();
sr = (legr + prgr) / 2;
if (n == 0)
rp = rasPorpoRazm84();
elseif (n == 1)
rp = rasPorpoRazm207();
end
le = length(sr);
j = 1; vo = 0; pw = 0; mx = 0;
for k=1:le
    if (rp(k) > 0)
    pw(j) = rp(k) * m2;
    vo(j) = sr(k)^3 / pw(j);
    mx = mx + vo(j);
    j = j + 1;
    end
end
l = mx^(1/3)
lb = l / (m2^(1 / 3))
nT=length(T); 
lamvz=0;
for k=1:nT
    lame = ktp(k);
    delta=lb*PoiskC(m2);
    lamvz(k) = DulKoefTep1(T(k), delta, lb, lame); 
end
te0=273.15;
T=T - te0;
%lamvz=lamvz'
%pl=plot(T,lamvz,'-k');
%set(pl,'LineWidth',2); hold on; grid on;
%xlabel({'Температура, °C'}); 
%ylabel({'КТП твердой фазы, Вт/(м*К)'}); 
%title({'График зависимости КТП ТФ от T - Взаимопроникающие компоненты'});
lam = 0;
end

function du1 = DulKoefTep1(T, de, L, lame)
lam1a = 1e5;
lam1b = 0;
ep = 1e-7;
lam2 = opredTeploprovVozd(T);
k=0; 
ra=abs(lam1a-lam1b);
    R3=1/lam2/de;
    R4=L/lam2/((L-de)^2);
    R=1/lame/L;
while ( ra > ep)
    lam1c = (lam1a + lam1b) / 2;
    R1a=L/lam1a/de^2;
    R1b=L/lam1b/de^2;
    R1c=L/lam1c/de^2;
    R2a=1/lam1a/(L-de);
    R2b=1/lam1b/(L-de);
    R2c=1/lam1c/(L-de);
    fa = 1/R1a+2/(R2a+R3)+1/R4-1/R;
    fb = 1/R1b+2/(R2b+R3)+1/R4-1/R;
    fc = 1/R1c+2/(R2c+R3)+1/R4-1/R;
    if (fc*fb > 0) 
        if (fa*fc < 0) 
            lam1b=lam1c; 
        end
    end
    if (fc*fa > 0) 
        if (fb*fc < 0) 
            lam1a=lam1c; 
        end
    end
    k=k+1;
    if (k>1e2) 
        break; 
    end
    ra=abs(lam1a-lam1b);
end
du1 = lam1c;
end

function du2 = PoiskC(m2)
C1a = 1e3;
C1b = 1e-8;
ep = 1e-7;
k=0; ra=1e2;
while ( ra > ep)
    C1c = (C1a + C1b) / 2;
    fa = 2*(C1a^3)-3*(C1a^2)+1-m2;
    fb = 2*(C1b^3)-3*(C1b^2)+1-m2;
    fc = 2*(C1c^3)-3*(C1c^2)+1-m2;
    if (fc*fb > 0) 
        if (fa*fc < 0) 
            C1b=C1c; end; 
    end;
    if (fc*fa > 0) 
        if (fb*fc < 0) 
            C1a=C1c; end;
    end;
    k=k+1;
    if (k>1e2) 
        break; 
    end;
    ra=abs(C1a-C1b);
end
du2 = C1c;
end