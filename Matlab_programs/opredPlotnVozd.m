function opr = opredPlotnVozd(T)
te0=273.15;
te = arrTempAir()+te0;
lamb = arrPlotnAir();
n = length(te);
f = 1; p=0;
for k=1:n
    if ((te(k)>T) && (f > 0))
        p = k;
        f = 0;
    end
        if ((te(k)==T) && (f > 0))
        p = k;
        f = 0;
        end
end
lam = 0;
if (f == 0)
if (p > 1)
lam = lamb(p-1) + (lamb(p) - lamb(p-1)) * (T - te(p-1)) / (te(p) - te(p-1));
end
if (p==1)
    ko = (lamb(2) - lamb(1))/(te(2) - te(1));
    lam = lamb(1) + ko * (T - te(1));
end
if (p==0)
    ko = (lamb(n) - lamb(n - 1))/(te(n) - te(n - 1));
    lam = lamb(n) + ko * (T - te(n));
end
end
if (f>0)
    ko = (lamb(n) - lamb(n - 1))/(te(n) - te(n - 1));
    lam=lamb(n) + ko * (T - te(n));
end
opr = lam;
end