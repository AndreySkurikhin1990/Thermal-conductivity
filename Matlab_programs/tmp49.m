%рассчитывает ослабление интегральной степени черноты для SiO2, Al2O3
function [ t ] = tmp49(vy,te)
format long g;
n=length(te);
scma=[1,0.940,0.870,0.801,0.736,0.676,0.635,0.590,0.567,0.543,0.530,0.525,0.515,0.507]; %степень черноты магнезита
mdma=[0.98,0.02]; %MgO - магнезит
%-------
scsh=[1,0.976,0.949,0.905,0.859,0.812,0.774,0.737,0.709,0.681,0.661,0.639,0.626,0.620]; %степень черноты шамота
scksh=[1,0.980,0.951,0.920,0.883,0.853,0.821,0.790,0.767,0.746,0.730,0.715,0.705,0.692]; %степень черноты корундошамота
scktik=[1,0.983,0.936,0.867,0.819,0.721,0.659,0.593,0.541,0.490,0.453,0.429,0.403,0.384]; %степень черноты каолинового теплоизоляционного кирпича
scmu=[1,0.984,0.941,0.882,0.813,0.751,0.695,0.641,0.594,0.558,0.530,0.499,0.479,0.462]; %степень черноты муллита
sckr=[1,0.984,0.953,0.917,0.854,0.808,0.756,0.711,0.578,0.523,0.495,0.468,0.448,0.429]; %степень черноты кремнезема
%-------
mdsh=[0.56,0.396]; mdsh=[mdsh(1)/(mdsh(1)+mdsh(2)),mdsh(2)/(mdsh(1)+mdsh(2))]; mdsh=mdsh'; %SiO2, Al2O3 - шамот
mdksh=[0.28,0.70]; mdksh=[mdksh(1)/(mdksh(1)+mdksh(2)),mdksh(2)/(mdksh(1)+mdksh(2))]; mdksh=mdksh'; %SiO2, Al2O3 - корундошамот
mdktik=[0.57,0.4]; mdktik=[mdktik(1)/(mdktik(1)+mdktik(2)),mdktik(2)/(mdktik(1)+mdktik(2))]; mdktik=mdktik'; %SiO2, Al2O3 - каол. т/и кирпич
mdmu=[0.28,0.72]; %SiO2, Al2O3 - муллит
mdkr=[0.985,0.01]; mdkr=[mdkr(1)/(mdkr(1)+mdkr(2)),mdkr(2)/(mdkr(1)+mdkr(2))]; mdkr=mdkr'; %SiO2, Al2O3 - кремнезем
for j=1:n
kor=0; k=1;
kor=RasKor(mdsh(1),mdsh(2),mdksh(1),mdksh(2),scsh(j),scksh(j),kor,k); k=k+1;
kor=RasKor(mdsh(1),mdsh(2),mdktik(1),mdktik(2),scsh(j),scktik(j),kor,k); k=k+1;
kor=RasKor(mdsh(1),mdsh(2),mdmu(1),mdmu(2),scsh(j),scmu(j),kor,k); k=k+1;
kor=RasKor(mdksh(1),mdksh(2),mdktik(1),mdktik(2),scksh(j),scktik(j),kor,k); k=k+1;
kor=RasKor(mdksh(1),mdksh(2),mdmu(1),mdmu(2),scksh(j),scmu(j),kor,k); k=k+1;
kor=RasKor(mdktik(1),mdktik(2),mdmu(1),mdmu(2),scktik(j),scmu(j),kor,k); k=k+1;
kor=RasKor(mdsh(1),mdsh(2),mdkr(1),mdkr(2),scsh(j),sckr(j),kor,k); k=k+1;
kor=RasKor(mdksh(1),mdksh(2),mdkr(1),mdkr(2),scksh(j),sckr(j),kor,k); k=k+1;
kor=RasKor(mdktik(1),mdktik(2),mdkr(1),mdkr(2),scktik(j),sckr(j),kor,k); k=k+1;
kor=RasKor(mdmu(1),mdmu(2),mdkr(1),mdkr(2),scmu(j),sckr(j),kor,k); k=k+1;
usr=UsredMas(kor,k);
msio(j)=usr(1);
malo(j)=usr(2);
end
msio=msio';
malo=malo';
for k=1:length(scma)
    scma(k)=(scma(k)-msio(k)*mdma(2))/mdma(1);
end
scma=scma';
if (vy==0)
    t=scma;
elseif (vy==1)
    t=msio;
elseif (vy==2)
    t=malo;
end
end

function [ ko ] = RasKor(a11,a12,a21,a22, b1,b2,kor,k)
mas(1,1)=a11; mas(1,2)=a12; 
mas(2,1)=a21; mas(2,2)=a22; 
b(1)=b1; b(2)=b2; 
x=inv(mas)*b';
if ((x(1)>=0) && (x(2)>=0) && (x(1)<=1) && (x(2)<=1))
    kor(k,1)=x(1); kor(k,2)=x(2);
else
    kor(k,1)=0; kor(k,2)=0;
end
ko=kor;
end

function [ usr ] = UsredMas(kor,k)
s1=0; s2=0; q=0;
for j=1:k-1
    if ((kor(j,1)>0) && (kor(j,2)>0))
        s1=s1+kor(j,1); 
        s2=s2+kor(j,2); 
        q=q+1;
    end
end
s1=s1/q;
s2=s2/q;
usr=[s1,s2];
end

function opr = opredStSch(te,lamb,T)
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

function t = tmp49_1()
T=573.15:1e2:1373.15;
for k=1:length(T)
    msiom(k)=opredStSch(te,msio,T(k));
    mscma(k)=opredStSch(te,scma,T(k));
    malom(k)=opredStSch(te,malo,T(k));
end
msiom=msiom';
mscma=mscma';
malom=malom';
t=0;
end