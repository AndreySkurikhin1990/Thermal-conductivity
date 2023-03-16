function [ oktptks ] = opredKTPTverKarkSha(tem,vysh)
if (vysh == 0)
por = (20+22)*1e-2/2; %пористость шамота
koef=[-0.435e-9,0.685e-6,0.134e-3,0.725]; 
end
if (vysh == 2)
por=(30+33)/1e2/2;
koef=[-0.377e-9,0.918e-6,-0.338e-3,0.77];
end
for k=1:length(tem)
    ts=tem(k);
laef=polyval(koef,ts);
laefm(k)=laef;
end
laefm=laefm'
srp=SreRazPorSha(); x0=srp/(por^(1/3))/2;
te0 = 273.15; tem = tem + te0;
cvp=1+17; cve=cvp+12; cved=cve+1; na = 1;
for k=1:length(tem)
    srk(k)=0;
end
wsi=55.98e-2; wal=39.64e-2; wmg=0.01e-2;
srk = DulnevZernSha(wmg,wal,wsi,por,tem,na,laefm);
f = 1;
for k=1:length(srk)
    if ((srk(k) == 0) && (f > 0))
        f = 0;
    end
end
if (f == 0)
sr = 0;
else
    for j=1:length(srk)
        sr(na,j)=srk(j);
    end
end
na = na + 1;
for k=na:cved
    %disp(k);
        for j=1:length(tem)
            srk(k)=0; 
        end
    if (k<=cvp)
        for j=1:length(tem)
            srk(j)=0;
        end
        if (k>0)
        srk = opredTvChaKoeTepSrInSpSha(k, por, tem, laefm, x0);
        end
    elseif ((k<=cve) && (k>cvp))
        srk = DulnevKoefTepShamot(k - cvp, por, tem, te0, vysh,(cve-cvp)/2);
    elseif ((k<=cved) && (k>cve))
        srk = VasFraySha(k - cvp, por, tem, te0, vysh);
    end
    for j=1:length(tem)
        sr(k,j)=srk(j);
    end
    %disp(sr(k,:));
end
q = 1; no = 0;
for j=1:cve
    f = 1;
    for k=1:length(tem)
        if ((sr(j,k) == 0) && (f > 0))
            f = 0;
        end
    end
    if (f > 0)
        no(q) = j;
        q = q + 1;
    end
end
q=length(no);
for j=1:length(tem)
s = 0;
for k=1:q
    s = s + sr(no(k),j);
end
srzn(j) = s / q;
end
for j=1:length(no)
srk=0;
    for k=1:length(tem)
    srk(k)=sr(no(j),k);
    end
    k=no(j)
srk=srk'
end
oktptks=srzn;
end

function [ t ] = VasFraySha(tem,n)
te=[1e2,18e1,22e1,3e2,4e2,6e2,8e2,1e3,12e2]; 
lsk=[1.52,2.25,2.5,2.85,3.24,3.28,3.46,4.06,4.98]; 
ke=length(tem);
for k=1:length(te)
    vfs(k)=opredKTPTKTochSha(lsk, te, tem(k), ke);
end
t=vfs;
end

function srp = SreRazPorSha()
rpr(1) = 0;
raspr=[21.8,21.75,21.625,21.25,20.75,20.25,19.75,19,18.25,17.5,12,9.875,8.75,8.25,7.25,6.25,5,4.25,3.25,0.875,0.325,0];
for k=2:length(raspr)
    rpr(k)=(raspr(k-1)-raspr(k))*1e2/raspr(1);
end
prgr=[0.05,0.1,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7,8,9,10,20,30,40,50]*1e-6; %в метрах
legr(1)=0;
for k=2:length(prgr)
    legr(k)=prgr(k-1);
end
sr = (legr + prgr) / 2; %средний размер пор
le = length(sr);
srpo = 0;
for k=1:le
    if (rpr(k) > 0)
    srpo = srpo + rpr(k) * sr(k) / 1e2; %объемная доля поры заданного размера в полном (во всем) объеме, в процентах
    end
end
srp = srpo;
end