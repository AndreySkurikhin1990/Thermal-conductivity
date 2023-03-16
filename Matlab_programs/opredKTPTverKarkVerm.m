function [ oktptks ] = opredKTPTverKarkVerm(vyfv, tem, vysv, vyuv, vmivmf, y0)
te0=273.15; ko=1e-2; wsio=368e-3; walo=132e-3; wmgo=217e-3; kk=1e-3; %qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); 
vybvesch=1; vybmar=0; vpmf=0; vpkf=0; 
por=novNapMas(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf);
if (vyfv==0) %por = 0.6635; %пористость вермикулита фракции 2-0,7 мм
    r = (2+0.7)*kk/2; walo=14*ko; wsio=45.8*ko; wmgo=29*ko;
elseif (vyfv==1) %por=0.5575; %пористость вермикулита фракции 8-4 мм
    r = (8+4)*kk/2; walo=12.9*ko; wsio=46.7*ko; wmgo=28.2*ko;
end
    k=0; tholvs = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf, tem);
	k=k+1; tgorvs = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf, tem);
	k=k+1; qobvs = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf, tem);
	k=k+1; ektpvs = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf, tem);
	k=k+1; tsredvers = napMasEKTPVerNac(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf, k, y0, vyuv, vmivmf, tem);
    cemv=length(tsredvers);
	k=0; tholv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv);
	k=k+1; tgorv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv);
	k=k+1; qobv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv);
	k=k+1; temvs = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv)
	k=k+1; laefm = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv)
lavo=0; dlma=length(temvs); nT=dlma;
for k=1:nT
    lavo(k)=opredTeploprovVozd(temvs(k));
end
lvmi=min(lavo);
k=1; tm=0; tm=rasPorpoRazVer(vyfv,k);
srp=tm(1); %k=2; tm=0; tm=rasPorpoRazVer(vyfv,k); marp=tm(1)
k=3; rp=0; rp=rasPorpoRazVer(vyfv, k);
s=sum(rp);
k=5; sr=0; sr=rasPorpoRazVer(vyfv, k);
k=k+1; legr=0; legr=rasPorpoRazVer(vyfv, k);
k=k+1; prgr=0; prgr=rasPorpoRazVer(vyfv, k);
na=1; cvp=na+21; cve=cvp+7*2; cved=cve+4;
q=1; stchm = DulnevZernVer(por, temvs, laefm, wmgo, wsio, walo, vyfv, lavo, q);
for k=na:cved
    if (k==na)
        q=0; srk = DulnevZernVer(por, temvs, laefm, wmgo, wsio, walo, vyfv, lavo, q);
    elseif ((k<=cvp) && (k>na)) %1..22
        srk = opredTvChaKoeTepSrInSpSha(k-na, por, temvs, te0, laefm, r, lavo);
    elseif ((k<=cve) && (k>cvp)) %23..29, 30..36
        srk=DulnevKoefTepVermN(k-cvp-1, temvs, laefm, srp, abs(cve-cvp)/2, rp, legr, prgr, sr, por, lavo);
    elseif ((k<=cved)  && (k>cve)) %37..40
        srk=DopFunRasTveKarVer(temvs, k-cve-1, por, r, laefm, length(temvs), lavo, wmgo, wsio, walo, stchm);
    end
    for j=1:nT
        sr(k,j)=srk(j);
    end
    %k=k'; srk=srk';
end
q = 1; e=1e-6;
for j=1:cved
    f = 1;
    for k=1:dlma
        if (abs(sr(j,k)<e) && (f >0))
            f = -1; break;
        end
    end
    fl=1; 
    for k=1:dlma 
        if ((sr(j, k)<lvmi) && (fl>0)) 
            fl=-1; break;
        end
    end
	pf=1; 
    for k=2:dlma 
        if ((sr(j, k)<sr(j,k-1)) && (pf>0)) 
            pf=-1; break;
        end
    end
    if ((f > 0) && (fl>0) && (pf>0))
        no(q) = j;
        q = q + 1;
    end
end
q=length(no); 
for k=1:q
    j=no(k);
    srk=0;
        for j=1:dlma
        srk(j)=sr(no(k),j);
        end
    srk=srk;
end
no=no;
q=length(no); e=1e-6;
srzn=zeros(1, nT);
for j=1:nT
s = 0;
for k=1:q
    t=abs(sr(no(k),j));
    if (t>e)
    s = s + t;
    end
end
srzn(j) = s / q;
end
n=length(tem); hf=1e1;
szktptk=zeros(1, n);
for k=1:n
    f=-1;
for j=1:length(tsredvers)
    if (abs(tsredvers(j)-tem(k))<hf)
        f=1; x=j; break;
    end
end
if (f>0)
szktptk(k)=srzn(x);
else
szktptk(k)=0;
end
end
oktptks=szktptk;
end

function [ vyma ] = vydelPol(temvcs, temvhs, qos, ektpvs, temvss, v, n)
tev0 = 273.15; templa = 134.0*1e1 + tev0; q = 0; e=1e-6;
for k = 1:n
    if ((temvcs(k)>e) && (temvhs(k)>e) && (qos(k)>e) && (ektpvs(k)>e) && (temvss(k)>e) && (temvhs(k)<templa)) 
        q=q+1; 
    end
end
qn = q; hf = 1e0; nf = 0;
        for k = 1:qn
            nf = nf + hf;
        end
	q = 1; 
for k=1:n
if ((temvcs(k)>e) && (temvhs(k)>e) && (qos(k)>e) && (ektpvs(k)>e) && (temvss(k)>e) && (temvhs(k)<templa)) 
                switch (v)
                case (0) 
                    vm(q) = temvcs(k);
                case (1) 
                    vm(q) = temvhs(k);
                case (2) 
                    vm(q) = qos(k);
                case (3) 
                    vm(q) = temvss(k);
                case (4) 
                    vm(q) = ektpvs(k);
                end
            q=q+1;
end
end
	vyma=vm;
end

function t = vyvodmassiva(sr, lentem, nomer)
w=nomer;
    for j=1:lentem
        srk(j)=sr(w,j);
    end
    srk=srk';
t=0;
end