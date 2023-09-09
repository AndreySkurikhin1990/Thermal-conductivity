function [ oktptks ] = opredKTPTverKarkVerm(vyfv, tem, vysv, vyuv, vmivmf, y0)
te0=273.15; ko=1e-2; kk=1e-3; wsio=368*kk; walo=132*kk; wmgo=217*kk; %qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); 
vybvesch=1; vybmar=0; vpmf=0; vpkf=0; hf=1e0;
por=novNapMas(vybvesch, vybmar, vyfv, vysv, vpmf, vpkf);
if (vyfv==0) %por = 0.6635; %пористость вермикулита фракции 2-0,7 мм
    r = (2+0.7)*kk/2; walo=14*ko; wsio=45.8*ko; wmgo=29*ko;
elseif (vyfv==1) %por=0.5575; %пористость вермикулита фракции 8-4 мм
    r = (8+4)*kk/2; walo=12.9*ko; wsio=46.7*ko; wmgo=28.2*ko;
end
    k=0; tholvs = napMasEKTPVerNac(vyfv, vysv, k, y0, vyuv, vmivmf, tem);
	k=k+1; tgorvs = napMasEKTPVerNac(vyfv, vysv, k, y0, vyuv, vmivmf, tem);
	k=k+1; qobvs = napMasEKTPVerNac(vyfv, vysv, k, y0, vyuv, vmivmf, tem);
	k=k+1; ektpvs = napMasEKTPVerNac(vyfv, vysv, k, y0, vyuv, vmivmf, tem);
	k=k+1; tsredvers = napMasEKTPVerNac(vyfv, vysv, k, y0, vyuv, vmivmf, tem);
    cemv=length(tsredvers);
    tholv = tholvs;
	tgorv = tgorvs;
	qobv = qobvs;
	temvs = tsredvers;
	laefm = ektpvs
	%k=0; tholv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv);
	%k=k+1; tgorv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv);
	%k=k+1; qobv = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv);
	%k=k+1; temvs = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv);
	%k=k+1; laefm = vydelPol(tholvs, tgorvs, qobvs, ektpvs, tsredvers, k, cemv);
lavo=0; dlma=length(temvs); nT=dlma;
for k=1:nT
    lavo(k)=opredTeploprovVozd(temvs(k));
end
lvmi=min(lavo);
k=1; tm=0; tm=rasPorpoRazVer(vyfv,k);
srp=tm(1) %k=2; tm=0; tm=rasPorpoRazVer(vyfv,k); marp=tm(1)
k=3; rp=0; rp=rasPorpoRazVer(vyfv, k);
s=sum(rp);
k=5; sr=0; sr=rasPorpoRazVer(vyfv, k);
k=k+1; legr=0; legr=rasPorpoRazVer(vyfv, k);
k=k+1; prgr=0; prgr=rasPorpoRazVer(vyfv, k);
na=1; cvp=na+21; cve=cvp+7*2; cved=cve+4;
q=1; stchm = DulnevZernVer(por, temvs, laefm, wmgo, wsio, walo, vyfv, lavo, q);
e=1e-8;
for k=na:cved
    if (k==na)
        j=0; srk = DulnevZernVer(por, temvs, laefm, wmgo, wsio, walo, vyfv, lavo, j);
    elseif ((k<=cvp) && (k>na)) %1..22
        srk = opredTvChaKoeTepSrInSpSha(k-na, por, temvs, laefm, r, lavo);
    elseif ((k<=cve) && (k>cvp)) %23..29, 30..36
        srk=DulnevKoefTepVermN(k-cvp-1, temvs, laefm, srp, abs(cve-cvp)/2, rp, legr, prgr, sr, por, lavo);
    elseif ((k<=cved)  && (k>cve)) %37..40
        srk=DopFunRasTveKarVer(temvs, k-cve-1, por, r, laefm, length(temvs), lavo, stchm);
    end
        srk=proverka(srk, laefm);
    if (sum(srk)>e)
        k=k'
        srk=srk
        no(q)=k;
        for j=1:length(srk)
        srktp(q,j)=srk(j);
        end
        q=q+1;
    end
end
srktp=srktp
n=length(temvs);
srk=zeros(1,n);
no=no
for j=1:n
    s=0; t=0;
for k=1:length(no)
    if (srktp(k,j)>e)
        s=s+srktp(k,j);
        t=t+hf;
    end
end
if (t>e)
srk(j)=s/t;
end
end
srk=srk
oktptks=srk;
end

function [ vm ] = proverka(srk, laefm)
n=length(srk);
f=1;
e=1e-8;
for k=1:n
    r=srk(k);
    t=laefm(k);
    if (r<t)
        f=-1;
    end
    if (k>1)
        p=srk(k-1);
    if ((abs(p-r)<e) || (p>r))
        f=-1;
    end
    end
    if (f<0)
        break;
    end
end
if (f<0)
    for k=1:n
        srk(k)=0;
    end
end
vm=srk;
end