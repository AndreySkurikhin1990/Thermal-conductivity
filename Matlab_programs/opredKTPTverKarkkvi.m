function [ oktptks ] = opredKTPTverKarkkvi(vmkvi, tem)
y0=3e1*1e-3; srp=rasPorpoRazkvi(vmkvi,1); e=1e-8; hf=1e0; %qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); 
por=oprporkvi(vmkvi);
lentem=length(tem);
koef=poiskoefktpkvi(vmkvi);
laefm=zeros(1,lentem);
w=1; x=w+1;
for k=1:lentem
ts=tem(k); laef=koef(x)*ts+koef(w); laefm(k)=laef;
end
laefm=laefm;
    k=0; tholv = napMasEKTPkviNac(vmkvi, k, y0, tem);
	k=k+1; tgorv = napMasEKTPkviNac(vmkvi, k, y0, tem)
	k=k+1; qobv = napMasEKTPkviNac(vmkvi, k, y0, tem);
	k=k+1; laefm = napMasEKTPkviNac(vmkvi, k, y0, tem)
	k=k+1; temvs = napMasEKTPkviNac(vmkvi, k, y0, tem)
    tem=temvs;
    lentem=length(tem);
x0=srp/(por^(1/3))/2;
srk=zeros(1,lentem);
k=4; q=0; tkusci=oprStCherkvi(vmkvi,tem,k,q); %4
k=k-1; kusci=oprStCherkvi(vmkvi,tem,k,q); %3
k=k-1; stchk=oprStCherkvi(vmkvi,tem,2,q); %2
j=1; wsio=stchk(j); j=j+1; walo=stchk(j); j=j+1; wmgo=stchk(j);
%stchk=stepCherMas(vmkvi); 
k=k-1; q=q+1; stchk=oprStCherkvi(vmkvi,tem,k,q) %1
na=1; cvp=na+21; cve=cvp+14; cved=cve+na; f=(cve-cvp)/2; 
lavo=zeros(1,lentem);
sr=zeros(lentem,cved);
for k=1:lentem
    lavo(k)=opredTeploprovVozd(tem(k));
end
tena=3e2; dtos=1e2; q=1;
for k=na:cved
if (k==na)
    %k=k+1; continue;
    srk=DulnevZernkvi(por, tem, srk, na, laefm, lavo, wsio, walo, wmgo, x0, lentem, tena, dtos, length(kusci), kusci, tkusci, stchk,x0);
elseif ((k<=cvp) && (k>na))
    %k=k+1; continue;
        srk=opredTvChaKoeTepSrInSpSha(k-na, por, tem, laefm, x0, lavo);
elseif ((k<=cve) && (k>cvp))
        srra=rasPorpoRazkvi(vmkvi, 5); 
        prgr=rasPorpoRazkvi(vmkvi, 7); 
        legr=rasPorpoRazkvi(vmkvi, 6); 
        raspr=rasPorpoRazkvi(vmkvi, 3);
        tm=rasPorpoRazkvi(vmkvi, 1); w=1; srp=tm(w); %dmi=length(srra);
        srk=DulnevKoefTepVermN(k-cvp-1, tem, laefm, srp, f, raspr, legr, prgr, srra, por, lavo);
elseif ((k<=cved) && (k>cve))
        %k=k+1; continue;
        srk=VasFrayObsh(laefm, lavo, por);
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
n=length(tem);
no=no
%qm=size(srktp);
%q=qm(1);
%lentem=qm(2);
%for k=1:q
    %for j=1:lentem
        %srk(j)=srktp(k,j);
    %end
%srk=neBoleeTempMaxPrim(vmkvi, tgorv, srk)
%for j=1:lentem
        %srktp(k,j)=srk(j);
    %end
%end
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

function vm = poiskoefktpkvi(v)
tek0=273.15; t1=25e0+tek0; t2=5e2+tek0; dt=t2-t1;
if (v==4)
kktp(2) = 0.00016; kktp(1) = 0.14;
ktp1=kktp(2)*(t1-tek0)+kktp(1);
ktp2=kktp(2)*(t2-tek0)+kktp(1);
kn=(ktp2-ktp1)/dt; kktp(2) = kn; 
kn=ktp2-kn*t2; kktp(1) = kn;
elseif (v==5)
kn=(0.165-0.105)/dt; kktp(2)=kn; 
kn=kn*t2; kktp(1)=0.165-kn; %выбор
elseif (v==6)
kn=(0.196-0.12)/dt; kktp(2) = kn; 
kn=kn*t2; kktp(1)=0.196-kn; %Двухслойные %выбор
elseif (v==7)
kn=(0.251-0.16)/dt; kktp(2) = kn; 
kn=kn*t2; kktp(1) = 0.251-kn; %Двухслойные %выбор
elseif (v==8)
kn=(0.226-0.16)/dt; kktp(2) = kn; 
kn=kn*t2; kktp(1) = 0.226-kn; %выбор
elseif (v==9)
kn=(0.23-0.195)/dt; kktp(2) = kn; 
kn=kn*t2; kktp(1) = 0.23-kn; %выбор
elseif (v==10)
kn=(0.287-0.25)/dt; kktp(2) = kn; 
kn=kn*t2; kktp(1) = 0.287-kn; %выбор
end
vm=kktp;
end

function por = oprporkvi(vmkvi)
if (vmkvi==4) 
    por=51.74*1e-2;
elseif (vmkvi==5)
    por=52.07*1e-2;
elseif (vmkvi==6)
    por=51.51*1e-2;
elseif (vmkvi==7)
    por=39.75*1e-2;
elseif (vmkvi==8)
    por=40.85*1e-2;
elseif (vmkvi==9)
    por=39.37*1e-2;
elseif (vmkvi==10)
    por=36.07*1e-2;
else
    disp('Net takoy marki KVI');
    por=0;
end
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

function [ sc ] = stepCherMas(vm)
switch (vm)
    case (4)
        sc=[0.815850842738303
         0.751935026308916
         0.690398216392302
         0.630163621003622
         0.585433185393359
         0.544917270569165
          0.50458367532262
         0.463026075798077
         0.432368627097673
           0.4090144726616
         0.352125381721003];
    case (5)
        sc=[0.841017456391689
         0.784574321266392
         0.727941711607636
         0.664273314949933
         0.614124317530928
         0.569812494469711
         0.511602150092411
         0.456866097555998
          0.41932947953891
          0.39015678678099
         0.343017080430017];
     case (6)
         sc=[0.814038134792372
         0.749805727617492
         0.688162233432253
         0.628814696581489
         0.584950064071956
         0.545042113975324
         0.507979433225087
         0.468553984501352
         0.438636494209279
         0.415601118373094
         0.355787643005069];
    case (7)
        sc=[0.847739114711142
         0.787493743149225
         0.727399894992161
         0.660676169923482
         0.612528598638054
         0.571613735430857
         0.514741539895165
         0.460186231436014
         0.422480824310932
         0.392949782733817
           0.3447203546461];
    case (8)
        sc=[0.848294157676158
         0.788234557319249
          0.72825338692151
         0.661383330038131
         0.613055772560333
         0.572016998265028
          0.51448695473577
         0.459474289218746
         0.421566580226359
         0.391892947658496
         0.344141665104063];
    case (9)
        sc=[0.84833320745557
         0.788319515106167
         0.728404735799144
          0.66179040353002
         0.613666857283505
         0.572773309486096
          0.51610293486519
         0.461609773548604
         0.423822867029093
         0.394129367552641
         0.345440976971307];
    case (10)
     sc=[0.848974523047907
         0.789168484275071
         0.729397993484367
         0.662763462068472
          0.61214319661852
         0.568273747030334
         0.513080598115416
         0.461307818621303
         0.424096756846662
         0.394246558415465
         0.345552952553887];
end
end