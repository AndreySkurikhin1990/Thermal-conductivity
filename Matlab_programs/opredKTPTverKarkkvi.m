function [ oktptks ] = opredKTPTverKarkkvi(vmkvi,tem)
te0=273.15; srp=rasPorpoRazkvi(vmkvi,1); ep=1e-8;
y0=30e-3; %qo=PoiskZavVelTemVer(0); te=PoiskZavVelTemVer(1); 
por=oprporkvi(vmkvi);
lentem=length(tem);
koef=poiskoefktpkvi(vmkvi);
for k=1:lentem
ts=tem(k); laef=koef(2)*ts+koef(1); laefm(k)=laef;
end
x0=srp/(por^(1/3))/2;
for k=1:lentem
    srk(k)=0;
end
dkosc=oprStCherkvi(vmkvi,tem,5,0);
dkosct=oprStCherkvi(vmkvi,tem,6,0);
kusci=oprStCherkvi(vmkvi,tem,3,0);
tkusci=oprStCherkvi(vmkvi,tem,4,0);
stchk=oprStCherkvi(vmkvi,tem,2,0);
wsio=stchk(1); walo=stchk(2); wmgo=stchk(3);
stchk=0; stchk=oprStCherkvi(vmkvi,tem,1,1)
na=1; cvp=na+21; cve=cvp+14; cved=cve+na; f=(cve-cvp)/2; 
s=0; p=0; ep=1e-4;
for k=1:lentem
    lavo(k)=opredTeploprovVozd(tem(k));
    for j=1:cved
        sr(j,k)=0;
    end
end
tena=3e2; dtos=1e2; q=1; j=1; ht=1;
for k=na:cved
if (k==na)
    %k=k+1; continue;
    srk=DulnevZernkvi(por, tem, srk, na, laefm, lavo, wsio, walo, wmgo, x0, lentem, tena, dtos, length(kusci), kusci, tkusci, stchk,x0);
elseif ((k<=cvp) && (k>na))
    %k=k+1; continue;
        srk=opredTvChaKoeTepSrInSpSha(k-na, por, tem, te0, laefm, x0);
elseif ((k<=cve) && (k>cvp))
        srra=rasPorpoRazkvi(vmkvi, 5); 
        prgr=rasPorpoRazkvi(vmkvi, 7); 
        legr=rasPorpoRazkvi(vmkvi, 6); 
        raspr=rasPorpoRazkvi(vmkvi, 3);
        tm=rasPorpoRazkvi(vmkvi, 1); srp=tm(1); %dmi=length(srra);
        srk=DulnevKoefTepVermN(k-cvp-1, por, tem, te0, laefm, srp, 0, f, 1, raspr, legr, prgr, srra);
elseif ((k<=cved) && (k>cve))
        %k=k+1; continue;
        srk=VasFrayObsh(laefm, lavo, por, lentem);
end
for w=1:lentem
sr(j,w)=srk(w); 
end
j=j+1;
end
lvmi=min(lavo);
for k=1:lentem
	p=0; s=0; q=1;
for j=1:cved
        if (abs(sr(j,k)>lvmi)) 
            s=s+sr(j,k); 
            p=p+ht; 
            no(q)=j; 
            q=q+1;
        end
end
        srzn(k)=s/p; 
end
q=length(no);
for k=1:q
    srk=0; w=no(k); disp(w); 
    for j=1:lentem
        srk(j)=sr(w,j);
    end
    %srk=srk'
end
no=no'
if (proverkananuli(srzn,lvmi)>0)
srzn=srzn'
end
%for k=29:30
    %q=vyvodmassiva(sr, lentem, k); 
%end
oktptks=srzn;
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

function p = proverkananuli(m,minim)
f=1;
for k=1:length(m)
    if (abs(m(k))<abs(minim))
    f=-1;
    break;
    end
end
p=f;
end