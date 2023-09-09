function [ vm ] = rasPorpoRazITOM(vybitom, ide)
format longg;
rpn=1e0; de=1e0; h=1e-6;
k=0; hps=RazPorItom(vybitom, k);
k=1; srp=hps(k);
k=3; por=hps(k);
k=1; raspr0=RazPorItom(vybitom, k);
k=3; prgr0=RazPorItom(vybitom, k);
k=4; legr0=RazPorItom(vybitom, k);
ma=max(prgr0);
if (ma<1e1)
    ko=1e6; 
else
    ko=1e0;
end
prgr0=ko*prgr0;
legr0=ko*legr0;
k=0; prgr=sortirovka(prgr0, legr0, raspr0, k);
k=k+1; legr=sortirovka(legr0, legr0, raspr0, k);
k=k+1; raspr=sortirovka(prgr0, legr0, raspr0, k);
k=0; raspr1=rasrasporporaz0(prgr, legr, raspr, rpn, k, de);
k=k+1; prgr1=rasrasporporaz0(prgr, legr, raspr, rpn, k, de);
rprgr=length(prgr1);
ma=max(prgr1);
if (rprgr<ma)
    h=prgr1(rprgr)-prgr1(rprgr-1);
    for k=rprgr:ma
        raspr1(k)=0;
        prgr1(k)=prgr1(k-1)+h;
    end
end
k=1; legr1(k)=0;
for k=2:rprgr
    legr1(k)=prgr1(k-1);
end
%k=k+1; legr1=rasrasporporaz0(ko*prgr, ko*legr, raspr, rpn, k, de);
ko=1e0;
if (max(prgr1)>1e1)
    ko=1e-6;
end
prgr1=ko*prgr1;
legr1=ko*legr1;
qg=length(legr1);
r=0; s=0; sp=0;
k=1; legr1(k)=0;
for k=1:qg
    p=prgr1(k); 
    l=legr1(k);
    h=p-l;
    srra(k)=(p+l)/2e0; 
    s=s+srra(k)*raspr1(k); 
    sp=sp+raspr1(k);
    r=r+h; 
end
s=s/sp;
switch (ide)
    case 0
        mu=hps;
    case 1
        k=1; mu(k)=s; %средний размер пор, м
    case 2
        k=1; mu(k)=r; %максимальный размер  пор, м
    case 3
        r=sum(raspr1);
        raspr1=raspr1/r;
        mu=raspr1; %распределение пор, %
    case 4
        mu=rp; %распределение пор, см3
    case 5
        mu=srra; %средний размер (массив)
    case 6
        mu=legr1; %лева€ граница
    case 7
        mu=prgr1; %права€ граница
    otherwise
        mu=[0];
end
vm=mu;
end

function [ t ] = rasrasporporaz0(prgr00, legr00, raspr0, rpn, vyb, de)
ht=1e0; e=ht/1e1; e1=rpn+e; n=length(prgr00);
s=0; 
for k=1:n
    if (prgr00(k)<e1) 
        p=k; s=s+raspr0(k); %размер пор до 1 мкм
    end
end
k=1; 
r=legr00(k);
g=prgr00(k);
x=raspr0(k);
k=k+1;
t=raspr0(k);
if (r>0)
    ko=(t-x)/(g-r);
    s=s+r*ko;
end
s=s';
q=1; rprm(q)=s; q=q+1;
for k=p:n 
    r=prgr00(k)-legr00(k);
    if ((prgr00(k)>e1) && (r<e1))
		rprm(q)=raspr0(k);
		q=q+1;
    end
    if (r>e1)
        m=raspr0(k)/r;
		t=r; j=0; 
        while (t>e) 
            j=j+1; 
            t=t-ht;
        end
            jk=j;
        for j=1:jk
		rprm(q)=m;
		q=q+1; 
        end
    end
    n=q;
end
prgr01=rpn; k=1; legrm(k)=min(legr00); qg=length(rprm);
for k=1:qg
    prgrm(k)=prgr01;
if (k>1) 
    legrm(k)=prgrm(k-1);
end
    ras01(k)=rprm(k);
    prgr01=prgr01+de;
end
switch (vyb)
    case (0) 
    t=ras01;
    case (1) 
    t=prgrm; 
    case (2) 
    t=legrm;
    otherwise t=[0];
end
end

function [ vm ] = sortirovka(prgr, legr, raspr, vy)
q=length(prgr);
while (q>1)
k=1; x=k; mp=prgr(k); ml=legr(k); mr=raspr(k);
for k=2:q
if (mp<prgr(k))
    mp=prgr(k);
    ml=legr(k);
    mr=raspr(k);
    x=k;
end
end
t=prgr(q);
g=legr(q);
s=raspr(q);
prgr(q)=mp;
legr(q)=ml;
raspr(q)=mr;
prgr(x)=t;
legr(x)=g;
raspr(x)=s;
q=q-1;
end
switch (vy)
    case (0)
        vm=prgr;
    case (1)
        vm=legr;
    case (2)
        vm=raspr;
    otherwise
        vm=[0];
end
end

function [ rt ] = RazPorItom(no, vy) %расчет размера пор %выбор номера »“ќћ
ko=1e-2; po=[(80+82)/2,(75+78)/2,(65+68)/2,(62+65)/2]*ko; %пористость
kk=1e-6; rapo=[150,95,85,75,65,55,45,35,25,15,7.5,4,2,0.75,0.3,0.055,0.0055]*kk; %в метрах
por=po(no+1);
switch no
    case (0)
    rpr = [0.57,0.3,0.6,0.87,1.35,2.07,3.72,3.81,5.38,7.6,9.67,10.87,34.68,10.78,6.25,1.47,0]*ko; %распределение частиц по размерам - 17
    case (1)
    rpr = [0.26,0.15,0.26,0.22,0.81,0.81,1.88,3.95,5.54,7.35,7.09,9.01,34.9,13.59,7.5,2.14,4.54]*ko;
    case (2)
    rpr = [0.4,0.09,0.44,0.22,0.66,1.02,1.33,2.66,4.07,10.71,12.17,11.29,35.06,11.24,7.13,1.51,0]*ko;
    case (3)
    rpr = [0.23,0.19,0.04,0.61,0.23,1.03,0.8,2.47,5.66,10.87,14.18,12.5,32.61,11.59,5.25,1.75,0]*ko;
end
rapo=fliplr(rapo);
rpr=fliplr(rpr);
k=1; legr(k)=0.0;
n=length(rapo);
for k=1:(n-1)
prgr(k)=(rapo(k)+rapo(k+1))/2;
end
for k=2:n
legr(k)=prgr(k-1);
end
srp=0;
for k=1:length(rpr)
    srp = srp + rpr(k) * rapo(k);
end
srp=srp/sum(rpr);
srts = (1 - por) * srp / por; %средний размер твердого скелета
switch (vy)
    case (0)
        rt=[srp,srts,po];
    case (1)
        rt=rpr;
    case (2)
        rt=rapo;
    case (3)
        rt=prgr;
    case (4)
        rt=legr;
end
end