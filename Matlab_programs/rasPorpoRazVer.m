function [ vm ] = rasPorpoRazVer(vfv, ide)
format longg;
rpn=1e0; de=1e0; h=1e-6;
if (vfv==0)
raspr0=rasPorpoRazm207(); 
elseif (vfv==1)
    raspr0=rasPorpoRazm84(); 
end
prgr0=PravGranPorom(); 
legr0=LevGranPorom(); 
ma=max(prgr0);
if (ma<1e1)
    ko=1e6; 
else
    ko=1e0;
end
prgr0=ko*prgr0;
legr0=ko*legr0;
k=0; prgr=sortirovka(prgr0, legr0, raspr0, k)';
k=k+1; legr=sortirovka(legr0, legr0, raspr0, k)';
k=k+1; raspr=sortirovka(prgr0, legr0, raspr0, k)';
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
    case 1
    k=1; mu(k)=s; %средний размер пор, м
    case 2
    k=1; mu(k)=r; %максимальный размер  пор, м
    case 3
    r=sum(raspr1);
    raspr1=raspr1/r;
    mu=raspr1'; %распределение пор, %
    case 4
    mu=rp'; %распределение пор, см3
    case 5
    mu=srra'; %средний размер (массив)
    case 6
    mu=legr1'; %левая граница
    case 7
    mu=prgr1'; %правая граница
    otherwise
    mu=[0];
end
vm=mu';
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