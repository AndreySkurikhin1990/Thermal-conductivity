%определяет показатель поглощения вермикулита
function t = tmp41()
TrRe1=RasFra60(); 
TrRe2=RasFra60_100(); 
TrRe3=RasFra100_150(); 
TrRe4=RasFra150_200(); 
dlv=1e6*dlivoln();
n=length(TrRe1);
for k=1:n
    al(k)=(TrRe1(k)+TrRe2(k)+TrRe3(k)+TrRe4(k));
    al(k)=al(k)/4;
    pp(k)=al(k)*dlv(k);
end
nadlvo=3;
kodlvo=10;
no=poisIndex(nadlvo,kodlvo,dlv);
kmdv=vydPodmas(no,dlv);
ppn=vydPodmas(no,pp);
n=length(kmdv);
d=(kmdv(n)-kmdv(1))
t=trapz(kmdv,ppn)/d
p=plot(dlv,pp,'-b');
set(p,'LineWidth',1); hold on; grid on; 
xlabel({'Длина волны, мкм'}); 
ylabel({'КП, мкм-1'}); 
title({'График зависимости КП от ДВ'});
end

function [ ras60 ] = RasFra60()
mkbr=[250.239 249.740 249.223 250.395 250.336 249.55]; 
mv=[0.283 0.464 0.812 0.22 0.547 0.777]; 
tol=[0.73 0.72 0.72 0.72 0.71 0.7];
rokbr=2.75; rov=0.49; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %в см3
vv(k)=mv(k)/(1e3*rov); %в см3
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %в мкм
end
ras60=SredGrafRass(1,xv);
end

function [ ras60_100 ] = RasFra60_100()
mkbr=[250 249.629 249.294 249.706 249.510 249.307 250.328 249.604 249.206]; 
mv=[0.255 0.539 0.809 0.295 0.517 0.756 0.36 0.534 0.843]; 
tol=[0.72 0.71 0.7 0.7 0.73 0.72 0.74 0.7 0.76]; 
rokbr=2.75; rov=0.52; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
end
ras60_100=SredGrafRass(2,xv);
end

function [ ras100_150 ] = RasFra100_150()
mkbr=[249.913 249.607 249.218 249.929 249.695 249.306 250.405 249.625 249.348];
mv=[0.315 0.473 0.709 0.293 0.528 0.83 0.27 0.493 0.764]; 
tol=[0.74 0.74 0.72 0.72 0.71 0.7 0.78 0.73 0.76]; 
rokbr=2.75; rov=0.53; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
end
ras100_150=SredGrafRass(3,xv);
end

function [ ras150_200 ] = RasFra150_200()
mkbr=[250.882 249.590 249.213 250.299 249.441 249.365];
mv=[0.320 0.533 0.849 0.223 0.502 0.797]; 
tol=[0.76 0.72 0.69 0.73 0.73 0.73]; 
rokbr=2.75; rov=0.56; n=length(mv);
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
end
ras150_200=SredGrafRass(4,xv);
end

function [ alsr ] = SredGrafRass(vy,xv)
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
switch vy
    case 1
Tal=privedkEdiPropus(TrkoVer5311()); Tal2=privedkEdiPropus(TrkoVer5312());Tal3=privedkEdiPropus(TrkoVer5313()); 
Tal4=privedkEdiPropus(TrkoVer53101()); Tal5=privedkEdiPropus(TrkoVer53102()); Tal6=privedkEdiPropus(TrkoVer53103()); 
alsre=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv);
    case 4
Tal=privedkEdiPropus(TrkoVer5681()); Tal2=privedkEdiPropus(TrkoVer5682()); Tal3=privedkEdiPropus(TrkoVer5683()); 
Tal4=privedkEdiPropus(TrkoVer5691()); Tal5=privedkEdiPropus(TrkoVer5692()); Tal6=privedkEdiPropus(TrkoVer5693());
alsre=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv);
    case 2
Tal=privedkEdiPropus(TrkoVer5421()); Tal2=privedkEdiPropus(TrkoVer5422()); Tal3=privedkEdiPropus(TrkoVer5423()); 
Tal4=privedkEdiPropus(TrkoVer5431()); Tal5=privedkEdiPropus(TrkoVer5432()); Tal6=privedkEdiPropus(TrkoVer5433());
Tal7=privedkEdiPropus(TrkoVer5441()); Tal8=privedkEdiPropus(TrkoVer5442()); Tal9=privedkEdiPropus(TrkoVer5443());
alsre=opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv);
    case 3
Tal=privedkEdiPropus(TrkoVer5551()); Tal2=privedkEdiPropus(TrkoVer5552()); Tal3=privedkEdiPropus(TrkoVer5553()); 
Tal4=privedkEdiPropus(TrkoVer5561()); Tal5=privedkEdiPropus(TrkoVer5562()); Tal6=privedkEdiPropus(TrkoVer5563());
Tal7=privedkEdiPropus(TrkoVer5571()); Tal8=privedkEdiPropus(TrkoVer5572()); Tal9=privedkEdiPropus(TrkoVer5573());
alsre=opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv);
end
alsr=alsre;
end

function [ prked ] = privedkEdiPropus(ar)
arr=ar;
p=length(arr);
for k=1:p
    if (arr(k)>1)
        arr(k)=1;
    end
    if (arr(k)<0)
        arr(k)=0;
    end
end
prked=arr;
end

function [ al ] = opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv)
nTal=length(Tal);
for k=1:nTal 
    Ta(k)=-log(Tal(k))/xv(1); Ta2(k)=-log(Tal2(k))/xv(2); Ta3(k)=-log(Tal3(k))/xv(3);
    Ta4(k)=-log(Tal4(k))/xv(4); Ta5(k)=-log(Tal5(k))/xv(5); Ta6(k)=-log(Tal6(k))/xv(6); 
end
alsr1=opredAlphSred2(Ta,Ta4);
alsr2=opredAlphSred2(Ta2,Ta5);
alsr3=opredAlphSred2(Ta3,Ta6);
al=(alsr1+alsr2+alsr3)/3;
end

function [ al ] =opredAlphSred2(Tal,Tal2)
alsr=0; nTal=length(Tal);
for k=1:nTal 
    alsr(k)=Tal(k)+Tal2(k);
end
al=alsr/2;
end

function [ al ] = opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv)
nTal=length(Tal);
for k=1:nTal
    Ta(k)=-log(Tal(k))/xv(1); Ta2(k)=-log(Tal2(k))/xv(2); Ta3(k)=-log(Tal3(k))/xv(3);
    Ta4(k)=-log(Tal4(k))/xv(4); Ta5(k)=-log(Tal5(k))/xv(5); Ta6(k)=-log(Tal6(k))/xv(6);
    Ta7(k)=-log(Tal7(k))/xv(7); Ta8(k)=-log(Tal8(k))/xv(8); Ta9(k)=-log(Tal9(k))/xv(9); 
end
alsr1=opredAlphSred3(Ta,Ta4,Ta7);
alsr2=opredAlphSred3(Ta2,Ta5,Ta8);
alsr3=opredAlphSred3(Ta3,Ta6,Ta9);
al=(alsr1+alsr2+alsr3)/3;
end

function [ al ] = opredAlphSred3(Tal,Tal2,Tal3)
alsr=0; nTal=length(Tal);
for k=1:nTal 
    alsr(k)=Tal(k)+Tal2(k)+Tal3(k);
end
al=alsr/3;
end

function [no] = poisIndex(nadlvo,kodlvo,dv)
n=length(dv); fn=1; fk=1;
for k=1:n
    if ((dv(k)>kodlvo) && (fk>0))
        fk=0;
        nok=k;
    end
    if ((dv(k)>nadlvo) && (fn>0))
        fn=0;
        non=k;
    end
end
no=[non nok];
end

function [kmdv] = vydPodmas(no,x)
non=no(1);
nok=no(2);
q=1;
for k=non:nok
    xk(q)=x(k);
    q=q+1;
end
kmdv=xk;
end