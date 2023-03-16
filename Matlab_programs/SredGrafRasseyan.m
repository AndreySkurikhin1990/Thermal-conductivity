function [ alsr ] = SredGrafRasseyan()
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
format long g;
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; Tal8=0; Tal9=0;
dlv=dlvoVer53101(); Tal=privedkEdiPropus(TrkoVer5311()); Tal2=privedkEdiPropus(TrkoVer5312());
Tal3=privedkEdiPropus(TrkoVer5313()); Tal4=privedkEdiPropus(TrkoVer53101()); 
Tal5=privedkEdiPropus(TrkoVer53102()); Tal6=privedkEdiPropus(TrkoVer53103()); 
mkbr=[250.239 249.740 249.223 250.395 250.366 249.55]; 
mv=[0.283 0.464 0.812 0.22 0.547 0.777]; 
tol=[0.73 0.72 0.72 0.72 0.71 0.7];
rokbr=2.75; rov=0.49; n=length(mv);
xv=0; vv=0; vkbr=0;
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); %в см3
vv(k)=mv(k)/(1e3*rov); %в см3
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; %в мкм
end
alsr1=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv);
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; Tal8=0; Tal9=0;
Tal=privedkEdiPropus(TrkoVer5421()); Tal2=privedkEdiPropus(TrkoVer5422()); 
Tal3=privedkEdiPropus(TrkoVer5423()); Tal4=privedkEdiPropus(TrkoVer5431());
Tal5=privedkEdiPropus(TrkoVer5432()); Tal6=privedkEdiPropus(TrkoVer5433());
Tal7=privedkEdiPropus(TrkoVer5441()); Tal8=privedkEdiPropus(TrkoVer5442());
Tal9=privedkEdiPropus(TrkoVer5443());
mkbr=[250 249.629 249.294 249.706 249.510 249.307 250.328 249.604 249.206]; 
mv=[0.255 0.539 0.809 0.295 0.517 0.756 0.36 0.534 0.843]; 
tol=[0.72 0.71 0.7 0.7 0.73 0.72 0.74 0.7 0.76]; 
rokbr=2.75; rov=0.52; n=length(mv);
xv=0; vv=0; vkbr=0;
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
end
alsr2=opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv);
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; Tal8=0; Tal9=0;
Tal=privedkEdiPropus(TrkoVer5551()); Tal2=privedkEdiPropus(TrkoVer5552()); 
Tal3=privedkEdiPropus(TrkoVer5553()); Tal4=privedkEdiPropus(TrkoVer5561());
Tal5=privedkEdiPropus(TrkoVer5562()); Tal6=privedkEdiPropus(TrkoVer5563());
Tal7=privedkEdiPropus(TrkoVer5571()); Tal8=privedkEdiPropus(TrkoVer5572());
Tal9=privedkEdiPropus(TrkoVer5573());
mkbr=[249.913 249.607 249.218 249.929 249.695 249.306 250.405 249.625 249.348];
mv=[0.315 0.473 0.709 0.293 0.528 0.83 0.27 0.493 0.764]; 
tol=[0.74 0.74 0.72 0.72 0.71 0.7 0.78 0.73 0.76];
rokbr=2.75; rov=0.53; n=length(mv);
xv=0; vv=0; vkbr=0;
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
end
alsr3=opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv);
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=0; Tal2=0; Tal3=0; Tal4=0; Tal5=0; Tal6=0; Tal7=0; Tal8=0; Tal9=0;
Tal=privedkEdiPropus(TrkoVer5681()); Tal2=privedkEdiPropus(TrkoVer5682());
Tal3=privedkEdiPropus(TrkoVer5683()); Tal4=privedkEdiPropus(TrkoVer5691());
Tal5=privedkEdiPropus(TrkoVer5692()); Tal6=privedkEdiPropus(TrkoVer5693());
mkbr=[250.882 249.590 249.213 250.299 249.441 249.365];
mv=[0.320 0.533 0.849 0.223 0.502 0.797]; 
tol=[0.76 0.72 0.69 0.73 0.73 0.73]; 
rokbr=2.75; rov=0.56; n=length(mv); 
xv=0; vv=0; vkbr=0;
for k=1:n
vkbr(k)=mkbr(k)/(1e3*rokbr); 
vv(k)=mv(k)/(1e3*rov); 
xv(k)=(vv(k)/(vv(k)+vkbr(k)))*tol(k)*1e3; 
end
alsr4=opredAlphSred6(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,xv);
alsre=0;
nTal=length(alsr4);
for k=1:nTal
    alsre(k)=(alsr1(k)+alsr2(k)+alsr3(k)+alsr4(k));
end
alsre=alsre/4;
k=ZapisFile(1e6*alsre);
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
    j=1; Tal(k)=-log(Tal(k))/xv(j); j=j+1; 
    Tal2(k)=-log(Tal2(k))/xv(j); j=j+1;
    Tal3(k)=-log(Tal3(k))/xv(j); j=j+1;
    Tal4(k)=-log(Tal4(k))/xv(j); j=j+1;
    Tal5(k)=-log(Tal5(k))/xv(j); j=j+1;
    Tal6(k)=-log(Tal6(k))/xv(j); j=j+1;
end
alsr=0; s1=0;
for k=1:nTal 
    s1=Tal(k)+Tal2(k)+Tal3(k);
    s2=Tal4(k)+Tal5(k)+Tal6(k);
    alsr(k)=s1+s2; s1=0; s2=0;
end
al=alsr/6;
end

function [ al ] = opredAlphSred9(Tal,Tal2,Tal3,Tal4,Tal5,Tal6,Tal7,Tal8,Tal9,xv)
nTal=length(Tal);
for k=1:nTal
   j=1; al(k)=-log(Tal(k))/xv(j); j=j+1;
    al2(k)=-log(Tal2(k))/xv(j); j=j+1;
    al3(k)=-log(Tal3(k))/xv(j); j=j+1;
    al4(k)=-log(Tal4(k))/xv(j); j=j+1;
    al5(k)=-log(Tal5(k))/xv(j); j=j+1;
    al6(k)=-log(Tal6(k))/xv(j); j=j+1;
    al7(k)=-log(Tal7(k))/xv(j); j=j+1;
    al8(k)=-log(Tal8(k))/xv(j); j=j+1;
    al9(k)=-log(Tal9(k))/xv(j);
end
alsr=0; s1=0; s2=0;
for k=1:nTal 
    s1=al(k)+al2(k)+al3(k)+al4(k)+al5(k);
    s2=al6(k)+al7(k)+al8(k)+al9(k);
    alsr(k)=s1+s2; s1=0; s2=0;
end
al=alsr/length(xv);
end

function t = ZapisFile(massi)
fid = fopen('Koefficient_pogloscheniya-1.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end

function uv = usrednen(T,usv,npp,dv)
PP=6.6260755e-34;
PB=1.380658e-23;
c0=299792458;
c1=2*pi*PP*c0^2;
c2=PP*c0/PB;
iz0=0; iz1=0; dl=0;
c1t=c1; c2t=c2;
for j=1:length(dv)
    c1t=c1t/(npp(j)^2);
    c2t=c2t/npp(j);
    la=dv(j)/npp(j);
    dl(j)=la;
iz0(j)=c1t/((la^5)*(exp(c2t/(la*T))-1));
iz1(j)=usv(j)*iz0(j); %усредняемая величина
c1t=c1;
c2t=c2;
end
chi=trapz(dl,iz1);
zna=trapz(dl,iz0);
uv=chi/zna; %усредненная величина
end