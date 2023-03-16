function tm = tmp25()
format long g;
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=privedkEdiPropus(TrkoVer5681()); Tal2=privedkEdiPropus(TrkoVer5682());
mkbr=250.882; mv=0.320; tol=0.76; rokbr=2.75; rov=0.56; 
mkbr2=249.590; mv2=0.533; tol2=0.72; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
Tal3=privedkEdiPropus(TrkoVer5683()); Tal4=privedkEdiPropus(TrkoVer5691());
mkbr3=249.213; mv3=0.849; tol3=0.69; 
mkbr4=250.299; mv4=0.223; tol4=0.73; rov3=0.56; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
Tal5=privedkEdiPropus(TrkoVer5692()); Tal6=privedkEdiPropus(TrkoVer5693());
mkbr5=249.441; mv5=0.502; tol5=0.73; rokbr=2.75; rov5=0.56; 
mkbr6=249.365; mv6=0.797; tol6=0.73; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
nTal=length(Tal);
for k=1:nTal
    Tal(k)=-log(Tal(k))/xv;
    Tal2(k)=-log(Tal2(k))/xv2;
    Tal3(k)=-log(Tal3(k))/xv3;
    Tal4(k)=-log(Tal4(k))/xv4;
    Tal5(k)=-log(Tal5(k))/xv5;
    Tal6(k)=-log(Tal6(k))/xv6; 
end
alsr4=0;
for k=1:nTal 
    alsr4(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k)+Tal6(k);
end
alsr4=alsr4/6;
dl = 15e-1:5e-1:15e0;
nd=length(dl);
f=1;
wl = dlvoVer53101();
j=1;
kopo=0;
n = length(wl);
for k=1:n
    wl(k)=1e4/wl(k);
end
for j=1:nd
    for k=1:n
    if (wl(k)>dl(j))
        if (f>0)
        kopo(j)=alsr4(k);
        f=0;
        end
    end
    end
    f=1;
end
kopo=1*kopo'
tm=0;
end

function [ prked ] = privedkEdiPropus(ar)
arr=ar;
p=length(arr);
for k=1:p
    if (arr(k)>1)
        arr(k)=1;
    end
end
prked=arr;
end