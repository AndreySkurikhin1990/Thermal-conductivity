function tm = tmp23()
format long g;
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=privedkEdiPropus(TrkoVer5421()); Tal2=privedkEdiPropus(TrkoVer5422()); 
mkbr=250; mv=0.255; tol=0.72; rokbr=2.75; rov=0.52; 
mkbr2=249.629; mv2=0.539; tol2=0.71; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
Tal3=privedkEdiPropus(TrkoVer5423()); Tal4=privedkEdiPropus(TrkoVer5431());
mkbr3=249.294; mv3=0.809; tol3=0.7;  
mkbr4=249.706; mv4=0.295; tol4=0.7; rov3=0.52; rov4=rov3;
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
Tal5=privedkEdiPropus(TrkoVer5432()); Tal6=privedkEdiPropus(TrkoVer5433());
mkbr5=249.510; mv5=0.517; tol5=0.73; rokbr=2.75; rov5=0.52; 
mkbr6=249.307; mv6=0.756; tol6=0.72; rov6=rov5;
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
Tal7=privedkEdiPropus(TrkoVer5441()); Tal8=privedkEdiPropus(TrkoVer5442());
Tal9=privedkEdiPropus(TrkoVer5443());
mkbr7=250.328; mv7=0.36; tol7=0.74; rov7=0.52;
mkbr8=249.604; mv8=0.534; tol8=0.7; rov8=rov7;
vkbr7=mkbr7/(1e3*rokbr); vv7=mv7/(1e3*rov7); 
vkbr8=mkbr8/(1e3*rokbr); vv8=mv4/(1e3*rov8);
mkbr9=249.206; mv9=0.843; tol9=0.76; rov9=0.52;
vkbr9=mkbr9/(1e3*rokbr); vv9=mv9/(1e3*rov9); 
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
xv7=(vv7/(vv7+vkbr7))*tol7*1e3; xv8=(vv8/(vv8+vkbr8))*tol8*1e3;
xv9=(vv9/(vv9+vkbr9))*tol9*1e3;
nTal=length(Tal);
for k=1:nTal
    Tal(k)=-log(Tal(k))/xv;
    Tal2(k)=-log(Tal2(k))/xv2;
    Tal3(k)=-log(Tal3(k))/xv3;
    Tal4(k)=-log(Tal4(k))/xv4;
    Tal5(k)=-log(Tal5(k))/xv5;
    Tal6(k)=-log(Tal6(k))/xv6;
    Tal7(k)=-log(Tal7(k))/xv7;
    Tal8(k)=-log(Tal8(k))/xv8;
    Tal9(k)=-log(Tal9(k))/xv9; 
end;
alsr2=0;
for k=1:nTal 
    alsr2(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k);
    alsr2(k)=alsr2(k)+Tal6(k)+Tal7(k)+Tal8(k)+Tal9(k);
end
alsr2=alsr2/9;
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
        kopo(j)=alsr2(k);
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