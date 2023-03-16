function [ alsr ] = SredGrafKoAbsVer()
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=provnaOtrZna(AbsKoVer5311()); Tal2=provnaOtrZna(AbsKoVer5312()); 
Tal3=provnaOtrZna(AbsKoVer5313()); Tal4=provnaOtrZna(AbsKoVer53101()); 
Tal5=provnaOtrZna(AbsKoVer53102()); Tal6=provnaOtrZna(AbsKoVer53103()); 
for k=1:length(Tal)
    Tal(k)=-Tal(k); Tal(k)=10^Tal(k);
    Tal2(k)=-Tal2(k); Tal2(k)=10^Tal2(k);
    Tal3(k)=-Tal3(k); Tal3(k)=10^Tal3(k);
    Tal4(k)=-Tal4(k); Tal4(k)=10^Tal4(k);
    Tal5(k)=-Tal5(k); Tal5(k)=10^Tal5(k);
    Tal6(k)=-Tal6(k); Tal6(k)=10^Tal6(k);
end
mkbr=250.239; mv=0.283; tol=0.73; rokbr=2.75; rov=0.49; 
mkbr2=249.740; mv2=0.464; tol2=0.72; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); 
vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
mkbr3=249.223; mv3=0.812; tol3=0.72; 
mkbr4=250.395; mv4=0.22; tol4=0.72; rov3=0.49; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); 
vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
mkbr5=250.336; mv5=0.547; tol5=0.71; rokbr=2.75; rov5=0.49; 
mkbr6=249.55; mv6=0.777; tol6=0.7; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); 
vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
xv=(vv/(vv+vkbr))*tol*1e3; 
xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3;
xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; 
xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
nTal=length(Tal);
for k=1:nTal 
    Tal(k)=-log(Tal(k))/xv; 
    Tal2(k)=-log(Tal2(k))/xv2;
    Tal3(k)=-log(Tal3(k))/xv3;
    Tal4(k)=-log(Tal4(k))/xv4;
    Tal5(k)=-log(Tal5(k))/xv5;
    Tal6(k)=-log(Tal6(k))/xv6; 
end
alsr1=0;
for k=1:nTal 
    alsr1(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k)+Tal6(k);
    alsr1(k)=alsr1(k)/6;
end
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=provnaOtrZna(AbsKoVer5421()); Tal2=provnaOtrZna(AbsKoVer5422());
mkbr=250; mv=0.255; tol=0.72; rokbr=2.75; rov=0.52; 
mkbr2=249.629; mv2=0.539; tol2=0.71; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
Tal3=provnaOtrZna(AbsKoVer5423()); Tal4=provnaOtrZna(AbsKoVer5431()); 
mkbr3=249.294; mv3=0.809; tol3=0.7;  
mkbr4=249.706; mv4=0.295; tol4=0.7; rov3=0.52; rov4=rov3;
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
Tal5=provnaOtrZna(AbsKoVer5432()); Tal6=provnaOtrZna(AbsKoVer5433()); 
mkbr5=249.510; mv5=0.517; tol5=0.73; rokbr=2.75; rov5=0.52; 
mkbr6=249.307; mv6=0.756; tol6=0.72; rov6=rov5;
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
Tal7=provnaOtrZna(AbsKoVer5441()); Tal8=provnaOtrZna(AbsKoVer5442());
Tal9=provnaOtrZna(AbsKoVer5443());
for k=1:length(Tal)
    Tal(k)=-Tal(k); Tal(k)=10^Tal(k);
    Tal2(k)=-Tal2(k); Tal2(k)=10^Tal2(k);
    Tal3(k)=-Tal3(k); Tal3(k)=10^Tal3(k);
    Tal4(k)=-Tal4(k); Tal4(k)=10^Tal4(k);
    Tal5(k)=-Tal5(k); Tal5(k)=10^Tal5(k);
    Tal6(k)=-Tal6(k); Tal6(k)=10^Tal6(k);
    Tal7(k)=-Tal7(k); Tal7(k)=10^Tal7(k);
    Tal8(k)=-Tal8(k); Tal8(k)=10^Tal8(k);
    Tal9(k)=-Tal9(k); Tal9(k)=10^Tal9(k);
end
mkbr7=250.328; mv7=0.36; tol7=0.74; rov7=0.52;
mkbr8=249.604; mv8=0.534; tol8=0.7; rov8=rov7;
vkbr7=mkbr7/(1e3*rokbr); vv7=mv7/(1e3*rov7); 
vkbr8=mkbr8/(1e3*rokbr); vv8=mv8/(1e3*rov8);
mkbr9=249.206; mv9=0.843; tol9=0.76; rov9=0.52;
vkbr9=mkbr9/(1e3*rokbr); vv9=mv9/(1e3*rov9); 
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
xv7=(vv7/(vv7+vkbr7))*tol7*1e3; xv8=(vv8/(vv8+vkbr8))*tol8*1e3;
xv9=(vv9/(vv9+vkbr9))*tol9*1e3;
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
end
alsr2=0;
for k=1:nTal 
    alsr2(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k);
    alsr2(k)=alsr2(k)+Tal6(k)+Tal7(k)+Tal8(k)+Tal9(k);
    alsr2(k)=alsr2(k)/9;
end
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=provnaOtrZna(AbsKoVer5551()); Tal2=provnaOtrZna(AbsKoVer5552()); 
mkbr=249.913; mv=0.315; tol=0.74; rokbr=2.75; rov=0.53; 
mkbr2=249.607; mv2=0.473; tol2=0.74; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
Tal3=provnaOtrZna(AbsKoVer5553()); Tal4=provnaOtrZna(AbsKoVer5561());
mkbr3=249.218; mv3=0.709; tol3=0.72; 
mkbr4=249.929; mv4=0.293; tol4=0.72; rov3=0.53; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
Tal5=provnaOtrZna(AbsKoVer5562()); Tal6=provnaOtrZna(AbsKoVer5563()); 
mkbr5=249.695; mv5=0.528; tol5=0.71; rokbr=2.75; rov5=0.53; 
mkbr6=249.306; mv6=0.83; tol6=0.7; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); 
vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
Tal7=provnaOtrZna(AbsKoVer5571()); Tal8=provnaOtrZna(AbsKoVer5572()); 
mkbr7=250.405; mv7=0.27; tol7=0.78; rov7=0.53;
mkbr8=249.625; mv8=0.493; tol8=0.73; rov8=rov7;
vkbr7=mkbr7/(1e3*rokbr); vv7=mv7/(1e3*rov7); 
vkbr8=mkbr8/(1e3*rokbr); vv8=mv8/(1e3*rov8);
Tal9=provnaOtrZna(AbsKoVer5573());  
for k=1:length(Tal)
    Tal(k)=-Tal(k); Tal(k)=10^Tal(k);
    Tal2(k)=-Tal2(k); Tal2(k)=10^Tal2(k);
    Tal3(k)=-Tal3(k); Tal3(k)=10^Tal3(k);
    Tal4(k)=-Tal4(k); Tal4(k)=10^Tal4(k);
    Tal5(k)=-Tal5(k); Tal5(k)=10^Tal5(k);
    Tal6(k)=-Tal6(k); Tal6(k)=10^Tal6(k);
    Tal7(k)=-Tal7(k); Tal7(k)=10^Tal7(k);
    Tal8(k)=-Tal8(k); Tal8(k)=10^Tal8(k);
    Tal9(k)=-Tal9(k); Tal9(k)=10^Tal9(k);
end
mkbr9=249.348; mv9=0.764; tol9=0.76; rov9=0.53;
vkbr9=mkbr9/(1e3*rokbr); vv9=mv9/(1e3*rov9); 
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
xv7=(vv7/(vv7+vkbr7))*tol7*1e3; xv8=(vv8/(vv8+vkbr8))*tol8*1e3;
xv9=(vv9/(vv9+vkbr9))*tol9*1e3;
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
end
alsr3=0;
for k=1:nTal 
    alsr3(k)=Tal(k)+Tal2(k)+Tal3(k)+Tal4(k)+Tal5(k);
    alsr3(k)=alsr3(k)+Tal6(k)+Tal7(k)+Tal8(k)+Tal9(k);
    alsr3(k)=alsr3(k)/9;
end
%rom60_1853_1_10 = 0.49; ro60_100_1854_234=0.52; ro100_150_1855_567=0.53; ro150_200_1856_89 = 0.56; 
Tal=provnaOtrZna(AbsKoVer5681()); Tal2=provnaOtrZna(AbsKoVer5682()); 
mkbr=250.882; mv=0.320; tol=0.76; rokbr=2.75; rov=0.56; 
mkbr2=249.590; mv2=0.533; tol2=0.72; rov2=rov; 
vkbr=mkbr/(1e3*rokbr); vv=mv/(1e3*rov); vkbr2=mkbr2/(1e3*rokbr); vv2=mv2/(1e3*rov2);
Tal3=provnaOtrZna(AbsKoVer5683()); Tal4=provnaOtrZna(AbsKoVer5691()); 
mkbr3=249.213; mv3=0.849; tol3=0.69; 
mkbr4=250.299; mv4=0.223; tol4=0.73; rov3=0.56; rov4=rov3; 
vkbr3=mkbr3/(1e3*rokbr); vv3=mv3/(1e3*rov3); vkbr4=mkbr4/(1e3*rokbr); vv4=mv4/(1e3*rov4); t0=273.15;
Tal5=provnaOtrZna(AbsKoVer5692()); Tal6=provnaOtrZna(AbsKoVer5693()); 
for k=1:length(Tal)
    Tal(k)=-Tal(k); Tal(k)=10^Tal(k);
    Tal2(k)=-Tal2(k); Tal2(k)=10^Tal2(k);
    Tal3(k)=-Tal3(k); Tal3(k)=10^Tal3(k);
    Tal4(k)=-Tal4(k); Tal4(k)=10^Tal4(k);
    Tal5(k)=-Tal5(k); Tal5(k)=10^Tal5(k);
    Tal6(k)=-Tal6(k); Tal6(k)=10^Tal6(k);
end
mkbr5=249.441; mv5=0.502; tol5=0.73; rokbr=2.75; rov5=0.56; 
mkbr6=249.365; mv6=0.797; tol6=0.73; rov6=rov5; 
vkbr5=mkbr5/(1e3*rokbr); vv5=mv5/(1e3*rov5); vkbr6=mkbr6/(1e3*rokbr); vv6=mv6/(1e3*rov6);
xv=(vv/(vv+vkbr))*tol*1e3; xv2=(vv2/(vv2+vkbr2))*tol2*1e3; 
xv3=(vv3/(vv3+vkbr3))*tol3*1e3; xv4=(vv4/(vv4+vkbr4))*tol4*1e3;
xv5=(vv5/(vv5+vkbr5))*tol5*1e3; xv6=(vv6/(vv6+vkbr6))*tol6*1e3; 
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
    alsr4(k)=alsr4(k)/6;
end
alsre=0;
for k=1:nTal
    alsre(k)=(alsr1(k)+alsr2(k)+alsr3(k)+alsr4(k))/4;
end
k=ZapisFile(1e6*alsre);
alsr=alsre;
end

function [ obnul ] = provnaOtrZna(ar)
arr=ar;
p=length(arr);
for k=1:p
    if (arr(k)<0)
        arr(k)=0;
    end
end
obnul=arr;
end

function t = ZapisFile(massi)
fid = fopen('Koefficient_pogloscheniya_ver-2.txt','w');
p=length(massi);
for k=1:p
    fprintf(fid,'%0.20f\n',massi(k));
end
fclose(fid);
t=0;
end