function [ r ] = TmpReshZad()
format long g;
L1=1e-3;
L2=250e-6;
C1=1e3*1e-12;
C2=4e3*1e-12;
r1=10;
r2=5;
M=50e-6;
om=0;
C1=C1-50*1e-12;
om01=1/sqrt(L1*C1)
om02=1/sqrt(L2*C2)
C2=RootsFind(C2,L1,C1,om01,M,L2,r1,r2)
%C2=L1*C1/L2;
%C2=((om01^2)*L1*L2-r1*r2-L2/C1+(om01*M)^2);
%C22=(L1-1/om01^2/C1);
%C2=(C2/C22)^(-1)
sh=2e5;
na=om01-sh;
ko=om01+sh;
om=0;
om=na:1e0:ko;
n=length(om);
y1=0; 
%y2=0; y3=0;
for k=1:n
y1(k)=RascY(om(k),om01,L1,L2,M,C1,C2,r1,r2);
%y2(k)=RascY(om(k),om01,L1,L2,M,C1,C2+620*1e-12,r1,r2);
%y3(k)=RascY(om(k),om01,L1,L2,M,C1,C2-620*1e-12,r1,r2);
end
%pl=plot(om,y1,'-k',om,y2,'-b',om,y3,'-r');
pl=plot(om,y1,'-k');
set(pl,'LineWidth',2); 
hold on; grid on;
xlabel({'Циклическая частота, рад/с'}); 
ylabel({'Относительное значение силы тока'}); 
%title({'График зависимости безразмерной силы тока от круговой частоты'});
r=0;
end

function y = RascY(om,om01,L1,L2,M,C1,C2,r1,r2)
ch=(om01*M)^2+r1*r2;
zn=(r1*r2-L1*L2*om^2+L2/C1+L1/C2-1/om^2/C1/C2+(om*M)^2);
zn=zn^2+(om*L1*r2-r2/om/C1+om*r1*L2-r1/om/C2)^2;
y=ch*om/om01/zn;
end

function kor = RootsFind(as,L1,C1,om,M,L2,r1,r2)
del=3.9e3*1e-12; 
aa=as-del;
ab=as+del;
nh=1e3;
eps=1; e=1e-7; h=0;
%r1=r1
%r2=r2
%L2=L2
%M=M
%om=om
%L1=L1
%C1=C1
%as=as
while (eps>e)
ac=(aa+ab)/2;    
%fc=(om*L1-1/om/C1)-((om*M)^2)*(om*L2-1/om/ac)/(r2^2+(om*L2-1/om/ac)^2);
%fa=(om*L1-1/om/C1)-((om*M)^2)*(om*L2-1/om/aa)/(r2^2+(om*L2-1/om/aa)^2);
%fb=(om*L1-1/om/C1)-((om*M)^2)*(om*L2-1/om/ab)/(r2^2+(om*L2-1/om/ab)^2);
%fc=(om*M)^2;
%fc=fc-om*L1*L2
%fc=fc+L2/C1
%fc=fc+L1/ac
%fc=fc-1/C1/ac/om^2
%fa=(om*M)^2
%fa=fa-om*L1*L2
%fa=fa+L2/C1
%fa=fa+L1/aa
%fa=fa-1/C1/aa/om^2
%fb=(om*M)^2
%fb=fb-om*L1*L2
%fb=fb+L2/C1
%fb=fb+L1/ab
%fb=fb-1/C1/ab/om^2
fc=om*L1*r2-r2/om/C1+om*r1*L2-r1/om/ac;
fa=om*L1*r2-r2/om/C1+om*r1*L2-r1/om/aa;
fb=om*L1*r2-r2/om/C1+om*r1*L2-r1/om/ab;
if ((fa*fc)<0)
    if ((fb*fc)>0)
    ab=ac; 
    end
end
if ((fa*fc)>0) 
    if ((fb*fc)<0) 
    aa=ac; 
    end
end
eps=abs(fa-fb);
h=h+1;
if (h>nh) 
    break; 
end
end
%h=(om*M)^2
%h=h-om*L1*L2
%h=h+L1/ac
%h=h-1/C1/ac/om^2
kor=ac;
end