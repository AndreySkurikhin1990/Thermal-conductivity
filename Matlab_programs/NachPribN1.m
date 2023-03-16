function [ tem ] = NachPribN1(T1,T2,T3,x)
N=length(x);
h=x(N);
l=x(N)-x(N-1);
d=15e-3;
a11=d*d;
a12=d;
a21=h*h;
a22=h;
c=T1;
b1=T3-c;
b2=T2-c;
de=a11*a22-a21*a12;
de1=b1*a22-b2*a12;
de2=a11*b2-a21*b1;
a=de1/de;
b=de2/de;
te=0;
for k=1:N
    te(k)=a*(x(k)^2)+b*x(k)+c;
end
tem=te;
end