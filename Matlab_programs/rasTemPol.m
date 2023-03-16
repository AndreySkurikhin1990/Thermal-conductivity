function [ tem ] = rasTemPol(T0,a,h,mk,x)
n=mk(1);
q=mk(2);
ti=(T0-q)/n;
m=(mk(3)-n)/h;
tm1=(m*h^3/6+n*h^2/2)/a+q;
p=(mk(4)-tm1)/h;
lx=length(x);
for k=1:lx
    teta(k)=(m*x(k)+n)*ti+(m*x(k)^3/6+n*x(k)^2/2)/a+p*x(k)+q;
end
tem = teta;
end