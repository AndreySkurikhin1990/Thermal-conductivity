function [ ns ]  = PokazPrelomPribl(x, xs, ys)
p=length(x);
for k=1:(p-1)
	ko=(ys(k+1)-ys(k))/(xs(k+1)-xs(k)); 
	y(k)=ko*(x(k)-xs(k))+ys(k);
end
k=p-1;
ko=(ys(k)-ys(k-1))/(xs(k)-xs(k-1));
y(k)=(x(k)-xs(k-1))*ko+ys(k);
ko=(ys(p)-ys(p-1))/(xs(p)-xs(p-1));
y(p)=(x(p)-xs(p-1))*ko+ys(p);
ns=y;
end