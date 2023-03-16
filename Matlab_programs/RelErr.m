function reler = RelErr(l,lv,t1,t2,asp, kotepv, kotepvo, eps)
asr=asp/(asp+1);
%lamvo=polyval(kotepvo,temv);
for k=1:5 
    polj(k)=kotepv(k)/(6-k);
end
polj(6)=0;
jt1=polyval(polj,t1)/l;
jv=polyval(polj,t2)/l-jt1;
tepf=0.24;
dt=(t2-t1)/(0.5*l/lv);
%--------------------
mr=1;
for n=4:2:201
ddt=dt/(n-1);
temz=0;
for k=1:n 
    temz(k)=t2-(k-1)*ddt;
end
tol=0;
for k=1:(n-1) 
    if rem(k,2)==1
    tol(k)=(asr)*lv/(n-1);
    end
    if rem(k,2)==0
    tol(k)=(1-asr)*lv/(n-1);
    end
end
    a=0;
for k=1:2:(n-1)
    a(k,k)=tepf/tol(k);
    a(k,k+1)=-a(k,k);
end
for k=2:2:(n-1)
    a(k,k)=0.5*(polyval(kotepvo,temz(k))+polyval(kotepvo,temz(k+1)))/tol(k);
    a(k,k+1)=-a(k,k);
end
a(n,:)=0;
sig=5.67037e-8; stch=0.735;
no=1.616;
ro=((no-1)/(no+1))^2;
tau=1-ro;
ast=0;
ast(1,1)=-eps;
ast(n,n)=eps;
d=0;
for k=2:2:(n-1)
    d=n-k;
    ast(k,k)=eps*(tau^d)*((1-ro*ro)^(d/2));
    ast(k,k+1)=-ast(k,k);
end
for k=1:5
polbi(k)=kotepv(k)/(6-k);
end
polbi(6)=0;
b=0;
for k=1:n
    b(k)=polyval(polbi,temz(k))/l-jt1;
end
c=temz';
m=0;
while sum(abs(c))>1e-8
f=0;
for d=1:n
    for k=1:n
        f(d,k)=a(d,k)+4*ast(d,k)*temz(k)^3;
    end
end
fsm=0;
for d=1:n
    for k=1:n
        fsm(d,k)=a(d,k)*temz(k)+ast(d,k)*temz(k)^4;
    end
end
mfi=0;
mfi=-sum(fsm')+b;
c=inv(f)*mfi';
temz=c'+temz;
m=m+1;
end
js=b(1:(n-1))*tol'/lv;
lamef=js*lv/(temz(1)-temz(n));
lam=polyval(kotepv,t2-0.5*dt);
rer(mr)=abs(lamef-lam)*100/lam;
mr=mr+1;
end
reler=rer;
end