function vf = porEmisEval()
vfs=squarPor();
vfr=rectPor();
vf=cylPor();
end

function kv = squarPor()
format long g;
F = 0;
rapo=0;
a=5e-1;
b=5e-1;
c=5e-1;
rapo(1) = a; rapo(2) = b; rapo(3) = c;
F(1) = view(38, 3, rapo);
rp=0;
w=a; l=b; h=c;
rp(1) = h; rp(2) = l; rp(3) = w;
F(2) = view(39, 3, rp);
F(3) = F(2);
l=a; w=b;
rpp=0;
rpp(1) = h; rpp(2) = l; rpp(3) = w;
F(4) = view(39, 3, rpp);
F(5) = F(4);
s=0;
for j=1:length(F)
    s=s+F(j);
end
s_sq=1*s;
F_sq=1*F';
kv=0;
end

function cy = cylPor()
A=3e-1;
R1=1e-1;
R2=R1;
rapo=0; F=0; rp=0;
rapo(1) = A; 
rapo(2) = R1; 
rapo(3) = R2;
F(1) = view(40, 3, rapo);
rp(1)=A;
rp(2)=R2;
F(2)=view(43, 2, rp);
s_cyl=F(1)+F(2);
F_cyl=1*F';
F(3)=view(42, 2, rp);
F(3)=(1-F(3))/2;
pr_cyl=pi*(R1^2)*F(2)-F(3)*2*pi*R1*A;
cy=0;
end

function rc = rectPor()
a=2e-2;
b=3e-2;
c=4e-2;
for k=1:6
    for j=1:6
        F(j,k)=0;
    end
end
rapo=0;
rapo(1) = a; rapo(2) = b; rapo(3) = c;
F(1,2) = view(38, 3, rapo);
w=a; l=b; h=c;
rapo(1) = h; rapo(2) = l; rapo(3) = w;
F(1,3) = view(39, 3, rapo);
F(1,4) = F(1,3);
l=a; w=b;
rapo(1) = h; rapo(2) = l; rapo(3) = w;
F(1,5) = view(39, 3, rapo);
F(1,6) = F(1,5);
w=a; h=c; l=b;
rapo(1) = h; rapo(2) = l; rapo(3) = w;
F(2,3) = view(39, 3, rapo);
F(2,4)=F(2,3);
w=a; h=c; l=b; l=a; w=b;
rapo(1) = h; rapo(2) = l; rapo(3) = w;
F(2,5) = view(39, 3, rapo);
F(2,6)=F(2,5);
rapo(1) = c; rapo(2) = b; rapo(3) = a;
F(3,4) = view(38, 3, rapo);
l=c; h=a; w=b;
rapo(1) = h; rapo(2) = l; rapo(3) = w;
F(3,5) = view(39, 3, rapo);
F(3,6)=F(3,5);
l=c; w=b; h=a;
rapo(1) = h; rapo(2) = l; rapo(3) = w;
F(4,5) = view(39, 3, rapo);
F(4,6) = F(4,5);
rapo(1) = a; rapo(2) = c; rapo(3) = b;
F(5,6) = view(38, 3, rapo);
A=0;
A(1)=a*b;
A(2)=A(1);
A(3)=b*c;
A(4)=A(3);
A(5)=a*c;
A(6)=A(5);
T(1)=5e2;
for k=2:5
    T(k)=52e1;
end
T(6)=54e1;
npp=Kramers_n();
for k=1:length(npp)
    ep(k)=epsilnu(npp(k));
end
sig=5.67e-8;
for k=1:6
    eps(k)=epssred(RasshDiapDlinVoln(),ep,T(k),npp);
    ho(k)=0;
    pin(k)=sig*eps(k)*T(k)^4;
    id(k)=1;
end
pout = GRAYDIFF(1, A, eps, ho, F, id, pin);
pout = 1*pout'
rc=checkViewFactors(b,F);
end

function ch = checkViewFactors(b,F)
[m n]=size(F);
for k=1:n
    s=0;
for j=1:m
    s=s+F(k,j);
end
s_sum=1*s;
end
F_i=1*F';
rpp=0;
%--------
nb=[2e-1 3e-1 5e-1 1e0 2e0 3e0 4e0 5e0 1e1 2e1]; lb=[1e-1 2e-1 4e-1 6e-1 1e0 2e0 4e0 6e0 1e1 2e1];
for j=1:length(nb)
for k=1:length(lb)
a=nb(j)*b; c=lb(k)*b;
l=b; h=a; w=c;
rpp(1) = h; rpp(2) = l; rpp(3) = w;
vfh=view(39, 3, rpp);
%disp(lb(k));
end
%disp(nb(j));
end
xb=[1e-1 2e-1 4e-1 6e-1 1e0 2e0 4e0 1e1]; yb=[1e-1 2e-1 4e-1 6e-1 1e0 2e0 4e0 6e0 1e1 2e1];
c=1; 
for j=1:length(xb)
for k=1:length(yb)
b=xb(j)*c; a=yb(k)*c;
rpp(1) = a; rpp(2) = b; rpp(3) = c;
vfh = view(38, 3, rpp);
%disp(yb(k));
end
%disp(xb(j));
end
ch=0;
end