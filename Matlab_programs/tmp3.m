function t = tmp3()
format longg;
t=tm33();
end

function t = tm31()
%tolvt=3e4; tol=1e3; t0=273.15; 
%T0=563+t0; Tvt=141+t0;
%Nt=tolvt/tol;
%Te=0; Te=tempvtst(T0,Tvt,tol,Nt)-t0; 
%koo=0; koo=0:tol:tolvt;
%disp(length(koo));
%disp(length(Te));
%b=plot(koo,Te,'-b');
%set(b,'LineWidth',3);
%xlabel('Координата, мкм');
%ylabel('Температура, К');
%title('График T(x)');
t=0;
end

function t = tm32()
%a=[5,10,6,10,8; 5,4,7,7,6; 1,2,7,1,5; 6,10,3,6,2; 8,10,10,9,10];
%b=[6,6,4,8,9];
%x=inv(a)*b'
a=[1,0,-1/4]; b=[0,2]; [b1,a1]=eqtflength(b,a);
[z,p,G]=tf2zp(b1,a1);
z=z
p=p
[r,p,c]=residuez(b1,a1);
r=r
c=c
for k=1:length(r)
resi(k)=r(k)*p(k);
end
resi=resi'
t=0;
end

function t = tm33()
x=2e1; n=2;
intx20=integroexpon(n,x);
intx21=1; chel=50:50:150;
for q=1:length(chel)
    s=1; zn=1;
for k=1:chel(q)
    t=1; z=1;
    for j=1:k
        t=t*(n+j-1);
        z=z*x;
    end
    zn=-zn; t=t'; z=z';
    s=s+zn*t/z;
end
s=s'
end
t=0;
end