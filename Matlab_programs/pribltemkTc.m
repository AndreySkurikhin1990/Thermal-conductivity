function [ tem ] = pribltemkTc(T,Th,Tc,x)
tet=0;
temp=T;
p=length(T);
h=x(p);
k1=(Th-Tc)/(x(1)-h);
k2=Tc-h*k1;
for k=1:p
    tet(k)=k1*x(k)+k2;
end
teta=0;
Tho=T(1);
Tco=T(p);
k1=(Tho-Tco)/(x(1)-h);
k2=Tco-h*k1;
for k=1:p
    teta(k)=k1*x(k)+k2;
end
det1=0;
det2=0;
det3=0;
for k=2:p
det1=(Th-teta(k));
det2=(Th-tet(k));
det3=det1-det2;
temp(k)=T(k)+det3;
end
temp(1)=Th;
tem=temp;
end