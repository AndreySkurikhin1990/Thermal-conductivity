function [ isdan ] = danPoTemTeplN1(tepv1,tepv2,tepv3,T1,T2,T3)
%Засыпка хаотичная, фракция 2-0,7, исходный
a11=T1^2;
a12=T1;
a13=1;
a21=T2^2;
a22=T2;
a23=1;
a31=T3^2;
a32=T3;
a33=1;
b1=tepv1;
b2=tepv2;
b3=tepv3;
de=a11*(a22*a33-a32*a23)-a12*(a21*a33-a31*a23)+a13*(a21*a32-a31*a22);
de1=b1*(a22*a33-a32*a23)-a12*(b2*a33-b3*a23)+a13*(b2*a32-b3*a22);
de2=a11*(b2*a33-b3*a23)-b1*(a21*a33-a31*a23)+a13*(a21*b3-a31*b2);
de3=a11*(a22*b3-a32*b2)-a12*(a21*b3-a31*b2)+b1*(a21*a32-a31*a22);
a=de1/de;
b=de2/de;
c=de3/de;
isdan=[a b c];
end