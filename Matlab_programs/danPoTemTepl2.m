function [ isdan ] = danPoTemTepl2(temvs,tepv)
%Засыпка вертикальная, исходный
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(5)+temvs(7)+temvs(11);
tem2=temvs(6)+temvs(8)+temvs(12);
tepv1=(tepv(5)*temvs(5)+tepv(7)*temvs(7)+tepv(11)*temvs(11))/tem1;
tepv2=(tepv(6)*temvs(6)+tepv(8)*temvs(8)+tepv(12)*temvs(12))/tem2;
tem1=tem1/3; tem2=tem2/3;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end