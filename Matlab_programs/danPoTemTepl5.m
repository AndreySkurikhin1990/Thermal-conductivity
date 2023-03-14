function [ isdan ] = danPoTemTepl5(temvs,tepv)
%Засыпка вертикальная, повторы
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(9)+temvs(13);
tem2=temvs(10)+temvs(14);
tepv1=(tepv(9)*temvs(9)+tepv(13)*temvs(13))/tem1;
tepv2=(tepv(10)*temvs(10)+tepv(14)*temvs(14))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end