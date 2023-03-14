function [ isdan ] = danPoTemTepl2071(temvs,tepv)
%Засыпка исходная, фракция 2-0,7 мм
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(1);
tem2=temvs(2);
tepv1=(tepv(1)*temvs(1))/tem1;
tepv2=(tepv(2)*temvs(2))/tem2;
tem1=tem1/1; tem2=tem2/1;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end