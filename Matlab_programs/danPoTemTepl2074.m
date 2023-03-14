function [ isdan ] = danPoTemTepl2074(temvs,tepv)
%Фракция 2-0,7 мм, после повторного обжига при 1000 °С
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(3);
tem2=temvs(4);
tepv1=(tepv(3)*temvs(3))/tem1;
tepv2=(tepv(4)*temvs(4))/tem2;
tem1=tem1; tem2=tem2;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end