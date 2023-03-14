function [ isdan ] = danPoTemTepl4(temvs,tepv)
%Засыпка вертикальная, после обжига при 1000 °С
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(1)+temvs(3);
tem2=temvs(2)+temvs(4);
tepv1=(tepv(1)*temvs(1)+tepv(3)*temvs(3))/tem1;
tepv2=(tepv(2)*temvs(2)+tepv(4)*temvs(4))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end