function [ isdan ] = danPoTemTepl6(temvs,tepv)
%Засыпка плоско-параллельная, после обжига при 1000 °С
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(25);
tem2=temvs(26);
tepv1=(tepv(25)*temvs(25))/tem1;
tepv2=(tepv(26)*temvs(26))/tem2;
tem1=tem1; tem2=tem2;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end