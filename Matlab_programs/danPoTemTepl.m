function [ isdan ] = danPoTemTepl(temvs,tepv)
%Засыпка плоско-параллельная, исходный
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(15)+temvs(17);
tem2=temvs(16)+temvs(18);
tepv1=(tepv(15)*temvs(15)+tepv(17)*temvs(17))/tem1;
tepv2=(tepv(16)*temvs(16)+tepv(18)*temvs(18))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end