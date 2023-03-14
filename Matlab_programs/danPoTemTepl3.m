function [ isdan ] = danPoTemTepl3(temvs,tepv)
%Засыпка плоско-параллельная, повторные измерения
tepv1=0; tepv2=0; tem1=0; tem2=0;
tem1=temvs(19)+temvs(21)+temvs(23);
tem2=temvs(20)+temvs(22)+temvs(24);
tepv1=(tepv(19)*temvs(19)+tepv(21)*temvs(21)+tepv(23)*temvs(23))/tem1;
tepv2=(tepv(20)*temvs(20)+tepv(22)*temvs(22)+tepv(24)*temvs(24))/tem2;
tem1=tem1/3; tem2=tem2/3;
k1=(tepv2-tepv1)/(tem2-tem1); k2=tepv2-k1*tem2; 
isdan=[k1 k2];
end