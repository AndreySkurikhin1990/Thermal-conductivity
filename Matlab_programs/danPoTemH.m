function [ isdan ] = danPoTemH(temvs,temvh)
%Засыпка плоско-параллельная, исходный
tho1=0; tho2=0; tem1=0; tem2=0;
tem1=temvs(15)+temvs(17);
tem2=temvs(16)+temvs(18);
tho1=(temvh(15)*temvs(15)+temvh(17)*temvs(17))/tem1;
tho2=(temvh(16)*temvs(16)+temvh(18)*temvs(18))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(tho2-tho1)/(tem2-tem1); k2=tho2-k1*tem2; 
isdan=[k1 k2];
end