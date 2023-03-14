function [ isdan ] = danPoPlotTeplPot(temvs,qv)
%Засыпка плоско-параллельная, исходный
qv1=0; qv2=0; tem1=0; tem2=0;
tem1=temvs(15)+temvs(17);
tem2=temvs(16)+temvs(18);
qv1=(qv(15)*temvs(15)+qv(17)*temvs(17))/tem1;
qv2=(qv(16)*temvs(16)+qv(18)*temvs(18))/tem2;
tem1=tem1/2; tem2=tem2/2;
k1=(qv2-qv1)/(tem2-tem1); k2=qv2-k1*tem2; 
isdan=[k1 k2];
end