function [ isdan ] = danPoTemC3(temvs,temvc)
%������� ������-������������, �������
tc1=0; tc2=0; tem1=0; tem2=0;
tem1=temvs(19)+temvs(21)+temvs(23);
tem2=temvs(20)+temvs(22)+temvs(24);
tc1=(temvc(19)*temvs(19)+temvc(21)*temvs(21)+temvc(23)*temvs(23))/tem1;
tc2=(temvc(20)*temvs(20)+temvc(22)*temvs(22)+temvc(24)*temvs(24))/tem2;
tem1=tem1/3; tem2=tem2/3;
k1=(tc2-tc1)/(tem2-tem1); k2=tc2-k1*tem2; 
isdan=[k1 k2];
end