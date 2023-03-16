%Определение КТП ТК шамота, ИТОМ, КВИ и вермикулита
function [ t ] = tmp55()
format longg; vybves=1; %vysh=2; 
te0=273.15; tem = 2e2:1e2:12e2; tem=tem+te0;
y0=3e1*1e-3; %vmivmf=1; vfv=1; vyuv=2; vysv=2;
nvyfv=[0, 1]; %фракции вермикулита
nvysv=[0, 1, 2]; %состояния вермикулита
nvmivmf=[0, 1, 2]; %стационарные методы измерений - 2019 и 2020
nvyuv=[1, 2]; %укладка вермикулита
b=length(tem); eps=1e-8;
for k=1:b
    s(k)=0;
    c(k)=0;
end
switch (vyves)
    case (1)
for kvfv=1:length(nvyfv)
for kvysv=1:length(nvysv)
for kvyuv=1:length(nvyuv)
for kvmivmf=1:length(nvmivmf)
vyfv=nvyfv(kvfv)
vysv=nvysv(kvysv)
vmivmf=nvmivmf(kvmivmf)
vyuv=nvyuv(kvyuv)
ktptkv=opredKTPTverKarkVerm(vyfv, tem, vysv, vyuv, vmivmf, y0)
for k=1:length(ktptkv)
    s(k)=s(k)+ktptkv(k);
    if (ktptkv(k)>eps)
    c(k)=c(k)+1;
    end
end
if (vyfv==1)
    break;
end
if (vysv>0)
    break;
end
%for vmkvi=4:10
%vmkvi=vmkvi'
%ktptkk=opredKTPTverKarkkvi(vmkvi,tem+te0);
%end
%t=ModelirTeploperItom(ktptks,ktptkv,tem);
%vmitom=0; t=RazPorItom(vmitom); 
end
if (vyfv==0)
    break;
end
end
end
end
    case (0)
        for vysh=0:10
            vysh=vysh
    ktptks=opredKTPTverKarkSha(tem, vysh)
    for k=1:length(ktptks)
    s(k)=s(k)+ktptks(k);
    if (ktptks(k)>eps)
    c(k)=c(k)+1;
    end
    end
        end
end
c=c;
for k=1:length(s)
    if (c(k)>eps)
        s(k)=s(k)/c(k);
    end
end
t=s;
end