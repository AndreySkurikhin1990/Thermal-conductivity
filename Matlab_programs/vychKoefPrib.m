function [ m ] = vychKoefPrib(timeAr, temHotAr, temColdAr)
p=length(timeAr);
for k=1:p
    x2(k)=timeAr(k)^2;
    xyh(k)=temHotAr(k)*timeAr(k);
    xyc(k)=temColdAr(k)*timeAr(k);
    x(k)=timeAr(k);
    y1(k)=temHotAr(k);
    y2(k)=temColdAr(k);
end
xyhs=0;
xycs=0;
x2s=0;
xs=0;
ysh=0;
ysc=0;
for k=1:p
    xyhs=xyhs+xyh(k);
    xycs=xycs+xyc(k);
    x2s=x2s+x2(k);
    xs=xs+x(k);
    ysh=ysh+y1(k);
    ysc=ysc+y2(k);
end
ah=(p*xyhs-xs*ysh)/(p*x2s-xs*xs);
bh=(ysh-ah*xs)/p;
ac=(p*xycs-xs*ysc)/(p*x2s-xs*xs);
bc=(ysc-ac*xs)/p;
m = [ah bh ac bc];
end