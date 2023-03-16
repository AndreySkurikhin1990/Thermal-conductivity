function [ kotepar ] = koeftepv(qv,y0,temvh,temvc,pt)
tepv=0;
for k=1:pt
tepv(k)=qv(k)*y0/(temvh(k)-temvc(k));
end
kotepar=tepv;
end