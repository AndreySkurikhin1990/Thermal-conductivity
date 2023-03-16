function tesr = sredtemp(temvh,temvc,pt)
ts=0;
for k=1:pt
    ts=ts+temvh(k)+temvc(k);
end
ts=ts/2/pt;
tesr = ts;
end