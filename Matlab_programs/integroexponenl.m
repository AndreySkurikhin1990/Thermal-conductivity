function integl = integroexponenl(n,x)
a=1-n;
if (a>0)
gaa=gamma(a);
gaaz=gammainc(x,a)*gaa;
gaaz=gaa-gaaz;
inte=gaaz/x^a;
inte=(inte+exp(-x)*x^a)/a;
integl=inte;
end
end