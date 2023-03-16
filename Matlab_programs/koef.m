function ak = koef(x, y, na)
for k=1:na
    r=x(k);
    for s=1:na
        d=na-s;
    ax(k,s) = r^d;
    end
end
ak = inv(ax)*y';
end