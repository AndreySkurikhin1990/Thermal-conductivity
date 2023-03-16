a = [50 25  20 10 5];
b=a;
for j=1:length(a)
    a(j)=a(j)/5;
end
b=b./5;
disp(a);
disp(b);