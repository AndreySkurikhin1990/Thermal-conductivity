%определяет номер билета и порядок выступления в случайном порядке
function t = tmp45()
Ni=1e5;
ce=4e1;
for k=1:ce
    a(k)=0;
    b(k)=k;
end
for k=1:Ni
    j=fix(ce*rand());
    a(j+1)=a(j+1)+1;
end
a=a/Ni;
ce=7;
q=1;
for k=1:ce
    c(k)=0;
end
f=1;
while (f>0)
for k=1:ce %q - номер выступления
    j=fix(ce*rand())+1; %номер студента
        if (c(j)>0)
        else
        c(j)=q;
        q=q+1;
        end
end
for k=1:ce-1
    for j=k+1:ce
        if (c(k)==c(j))
            c(j)=0;
        end
    end
end
f=0;
for k=1:ce
    if (c(k)==0) 
        f=1;
    end
end
end
c=c';
cs=8;
cb=4e1;
c=0;
for k=1:cs
    c(k)=0;
end
    f=1;
while (f>0)
for k=1:cs %номер студента
    j=fix(cb*rand())+1; %номер билета
        if (c(k)>0)
        else
        c(k)=j;
        end
end
f=0;
for k=1:cs-1
    for j=k+1:cs
        if (c(k)==c(j))
            c(j)=0;
        end
    end
end
f=0;
for k=1:cs
    if (c(k)==0) 
        f=1;
    end
end
end
c=c'
pl=plot(b,a,'-b');
set(pl,'LineWidth',3); hold on; grid on; 
xlabel({'Номер билета'}); 
ylabel({'Вероятность попадания'}); 
title({'График № 1'});
t=0;
end