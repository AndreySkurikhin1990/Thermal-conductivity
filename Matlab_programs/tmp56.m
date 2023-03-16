function t = tmp56()
format longg;
%t=provravnrasp();
di=2*erfcinv(1e-4);
mini=0; maxi=6e1; 
mo=(maxi+mini)/2 %МО
si=abs((mini-maxi)/di)/sqrt(2); %СКО
si=si*sqrt(2)
t=provnormrasp(mo,si,mini,maxi);
end

function t = provravnrasp()
%rng('shuffle','v5uniform');
%rng('default');
na=1; dl=1; ko=na+dl;
Ns=1e7; No=1e3; h=dl/No;
x=na:h:ko; lx=length(x);
for k=1:lx
    y(k)=0;
end
for k=1:Ns
    xk=na+dl*rand();
    no=poisknomera(x,lx,xk);
    if (no>0)
    y(no)=y(no)+1;
    end
end
%y=y/Ns;
for k=1:lx-1
    xn(k)=(x(k)+x(k+1))/2;
    yn(k)=(y(k)+y(k+1))/2;
end
pl=plot(xn,yn,'.b');
xlabel('x'); ylabel('y');
%ylim([0 1]);
set(pl,'LineWidth',2); hold on; grid on;
title({'y(x)'});
legend('y(x)','location','best')
t=0;
end

function t = provnormrasp(mo,si,mini,maxi)
%rng('shuffle','v4');
%rng(0);
%rng('default');
Ns=1e6; h=1e-3; r=2e0;
x=-(r*si):h:(si*r); x=x+mo; lx=length(x);
for k=1:lx
    y(k)=0;
end
b=0;
for k=1:Ns
    px=mo+si*randn();
    if ((px>maxi) || (px<mini))
        b=b+1;
    end
    no=poisknomera(x,lx,px);
    if (no>0)
    y(no)=y(no)+1;
    end
end
b=b*1e2/Ns
plo=trapz(x,y);
y=y/plo;
b=0;
for k=1:lx
    if (abs(x(k)-mo)>=si)
    b=b+y(k);
    end
end
b=b/Ns;
%y=y/Ns;
for k=1:lx-1
        if (x(k)>(mo+r*h))
    xn(k)=(x(k)+x(k+1))/2;
    yn(k)=(y(k)+y(k+1))/2;
        else
            if (x(k)<(mo-r*h))
            xn(k)=(x(k)+x(k+1))/2;
            yn(k)=(y(k)+y(k+1))/2;
            end
        end
end
yn=y; xn=x;
for k=1:lx-1
    if (abs(x(k)-mo)<abs(r*h))
    if (k>1)
        xn(k)=xn(k-1);
        yn(k)=yn(k-1);
    else
        xn(k)=xn(k+1);
        yn(k)=yn(k+1);
    end
    end
end
for k=1:lx
    v=(x(k)-mo)^2; u=si^2;
    ys(k)=exp(-v/u/2)/si/sqrt(2*pi);
end
t=postrgraf(xn,yn,x,ys,maxi,mini);
end

function t = postrgraf(xn,yn,x,ys,maxi,mini)
pl=plot(xn,yn,'.g',x,ys,'.k');
xlabel('Размер частиц, мкм'); ylabel('Плотность распределения, 1/мкм');
%ylim([0 1]);
xlim([mini maxi]);
set(pl,'LineWidth',2); hold on; grid on;
title({'y(x)'});
%legend('y(x)','exp(-x^2)','location','best');
t=0;
end

function n = poisknomera(x,lx,zx)
q=0; na=1; nb=lx; ra=abs(na-nb); h=1; nh=1e2;
while ((ra>1) && (h<nh))
    nc=round((na+nb)/2);
    fa=x(round(na))-zx;
    fb=x(round(nb))-zx;
    fc=x(round(nc))-zx;
    no=VybCAB(fa,fb,fc,na,nb,nc);
na=round(no(1)); nb=round(no(2)); nc=round(no(3));
ra=abs(na-nb); h=h+1;
end
    if (x(nc)<zx) 
        if (nc>1)
            q=nc-1;
        else q=1;
        end
    else q=nc;
    end
    if ((q==nc) && (q>1))
    %disp(nc-1); disp(x(nc-1));
    %disp(nc); disp(x(nc));
    end
    if ((q>nc) && (q<lx))
        %disp(q-1); disp(x(q-1));
        %disp(q); disp(x(q));
    end
n=q;
end

function [ vy ] = VybCAB(fa,fb,fc,xa,xb,xc)
if (((fa*fc)<0) && ((fb*fc)>0))
    xb=xc; 
end
if (((fa*fc)>0) && ((fb*fc)<0))
    xa=xc; 
end
vy = [xa,xb,xc];
end