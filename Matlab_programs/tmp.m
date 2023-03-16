format long g;
rad=0.3/2;
alp=150;
lamb=164;
len=0.24/2;
t0=400;
tl=15;
ro=2787;
cp=833;
tau1=0.4*60;
at=lamb/(ro*cp);
For1=at*tau1/rad^2;
Fol1=at*tau1/len^2;
Bir=rad*alp/lamb;
Bil=len*alp/lamb;
N=1e4; eps=1e-10; h=1e-2; j=1; nac=0; kon=0;
for k=0:N/h
    a=k*h; b=(k+1)*h;
    fa=besselj(0,a)/(besselj(1,a)*a)-1/Bir;
    fb=besselj(0,b)/(besselj(1,b)*b)-1/Bir;
if (fa*fb<0) nac(j)=a; kon(j)=b; j=j+1; end;
end;
r=length(nac);
mu=0;
for j=1:r
    a=nac(j);
    b=kon(j);
    t=1; f=1;
    while (abs(b-a)>eps)
    c=(a+b)/2;
    fa=besselj(0,a)/(besselj(1,a)*a)-1/Bir;
    fb=besselj(0,b)/(besselj(1,b)*b)-1/Bir;
    fc=besselj(0,c)/(besselj(1,c)*c)-1/Bir;
    if (fc==0) break; end; 
    if (fa==0) 
        c=a; break; end; 
    if (fb==0) 
        c=b; break; end;
    if (fa*fc<0) 
        if (fb*fc>0) b=c; end; end; 
    if (fa*fc>0)
        if (fb*fc<0) a=c; end; end;
    if (t>1e3) 
        f=0; break;
    end;
    t=t+1;
    end;
    if (f==1) 
        mu(j)=c; end;
end;
j=1; muk=0; 
for k=1:length(mu)
    a=mu(k);
    m=besselj(0,a)/(besselj(1,a)*a)-1/Bir;
    if (abs(m)<1e-6)   muk(j)=a;    j=j+1; else continue; end;
end
r=length(muk); d=0;
for k=1:r
    a=muk(k);
    d(k)=2*besselj(1,a)/(a*(besselj(0,a)^2+besselj(1,a)^2));
end
s=0; xc1=0; s1=0; xg1=rad;
for k=1:r
    a=muk(k);
    s=s+d(k)*besselj(0,a*xc1/rad)*exp(-(a^2)*For1);
    s1=s1+d(k)*besselj(0,a*xg1/rad)*exp(-(a^2)*For1);
end
nacl=0;
konl=0; j=1;
for k=0:N
    nacl(j)=k*pi+1e-8; 
    konl(j)=k*pi+pi/2-1e-8;
    j=j+1;
end;
r=length(nacl);
mul=0;
for j=1:r
    a=nacl(j);
    b=konl(j);
    t=1; f=1;
    while (abs(b-a)>eps)
    c=(a+b)/2;
    fa=a/Bil-cot(a);
    fb=b/Bil-cot(b);
    fc=c/Bil-cot(c);
    if (fc==0) break; end; 
    if (fa==0) 
        c=a; break; end; 
    if (fb==0) 
        c=b; break; end;
    if (fa*fc<0) 
        if (fb*fc>0) b=c; end; end; 
    if (fa*fc>0)
        if (fb*fc<0) a=c; end; end;
    if (t>1e4) 
        f=0; break;
    end;
    t=t+1;
    end;
    if (f==1) 
        mul(j)=c; end;
end;
j=1; mukl=0; r=length(mul); 
for k=1:r
    a=mul(k);
    m=a/Bil-cot(a);
    if (isinf(abs(m))) continue; end;
    if (abs(m)<1e-6)   
        mukl(j)=a
        j=j+1; 
    else continue; 
    end;
end;
r=length(mukl); cl=0;
for k=1:r
    a=mukl(k);
    cl(k)=2*sin(a)/(a+sin(a)*cos(a));
end
s2=0; x01=0; s3=0; xl2=len/2;
for k=1:r
    a=mukl(k);
    s2=s2+cl(k)*cos(a*x01/len)*exp(-(a^2)*Fol1);
    s3=s3+cl(k)*cos(a*xl2/len)*exp(-(a^2)*Fol1);
end
Nl=1e4; Nr=Nl;
kool=0:len/Nl:len;
koor=0:rad/Nl:rad;
r=length(mukl);
t=length(kool);
p=length(koor);
for j=1:t
st1=0; st2=0;
for k=1:r
    a=mukl(k);
    st1=st1+cl(k)*cos(a*kool(j)/len)*exp(-(a^2)*Fol1);
end
st2=st1; st1=s*st1; st2=s1*st2;
temxc1(j)=st1*(t0-tl)+tl;
temxg1(j)=st2*(t0-tl)+tl;
end
r=length(muk);
for j=1:p
st3=0; st4=0;
for k=1:r
    a=muk(k);
    st3=st3+d(k)*besselj(0,a*koor(j)/rad)*exp(-(a^2)*For1);
end
st4=st3; st3=s2*st3; st4=st4*s3;
temx01(j)=st3*(t0-tl)+tl;
temxl2(j)=st4*(t0-tl)+tl;
end
ti=0:0.1:1.5e3;
p=length(ti);
temt=0;
for k=1:p
Fort=at*ti(k)/rad^2;
Folt=at*ti(k)/len^2;
r=length(muk); s1=0;
for j=1:r
    a=muk(j);
    s1=s1+d(j)*besselj(0,0)*exp(-(a^2)*Fort);
end
r=length(mukl); s2=0;
for j=1:r
    a=mukl(j);
    s2=s2+cl(j)*cos(0)*exp(-(a^2)*Folt);
end
temt(k)=s1*s2*(t0-tl)+tl;
end
q1=plot(kool,temxc1,'-b');
set(q1,'LineWidth',3);
xlabel('Координата, м');
ylabel('Температура от x на оси, град С');
title('График T(x,0,tau1)');
grid on;
figure;
q2=plot(kool,temxg1,'-b');
set(q2,'LineWidth',3);
xlabel('Координата, м');
ylabel('Температура от x на границе цилиндра, град С');
title('График T(x,r0,tau1)');
grid on;
figure;
q3=plot(koor,temx01,'-b');
set(q3,'LineWidth',3);
xlabel('Координата, м');
ylabel('Температура от r в середине цилиндра, град С');
title('График T(0,r,tau1)');
grid on;
figure;
q4=plot(koor,temxl2,'-b');
set(q4,'LineWidth',3);
xlabel('Координата, м');
ylabel('Температура от r в торце цилиндра, град С');
title('График T(L/2,r,tau1)');
grid on;
figure;
q5=plot(ti,temt,'-b');
set(q5,'LineWidth',3);
xlabel('Время, с');
ylabel('Температура в центре от времени, град С');
title('График T(0,0,tau)');
grid on;