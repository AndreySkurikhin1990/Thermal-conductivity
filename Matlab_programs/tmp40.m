function t = tmp40()
phip=pi/6;
a=4; b=3;
xp=2; yp=1;
r=PoiskDlPutiLuchaOsn(phip,a,b,xp,yp)
t=0;
end
%phip - угол направления падения, xp, yp - координаты падения на эллипсе, a, b - полуоси
function dlot = PoiskDlPutiLuchaOsn(phip,a,b,xp,yp)
ro=0; ep=1e-6; Nit=1e3; ko=tan(phip); sc=yp-ko*xp; zn=1;
for j=1:2
zn=-zn; k=0; xa=min([0 zn*a]); xb=max([0 zn*a]); ra=abs(xa-xb); 
while ((ra>ep) && (k<Nit))
    xc=(xa+xb)/2; 
    fa=zn*a*sqrt(1-(xa/b)^2)-(ko*xa+sc); 
    fb=zn*a*sqrt(1-(xb/b)^2)-(ko*xb+sc); 
    fc=zn*a*sqrt(1-(xc/b)^2)-(ko*xc+sc); 
    ro=vybPoKo(fa,fb,fc,xa,xb,xc); 
    xa=ro(1); xb=ro(2); xc=ro(3); k=k+1;
    ra=abs(xa-xb);
end
xpe(j)=xc; ype(j)=ko*xc+sc;
end
v1=[xpe(1)-xp ype(1)-yp]; v2=[xpe(2)-xp ype(2)-yp];
xna=b*cos(phip); yna=a*sin(phip); 
if ((sign(v1(1))==sign(xna)) && (sign(v1(2))==sign(yna)))
    dlo=sqrt(v1(1)^2+v1(2)^2);
end
if ((sign(v2(1))==sign(xna)) && (sign(v2(2))==sign(yna)))
    dlo=sqrt(v2(1)^2+v2(2)^2);
end
dlot=dlo;
end
function [ v ] = vybPoKo(fa,fb,fc,a,b,c)
if ((fc*fb)>0) 
        if ((fa*fc)<0) 
            b=c; 
        end
    end
    if ((fc*fa)>0) 
        if ((fb*fc)<0) 
            a=c; 
        end
    end
    v=[a b c];
end