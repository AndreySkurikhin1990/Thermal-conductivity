%format long g;
dl=dlvoVer53101;alfs=SredGraf; npp=Kramers_n();lnpp=length(npp);te0=273.15;T0=1000+te0; Tk=235+te0;tol=33.5;y0=30;Tz=(T0-Tk)*abs(tol-y0)/y0/2+T0;
for j=1:lnpp
    dl(j)=1e-2/dl(j);pp(j)=dl(j)*alfs(j)/4/pi;emsic=aSiCemmis(dl(j));R2n(j)=((emsic(2)-1)^2+emsic(3))/((emsic(2)+1)^2+emsic(3));R2m(j)=1-emsic(1);epsa(j)=emmet(npp(j),pp(j));
end
R1=1-epssred(dl,epsa,Tz)
R2=kootsre(dl,R2m,Tz)
%R2=kootsre(dl,R2n,Tz);
%p=plot(dl*1e6,R2n*1e2,'-b');set(p,'LineWidth',2); hold on; grid on; xlabel({'Длина волны, мкм'}); ylabel({'Нормальный коэффициент отражения Rn, %'}); title({'График зависимости коэффициента отражения Rn при нормальном падении луча от длины волны'});
%np=real(nsreddv(dl,npp,Tz));alf=1/sprko(Tz);ppo=alf*(dl(lnpp)+dl(1))/2;R1=1-emmet(np,ppo);
%Две одинаковые параллельные пластины, расположенные друг напротив друга, длина a, ширина b, расстояние между ними c
a=114; b=114; c=33.5; X=a/c; Y=b/c; ab=a;h=c;
F12kpkp=(0.5*log((1+X^2)*(1+Y^2)/(1+X^2+Y^2))+X*sqrt(1+Y^2)*atan(X/sqrt(1+Y^2))+Y*sqrt(1+X^2)*atan(Y/sqrt(1+X^2))-X*atan(X)-Y*atan(Y))*2/pi/X/Y;
ah=a^2+h^2;F122kv=(2/pi)*(2*sqrt(ah)*atan(a/sqrt(ah))/a-2*h*atan(a/h)/a+0.5*log(ah^2/h^2/(ah+a^2))*(h/a)^2);
%Интегрирование по полосам вровень с цилиндром равной длины
H  = 65; L  = 150; R  = 4; 
%IAe0c(1)=H; IAe0c(2)=L; IAe0c(3)=R; 
HH = H/R; LL = L/R; X  = (1+HH)^2+LL^2; Y  = (1-HH)^2+LL^2; 
F12el0cyl = LL*(atan(LL/sqrt(HH^2-1))/LL+((X-2*HH)/sqrt(X*Y))*atan(sqrt(X*(HH-1)/(HH+1)/Y))-atan(sqrt((HH-1)/(HH+1))))/pi/HH;
%F12el0cyl=view(15,3,IAe0c);
%Интегрирование по полосам нулевого ряда
H=150; S=65; X=0;
HH = H/R; SS = S/R; XX = X/R; C=SS^2+XX^2; CC = sqrt(C); A = HH^2+C-1; B = HH^2-C+1; 
F12pol0cyl= (SS/C)*(1-(acos(B/A)-(sqrt(A^2+(2*HH)^2)*acos(B/A/CC)+B*asin(1/CC))/HH/2)/pi-A/HH/4);
%IAp0c(1)=H; IAp0c(2)=R; IAp0c(3)=S; IAp0c(4)=X; F12pol0cyl=view(28,4,IAp0c);
%Интегрирование по полосам первого ряда
X=21; HH = H/R; SS = S/R; XX = X/R; C=SS*SS+XX*XX; CC = sqrt(C); 
A = HH^2+C-1; B = HH^2-C+1; 
F12pol1cyl= (SS/C)*(1-(acos(B/A)-(sqrt(A^2+(2*HH)^2)*acos(B/A/CC)+B*asin(1/CC))/HH/2)/pi-A/HH/4);
%IAp1c=IAp0c; IAp1c(4)=X;F12pol1cyl=view(28,4,IAp1c);
%Интегрирование по полосам второго ряда
X=21*2; HH = H/R; SS = S/R; XX = X/R; C=SS*SS+XX*XX; CC = sqrt(C); 
A = HH*HH+C-1; B = HH*HH-C+1;
F12pol2cyl=(SS/C)*(1-(acos(B/A)-(sqrt(A^2+(2*HH)^2)*acos(B/A/CC)+B*asin(1/CC))/HH/2)/pi-A/HH/4);
%IAp2c=IAp1c; IAp2c(4)=X;F12pol2cyl=view(28,4,IAp2c);
X=0; HH = H/R;SS = S/R;XX = X/R;Y=10; YY=Y/R;B=XX^2+SS^2;A=XX^2+SS^2+YY^2;C=(HH-YY)^2;CB=(C-B+1)/(C+B-1);YBA=(YY^2-B+1)/(A-1);
s1=SS/B; s2=acos(YBA); s3=acos(CB); s4=YY*(A+1)*acos(YBA/sqrt(B))/sqrt((A-1)^2+(2*YY)^2);s5=sqrt(C)*(C+B+1)*acos(CB/sqrt(B))/sqrt((C+B-1)^2+4*C);s6=HH*acos(1/sqrt(B));
F12elXYcyl=s1-(s2+s3-s4-s5+s6)*s1/2/pi;
hx=10e-2;hy=10e-2; Y0=18; Yk=ab+Y0; Xa=0; Ya=0; Ya=Y0:hy:Yk; ly=length(Ya);X0=0;deX=21;
R=4;H=150; HH = H/R;S=65;SS = S/R;
for p=1:3
    Xa=X0-ab/2:hx:X0+ab/2; lx=length(Xa);
    %suy=[0 0 0 0];
for j=1:ly
    YY=Ya(j)/R; 
    %sux=[0 0 0 0];
for k=1:lx
XX = Xa(k)/R;B=XX^2+SS^2;A=XX^2+SS^2+YY^2;C=(HH-YY)^2;CB=(C-B+1)/(C+B-1);YBA=(YY^2-B+1)/(A-1);
s1=SS/B; s2=acos(YBA); s3=acos(CB); s4=YY*(A+1)*acos(YBA/sqrt(B))/sqrt((A-1)^2+(2*YY)^2);s5=sqrt(C)*(C+B+1)*acos(CB/sqrt(B))/sqrt((C+B-1)^2+4*C);s6=HH*acos(1/sqrt(B));        
F12xy=s1-(s2+s3-s4-s5+s6)*s1/2/pi; F11xy=1-F12xy; 
F21xy=hx*hy*F12xy/(2*pi*R*H); F22xy=1-F21xy;
Ma=0; 
Ma(1,1)=1-R1*F11xy; Ma(1,2)=-R2*F21xy; 
Ma(2,1)=-R1*F12xy; Ma(2,2)=1-R2*F22xy; 
Ma(3,3)=1-R1*F11xy;Ma(3,4)=-R2*F21xy; 
Ma(4,3)=-R1*F12xy; Ma(4,4)=1-R2*F22xy; 
sb(1)=F11xy; sb(2)=F12xy; 
sb(3)=F21xy; sb(4)=F22xy;
t=(inv(Ma)*sb')'; 
t(4)=(1-(1-R1)*abs(t(3)))/(1-R2); 
t(2)=(1-(1-R1)*abs(t(1)))/(1-R2);
for r=1:4  
    Fb(p,r,j,k)=t(r);
    %sux(r)=sux(r)+t(r)*hx; 
end;
    end
%for r=1:4  suy(r)=suy(r)+sux(r)*hy; end;
end
F3b=0;
for r=1:4
%Fb(p,r)=suy(r)/ab^2;
end
X0=X0+deX;
end
for p=1:3
    for r=1:4
sy=0;
for j=1:ly
    sx=0;
    for k=1:lx
       sx=sx+Fb(p,r,j,k)*hx;
    end
    sy=sy+sx*hy;
end
F3b(p,r)=sy/ab^2;
    end
end
disp(F3b);