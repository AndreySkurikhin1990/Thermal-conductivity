function t = RasExpSVodoy()
te0=27;
mnc=PolMasTemCen();
mn08=PolMasTemNaRas08();
n=length(mnc);
mtr=MasTemKhromelKopelTermop();
mnr=MasNaprKhromelKolpelTermop();
vr=PolMasVrem();
mtc=zeros(1,n);
mt08=zeros(1,n);
for k=1:n
    mtc(k)=opredTochZnach(mtr, mnr, mnc(k))+te0;
    mt08(k)=opredTochZnach(mtr, mnr, mn08(k))+te0;
end
kp1=koefPribLin(vr, mtc, mt08)
t=0;
end

function [ m ] = PolMasTemCen()
m=[0.36
0.36
0.36
0.36
0.36
0.36
0.36
0.36
0.36
0.37
0.37
0.37
0.37
0.38
0.38
0.39
0.39
0.4
0.4
0.41
0.41
0.42
0.43
0.43
0.44
0.44
0.45
0.46
0.47
0.48
0.49
0.5
0.51
0.52
0.53
0.53
0.54
0.55
0.56
0.56
0.57
0.58
0.58
0.6
0.6
0.61
0.62
0.63
0.63
0.64
0.65
0.66
0.66
0.67
0.68
0.69
0.7
0.71
0.72
0.73
0.73
0.74
0.75
0.76
0.77
0.78
0.79
0.8
0.81
0.81
0.82
0.82
0.83
0.83
0.84
0.84
0.85
0.86
0.86
0.87
0.87
0.88
0.88];
end

function [ m ] = PolMasTemNaRas08()
m=[0.32
0.32
0.32
0.32
0.34
0.35
0.37
0.39
0.41
0.43
0.45
0.47
0.48
0.49
0.51
0.53
0.54
0.57
0.57
0.59
0.61
0.61
0.62
0.63
0.65
0.66
0.69
0.69
0.69
0.72
0.72
0.72
0.73
0.75
0.76
0.78
0.78
0.78
0.79
0.8
0.81
0.82
0.82
0.83
0.84
0.84
0.85
0.86
0.86
0.87
0.87
0.87
0.89
0.9
0.9
0.92
0.91
0.92
0.92
0.93
0.94
0.95
0.95
0.95
0.96
0.96
0.97
0.97
0.98
0.98
0.99
1
1
1
1.01
1.01
1.02
1.03
1.03
1.03
1.04
1.04
1.04];
end

function [ m ] = PolMasVrem()
m=[0
10
20
30
40
50
60
70
80
90
100
110
120
130
140
150
160
170
180
190
200
210
220
230
240
250
260
270
280
290
300
310
320
330
340
350
360
370
380
390
400
410
420
430
440
450
460
470
480
490
500
510
520
530
540
550
560
570
580
590
600
610
620
630
640
650
660
670
680
690
700
710
720
730
740
750
760
770
780
790
800
810
820];
end

function [ m ] = MasTemKhromelKopelTermop()
m=[0
1
2
3
4
5
6
7
8
9
10
11
12
13
14
15
16
17
18];
end

function [ m ] = MasNaprKhromelKolpelTermop()
m=[0
0.063
0.127
0.19
0.254
0.318
0.382
0.446
0.51
0.575
0.639
0.704
0.768
0.833
0.898
0.963
1.028
1.093
1.159];
end

function ktp = opredTochZnach(ktptks, te, temp)
f=1; p=0; nn=1; n=length(te); e=1e-6; 
if ((temp>=te(nn)) && (temp<=te(n))) 
    for k = 1:n
        if ((te(k) >= temp) && (f>0)) 
            p = k; f = 0; 
            break;
        end
    end
if (p==0) 
    if (f==0) 
            p = 2; 
    else
            p = n; f = 0;
    end
end
elseif (temp<te(nn)) 
        p=2; f=0;
elseif (temp>te(n)) 
        p=n; f=0;
end
if ((f==0) && (p>nn))
		x2=te(p);
		x1=te(p - 1);
		dt = x2 - x1;
            if (abs(dt) > e)
                y2=ktptks(p);
                y1=ktptks(p - 1);
                b=y1;
			if ((n==p) && (temp>te(p))) 
                b=y2;
            elseif ((n==p) && (temp<=te(p))) 
                b=y1;
            end
                ko = (y2 - y1) / dt;
                ktp = b + ko*(temp - x1);
            end
        else
            ktp=0.0;
end
end
                
function [ ko ] = koefPrib(ktp, te)
    le=length(ktp);
	e=1e-15; yx2=0.0; yx=0.0; p=0.0; hf=1e0;
	x4=0.0; x3=x4; x2=x3;
    x=sum(te);
    y=sum(ktp);
	for k=1:le
        t=te(k);
        r=ktp(k);
		yx2 = yx2 + r * (t^2); 
		yx = yx + r * t; 
		x4 = x4 + (t^4); 
		x3 = x3 + (t^3); 
		x2 = x2 + (te^2); 
    end
	b = [yx2, yx, y];
	A=[x4,x3,x2; x3,x2,x; x2, x, le];
	ko=inv(A)*b;
end

function [ ko ] = koefPribLin(tv, tc, t08)
n=length(tv);
x2=zeros(1,n);
xy1=zeros(1,n);
xy2=zeros(1,n);
for k=1:n
    x2(k)=tv(k)^2;
    xy1(k)=tv(k)*tc(k);
    xy2(k)=tv(k)*t08(k);
end
sx2=sum(x2);
sx=sum(tv);
sy1=sum(tc);
sy2=sum(t08);
sxy1=sum(xy1);
sxy2=sum(xy2);
k1=(n*sxy1-sx*sy1)/(n*sx2-(sx^2));
k2=(n*sxy2-sx*sy2)/(n*sx2-(sx^2));
b1=(sy1-k1*sx)/n;
b2=(sy2-k2*sx)/n;
ko=[k1,k2,b1,b2];
end