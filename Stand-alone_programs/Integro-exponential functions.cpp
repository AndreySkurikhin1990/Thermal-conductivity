#include <stdio.h>
#include <conio.h>
#include <math.h>
double integexpon(int,double);
double integexpon(int n, double x)
{
int N=1e3, m, k;
double E=0, s=0, t, e;
    for (m=0; m<=N; m++) 
	{
		if (m==(n-1))
		 continue;  
        t=1;
        for (k=1; k<=m; k++)
		{ t=t*k; }
        s=s+pow((-x),m)/((m-n+1)*t);
	}
e=0.577215665057043;
t=1;
for (k=1; k<(n-1); k++)
	t=t*k; 
E=pow((-x),(n-1));
E=E/t;
t=0;
for (k=1; k<=(n-1); k++)
 t=t+1/k; 
psi=-e+t;
E=E*(-log(x)+psi);
E=E-s;
return(E);
}
void main()
{
double h=0.01, ko=0.1, ie;
int i, n=2, x=0;
for (i=0;i<ko;i++)
{
	ie=integexpon(n,x);
	printf("%lf",ie);
		x=x+h;
}
}