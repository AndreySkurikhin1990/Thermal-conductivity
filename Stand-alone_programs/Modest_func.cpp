#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <iostream>
using namespace std;
const double pi=acos(-1e0);
double *GRAYDIFFSPEC(int, double *, double *, double *, double *, double **, int *, double *, double *, int);
double *GRAYDIFF(int, double *, double *, double *, double **, int *, double *, double *, int);
double VIEW(int, int, double *, int);
double perpplates(double, double, double, double, double, double, double);
double parlplates(double, double, double, double, double, double, double);
double *GAUSS(double **, double *, int, double *);
double SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *, int);
double SeryeStenkiRasIzlDifPov(double, double, double, double *, double *, double *, double *, double *, int *, int);
double *SEMIGRAY(int, double *, double **, double **, double **, double ***, int *, double *, double *, double *, double *, int, int);
double *MetodGaussa(double **, double *, int, double *);
void main()
{
int f=6, k=2, *nost=new int[k];
double sigma=5.67e-8, w=0.4, h=0.3, l=1e6, ra=0.0, *RHO=new double[f];
double *A=new double[f], *T=new double[f], *HO=new double[f], *EPS=new double[f];
A[0]=w*l; A[2]=w*l; A[1]=h*l; A[3]=h*l; A[4]=h*w; A[5]=h*w; k=0; nost[k]=k; k++; nost[k]=2;
for (k=0; k<f; k++) HO[k]=0.0;
RHO[0]=0.7; RHO[2]=0.0; RHO[1]=0.2; RHO[3]=0.0; RHO[4]=0.0; RHO[5]=0.0;
EPS[0]=0.3; EPS[2]=0.3; EPS[1]=0.8; EPS[3]=0.8; EPS[4]=0.8; EPS[5]=0.8;
T[0]=1e3; T[2]=1e3; T[1]=6e2; T[3]=6e2; T[4]=6e2; T[5]=6e2;
ra=SeryeStenkiRasIzl(w, h, l, T, EPS, HO, RHO, A, nost, f); 
k=getchar();
}
double SeryeStenkiRasIzl(double ww, double hh, double ll, double *T, double *EPS, double *HOs, double *RHOs, double *A, int *ns, int chst) //N серых диффузно излучающих поверхностей с диффузным и зеркальным коэффициентами отражения
{ int N=chst, fl=0; double b=0.0, t=0.0, epsi=1e-10; //t=T[ns[0]]; b=T[ns[1]]; if (fabs(t-b)<epsi) return 0.0;
double **Fs=new double*[N], z=0.0, r=0.0, y=0.0; int i=0, j=0, f=3, k=0; 
double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p=NULL;
for (i=0; i<N; i++) { p=new double[N]; for (j=0; j<N; j++) p[j]=0.0; Fs[i]=p; }
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd=0; for (i=0; i<N; i++) id[i]=1; // Поверхности 1 и 3 (bottom and top)
    	for (i=0; i<N; i++) if (!(id[i])) PIN[i]=q[i]; else PIN[i]=sigma*pow(T[i], 4e0); //Заполнение массива PIN ПТП q и температурой T - перевод температур в ППИ, Fs(1,2)=F1-2, 
	iclsd=1; //угловые коэффициенты; конфигурация замкнута (iclsd=1), диагональные элементы не нужны //for (j=0; j<N; j++) cout << "RHOs ( " << j << " ) = " << RHOs[j] << endl;  for (j=0; j<N; j++) cout << "EPS ( " << j << " ) = " << EPS[j] << endl;
	k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(ww); //h,l,w
	k=39; i=0; j=1; Fs[i][j]=VIEW(k, f, arg, fl); //cout << "Fs(1,2) = " << Fs[0][1] << "\th = " << hh << "\tl = " << ll << "\tw = " << ww << "\tf = " << f << endl; //Fs(1,2)=F12
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh);     // a,b,c
	k=38; r=VIEW(k, f, arg, fl); 
	z=0.0; t=parlplates(ww, ww, 2.0*ww, ll, z, ll, hh); k=1; b=RHOs[k];
    i=0; j=2; Fs[i][j]=r+b*t; //Fs13=F13+rhos2*F1(2)3	
	t=perpplates(ww, 2.0*ww, z, hh, ll, z, ll); k=1; b=RHOs[k];
    i=0; j=3; k=1; Fs[i][j]=Fs[i][k]+b*t; //cout << "t = " << t << "\tb = " << b << endl; //Fs14=F14+rhos2*F1(2)4
    i=0; j=4; k=1; Fs[i][j]=0.5*(1.0-(1.0-RHOs[k])*Fs[i][k]-Fs[i][2]-Fs[i][3]); //смещены два прямоугольника на ww по оси X
    i=0; j=5; k=4; Fs[i][j]=Fs[i][k]; 
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); //h,l,w
    k=39; r=VIEW(k, f, arg, fl); 
	k=0; b=RHOs[k];
	t=perpplates(hh, 2.0*hh, z, ww, ll, z, ll);
	i=1; j=2; Fs[i][j]=r+b*t; //Fs23=F23+rhos1*F1(1)3 //смещены два прямоугольника на hh по оси X, т.к. у них разные размеры
	k=0; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ww); //a,b,c
    k=38; r=VIEW(k, f, arg, fl);
	k=0; b=RHOs[k];
	t=parlplates(hh, hh, 2.0*hh, ll, z, ll, ww);
	i=1; j=3; Fs[i][j]=r+b*t;  //Fs24=F24+rhos1*F2(1)4
    i=1; j=4; k=0; r=A[k]*Fs[k][i]/A[i];
	Fs[i][j]=0.5*(1.0-(1.0-RHOs[k])*r-Fs[i][2]-Fs[i][3]); //Fs25
    i=1; j=5; k=4; Fs[i][j]=Fs[i][k]; //Fs26
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=2.0*fabs(hh);  // a,b,c
	k=38; r=VIEW(k, f, arg, fl);
	k=1; b=RHOs[k];
	t=parlplates(ww, ww, 2.0*ww, ll, z, ll, 2.0*hh);
    i=2; j=0; Fs[i][i]=RHOs[j]*(r+b*t); //Fs33=rhos1*F3(1)3+rhos1*rhos2*F3(12+21)3 //излучение стенки саму на себя
	r=perpplates(z, ww, hh, 2.0*hh, ll, z, ll);
	k=1; b=RHOs[k];
	t=perpplates(ww, 2.0*ww, hh, 2.0*hh, ll, z, ll);
    i=2; j=3; k=0; Fs[i][j]=Fs[k][j]+RHOs[k]*(r+b*t); // Fs34=Fs1-4+rhos1*F3(1)-3+rhos1*rhos2*F3(12+21)-4
    i=2; j=4; k=0; r=A[k]*Fs[k][i]/A[i]; t=A[1]*Fs[1][i]/A[i];
	Fs[i][j]=0.5*(1.0-(1.0-RHOs[k])*r-(1.0-RHOs[1])*t-Fs[i][i]-Fs[i][3]);
    i=2; j=5; k=4; Fs[i][j]=Fs[i][k];
    k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=2.0*fabs(ww); //a,b,c
	k=38; r=VIEW(k, f, arg, fl);
	k=0; b=RHOs[k];
	t=parlplates(hh, hh, 2.0*hh, ll, z, ll, 2.0*ww);
    i=3; j=1; Fs[i][i]=RHOs[j]*(r+b*t); //Fs44=rhos2*F4(2)-4+rhos1*rhos2*F4(12+21)-4
    i=3; j=4; k=0; r=A[k]*Fs[k][i]; t=A[1]*Fs[1][i]; b=A[2]*Fs[2][i];
	Fs[i][j]=0.5*(1.0-((1.0-RHOs[k])*r+(1.0-RHOs[1])*t+b)/A[i]-Fs[i][i]);
    i=3; j=5; k=4; Fs[i][j]=Fs[i][k];
    i=4; j=5; k=0; r=A[k]*Fs[k][i]; t=A[1]*Fs[1][4]; b=A[2]*Fs[2][4]+A[3]*Fs[3][4];
	Fs[i][j]=0.5*(1.0-((1.0-RHOs[k])*r+(1.0-RHOs[1])*t+b)/A[4]); //Fs56 //for (j=0; j<N; j++) cout << "A ( " << j << " ) = " << A[j] << endl; //for (i=0; i<N; i++) { y=RHOs[i]+EPS[i]; cout << "Edinitsa = " << y << endl; } //for (i=0; i<N; i++) cout << "Temp ( " << i << " ) = " << T[i] << endl; //POUT=GRAYDIFFSPEC(iclsd, A, EPS, RHOs, HOs, Fs, id, PIN, POUT, N); //Решение системы уравнений методом Гаусса //for (i=0; i<N; i++) cout << POUT[i] << endl; 
	POUT=GRAYDIFFSPEC(iclsd, A, EPS, RHOs, HOs, Fs, id, PIN, POUT, N);
	sumq=0.0; //проверка полной ППИ
    	for (i=0; i<N; i++) //Вывод - перевод в температуры
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; 
			cout << "Surface " << i << "\tT [K] = " << T[i] << "\tq [W/m2] = " << q[i] << "\tQA [W] = " << QA[i]/ll << endl;
		} r=0.0; z=1.0; for (i=0; i<N; i++) if (ns[i]==i) { r=r+z*q[i]; z=-z; } 
		for (i=0; i<N; i++) { p=Fs[i]; delete []p; } delete []Fs; delete []PIN; 
		delete []POUT; delete []arg; delete []q; delete []QA; delete []id; 
		return r; }
double *GRAYDIFF(int iclsd, double *A, double *EPS, double *HOs, double **Fs, int *ID, double *PIN, double *POUT, int chst)
{
int N=chst, i, j;
double **qm=new double*[N], **em=new double*[N], **pm=new double*[N], *B=new double[N], *p, *idf=new double[N], hf=1e0, s, ikr, epsi=1e-10;
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; p=new double[N]; em[i]=p; p=new double[N]; pm[i]=p; }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<i; j++)
		  if (fabs(A[i])>epsi) 
			  Fs[i][j]=(A[j]/A[i])*Fs[j][i]; else Fs[i][j]=0.0; //Для закрытой конфигурации, необходимо рассчитать диагональные элементы из правила суммирования
	if (iclsd==1) { for (i=0; i<N; i++) {
		Fs[i][i]=1e0;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-Fs[i][j]; } } } }
	 for (i=0; i<N; i++) { // Заполнение матриц коэффициентов q и e
		for (j=0; j<N; j++) {
			if (i==j) ikr=1e0; else ikr=0.0; // Символ Кронекера delta_ij //cout << "i = " << i << "\tj = " << j << "\tk = " << ikr << endl;
			if (fabs(EPS[j])>epsi) 
				qm[i][j]=ikr/EPS[j]-(1e0/EPS[j]-1e0)*Fs[i][j];
			em[i][j]=ikr-Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Заполнение матрицы выходных коэффициентов POUT и RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(1e0-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(1e0-idf[j]))*PIN[j]; } }
	 POUT=MetodGaussa(pm, B, N, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=GAUSS(pm, B, N, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << "\t"; cout << endl;
	 double sumq = 0.0, *T=new double[N], sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	 for (i=0; i<N; i++) { if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0); else q[i]=POUT[i];
     QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; } //printf("Summa vsekh potokov = %0.10lf\n", sumq); 
	 delete []T; delete []q; delete []QA;
	 for (i=0; i<N; i++) { for (j=0; j<N; j++) printf("%0.6lf\t", Fs[i][j]); cout << endl; }
	 for (i=0; i<N; i++) { p=qm[i]; delete []p; p=em[i]; delete []p; p=pm[i]; delete []p; } 
	 delete []qm; delete []em; delete []pm; delete []B; delete []idf; return POUT;
}
double *GRAYDIFFSPEC(int iclsd, double *A, double *EPS, double *RHOs, double *HOs, double **Fs, int *ID, double *PIN, double *POUT, int chst)
{ int N=chst, i=0, j=0; 
double **qm=new double*[N], **em=new double*[N], **pm=new double*[N], *B=new double[N], *p, *idf=new double[N], hf=1e0, s, ikr, y, epsi=1e-10; // Вычисление недостающих угловых коэффициентов - левой нижней части матрицы из свойства взаимности
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; p=new double[N]; em[i]=p; p=new double[N]; pm[i]=p; }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<i; j++)
		  if (fabs(A[i])>epsi) 
			  Fs[i][j]=(A[j]/A[i])*Fs[j][i]; else Fs[i][j]=0.0; //Для закрытой конфигурации, необходимо рассчитать диагональные элементы из правила суммирования
	if (iclsd==1) { for (i=0; i<N; i++) {
		Fs[i][i]=1e0;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-(1e0-RHOs[j])*Fs[i][j]; } }
		if (RHOs[i]<1e0) Fs[i][i]=Fs[i][i]/(1e0-RHOs[i]); else RHOs[i]=0.0; } }
	 for (i=0; i<N; i++) { // Заполнение матриц коэффициентов q и e
		for (j=0; j<N; j++) {
			if (i==j) ikr=1e0; else ikr=0.0; // Символ Кронекера delta_ij //cout << "i = " << i << "\tj = " << j << "\tk = " << ikr << endl;
			if (fabs(EPS[j])>epsi) 
				qm[i][j]=ikr/EPS[j]-((1.0-RHOs[j])/EPS[j]-1.0)*Fs[i][j];
			em[i][j]=ikr-(1.0-RHOs[j])*Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Заполнение матрицы выходных коэффициентов POUT и RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(1e0-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(1e0-idf[j]))*PIN[j]; } } //cout << "Fs" << endl; for (i=0; i<N; i++) { for (j=0; j<N; j++) printf("%0.4lf\t",Fs[i][j]); cout << endl; } //for (i=0; i<N; i++) { y=0.0; for (j=0; j<N; j++) y=y+(1e0-RHOs[j])*Fs[i][j]; cout << "stroka ( " << i << " ) = " << y << endl; } //j=getchar();
	 POUT=MetodGaussa(pm, B, N, POUT); //for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=GAUSS(pm, B, N, POUT); for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl; POUT=reshMetObrMatr(pm, B, N, POUT); for (i=0; i<N; i++) cout << "pout = " << POUT[i] << endl; cout << endl;
	 for (i=0; i<N; i++) { for (j=0; j<N; j++) printf("%0.8lf\t", Fs[i][j]); cout << endl; }
	 double sumq = 0.0, *T=new double[N], sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	for (i=0; i<N; i++) { if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0); else q[i]=POUT[i];
    QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; } //printf("Summa = %0.10lf\n", sumq); 
	delete []T; delete []q; delete []QA;
	for (i=0; i<N; i++) { p=qm[i]; delete []p; p=em[i]; delete []p; p=pm[i]; delete []p; } 
	 delete []qm; delete []em; delete []pm; delete []B; delete []idf; 
	 return POUT; } //Инвертирование матрицы коэффициентов POUT и умножение на матрицу RHS, чтобы получить POUT
double *SEMIGRAY(int iclsd, double *A, double **EPS, double **RHOs, double **HOs, double ***Fs, int *ID, double *q, double *T, double *PIN, double *POUT, int chst, int chel)
{ int N=chst, i, j, k, *id1=new int[N];
double ikr, sigma=5.670E-8, *epsl=new double[N], *rhosl=new double[N], **FSl=new double*[N], *q1=new double[N], *HOs2=new double[N], *p; // Расчет плотности теплового потока  для внешнего излучения (для кажого диапазона)
for (j=0; j<N; j++) { p=new double[N]; FSl[j]=p; }
    for (i=0; i<N; i++) { id1[i]=1; PIN[i]=0.0; } // Установить значения eps, rhos и Fs для k-ого диапазона
	k=0; for (i=0; i<N; i++) { epsl[i]=EPS[k][i]; rhosl[i]=RHOs[k][i];
        for (j=i; j<N; j++) FSl[i][j]=Fs[k][i][j]; }
	GRAYDIFFSPEC(iclsd, A, epsl, rhosl, HOs[k], FSl, id1, PIN, q1, N); //Вычисление неизвестных плотностей потоков тепловой энергии и температур при q2=q-q1
    k=1; for (i=0; i<N; i++) {
        epsl[i]=EPS[k][i]; rhosl[i]=RHOs[k][i];
        for (j=i; j<N; j++) { 
			FSl[i][j]=Fs[k][i][j];
			FSl[j][i]=A[i]/A[j]*FSl[i][j]; }
        if (iclsd==1) {
			FSl[i][i]=1e0;
            for (j=0; j<N; j++) 
				if (j!=i)
                    FSl[i][i]=FSl[i][i]-(1e0-rhosl[j])*FSl[i][j]; }
        HOs2[i]=-q1[i]/EPS[k][i];
        for (j=0; j<N; j++) HOs2[i]=HOs2[i]+((1e0-rhosl[j])/EPS[2][j]-1e0)*FSl[i][j]*q1[j];
        if (!(ID[i])) PIN[i]=q[i];
		else PIN[i]=sigma*pow(T[i],4e0); }
    GRAYDIFFSPEC(iclsd, A, epsl, rhosl, HOs2, FSl, ID, PIN, POUT, N);
    for (i=0; i<N; i++) {
        if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0);      
        else q[i]=POUT[i]; } 
	delete []id1; delete []epsl; delete []rhosl; delete []q1; delete []HOs2;
	for (j=0; j<N; j++) { p=FSl[j]; delete []p; } delete []FSl;
	return POUT;
}
double VIEW(int NO, int NARG, double *ARG, int fl)
{ double VIEW=0.0, epsi=1e-15, e=epsi, ov=1e0/2e0, A=0.0, B=0.0, C=0.0; int k=0;
if ((NO<38) || (NO>39)) { cout << "Illegal value for number (NO =" << NO << ") " << endl; k=getchar(); exit(1); }
if (NO==38) { if (NARG!=3) {cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1);}
double X=0.0, Y=0.0, RTX=0.0, RTY=0.0, RT=0.0;
      k=0; A=ARG[k]; k++; B=ARG[k]; k++; C=ARG[k];
	  if (fabs(C)>e) { X=A/C; Y=B/C; } else { X=0.0; Y=0.0; }
	  RTX=pow(1.0+X*X, ov); RTY=pow(1.0+Y*Y, ov); RT=pow(1.0+X*X+Y*Y, ov);
      if ((fabs(X)>e) && (fabs(Y)>e))     
	  VIEW=(log(RTX*RTY/RT)+X*RTY*atan(X/RTY)+Y*RTX*atan(Y/RTX)-X*atan(X)-Y*atan(Y))*2.0/(pi*X*Y);
	  else VIEW=0.0; }
else if (NO==39) { if(NARG!=3) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1); }
	double H=0.0, L=0.0, W=0.0, HH=0.0, WW=0.0, W2=0.0, H2=0.0, HW2=0.0, HW=0.0, H12=0.0, W12=0.0, C1=0.0, C2=0.0, C3=0.0;
      k=0; H=ARG[k]; k++; L=ARG[k]; k++; W=ARG[k];
	  if (fabs(L)>e) { HH=H/L; WW=W/L; } else { HH=0.0; WW=0.0; } 
	  W2=WW*WW; H2=HH*HH; HW2=H2+W2; HW=sqrt(H2+W2); H12=H2+1.0; W12=W2+1.0;
	  if (HW>e) { C2=W2*(H12+W2)/W12/HW2; C3=H2*(H12+W2)/H12/HW2; } else { C2=0.0; C3=0.0; }
      C1=W12*H12/(H12+W2);  
	  double atw=0.0, ath=0.0, athw=0.0;
	  if (fabs(WW)>e) atw=atan(1e0/WW); else atw=pi/2e0; 
	  if (fabs(HH)>e) ath=atan(1e0/HH); else ath=pi/2e0; 
	  if (fabs(HW)>e) athw=atan(1e0/HW); else athw=pi/2e0;
      if (fabs(WW)>e) VIEW=(WW*atw+HH*ath-HW*athw+0.25*(log(C1)+W2*log(C2)+H2*log(C3)))/(pi*WW); else VIEW=0.0; //cout << "H = " << HH << "\tW = " << WW << "\tHW = " << HW << "\tpi = " << pi << "\tC1 = " << C1 << "\tC2 = " << C2 << "\tC3 = " << C3 << "atw = " << atw << "\tath = " << ath << "\tathw = " << athw << "\tV = " << VIEW << endl; 
} //if ((fl==13) || (fl==15) || (fl==4)) { cout << "UK = " << VIEW << "\t"; cout << "A = " << A << "\tB = " << B << "\tC = " << C << "\t"; }
return VIEW; }
double perpplates(double X1, double X2, double Y1, double Y2, double Z1, double Z2, double Z3)
{ int NARG=3, k=3, n=39, fl=0; double perppltf=0.0, A=0.0, F=0.0, *ARG=new double[k], epsi=1e-15, e=epsi;
      k=0; ARG[k]=Y2; k++; ARG[k]=Z3; k++; ARG[k]=X2;
      A=X2*Z3;
      k=39; fl++; if (fabs(A*Y2)>e) F=A*VIEW(k, NARG, ARG, fl); //cout << "F1 per = " << F << "\t";
	  //------
	  k=0; ARG[k]=Y1; A=X2*Z3; fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F2 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X1; A=X1*Z3; fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F3 per = " << F << "\t";
	  //------
	  k=0; ARG[k]=Y1; A=X1*Z3; fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F4 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=Z2; A=X1*Z2; fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F5 per = " << F << "\t";
	  //------
      k=0; ARG[k]=Y1; A=X1*Z2; fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F6 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X2; A=X2*Z2; fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F7 per = " << F << "\t";
	  //------
      k=0; ARG[k]=Y1; A=X2*Z2; fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F8 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=(Z3-Z1); A=X2*(Z3-Z1); fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F9 per = " << F << "\t";
	  //------
      k=0; ARG[k]=Y1; A=X2*(Z3-Z1); fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F10 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X1; A=X1*(Z3-Z1); fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F11 per = " << F << "\t";
	  //------
      k=0; ARG[k]=Y1; A=X1*(Z3-Z1); fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F12 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=(Z2-Z1); k++; ARG[k]=X2; A=X2*(Z2-Z1); fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F13 per = " << F << "\t";
	  //------
      k=0; ARG[k]=Y1; A=X2*(Z2-Z1); fl++;
	  if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F14 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X1; A=X1*(Z2-Z1); fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F15 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X1*(Z2-Z1); fl++;
      if(fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F16 per = " << F << "\n";
      if ((fabs(Z1)>e) && (fabs(X2-X1)>e)) perppltf=F/(2.0*(X2-X1)*Z1); else perppltf=0.0;
	  if (ARG) delete []ARG; 
	  return perppltf; }
double parlplates(double X1, double X2, double X3, double Y1, double Y2, double Y3, double C)
{ int NARG=3, j=0, k=0, n=38, fl=0; double parlpltf, A=0.0, F=0.0, *ARG=new double[NARG], epsi=1e-10, e=epsi;
      k=0; ARG[k]=X3; k++; ARG[k]=Y3; k++; ARG[k]=C; k=0; j=1; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=A*VIEW(n, NARG, ARG, fl); //cout << "F1 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F2 par = " << F << "\t";
	  //-----
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F3 par = " << F << "\t";
	  k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F4 par = " << F << "\t";
      //-----
	  k=0; ARG[k]=X2; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
	  if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F5 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F6 par = " << F << "\t";
	  //-----
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F7 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F8 par = " << F << "\t";
	  //-----
	  k=0; ARG[k]=X3-X1; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F9 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F10 par = " << F << "\t";
	  //-----
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F11 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F12 par = " << F << "\t";
	  //-----
	  k=0; ARG[k]=X2-X1; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F13 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F14 par = " << F << "\t";
	  //-----
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F15 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F16 par = " << F << "\n";
      if ((fabs(X1)>e) && (fabs(Y1)>e)) parlpltf=F/(4.0*(X1*Y1));
      if (ARG) delete []ARG; 
	  return parlpltf; }
double *GAUSS(double **A, double *B, int chst, double *X)
{ int N=chst, I, *L=new int[N], K, J, LK;
double *S=new double[N], SMAX, RMAX, R, XMULT, SUM, epsi=1e-10;
if ((!S) || (!L)) { cout << "No memory!" << endl; K=getchar(); exit(1); }
      for (I=0; I<N; I++)
			{ L[I]=I; SMAX=0.0; //L[I] - номер строки
			for (J=0; J<N; J++)
				if (fabs(A[I][J])>SMAX) 
					SMAX=fabs(A[I][J]);
			S[I]=SMAX; } //S[I] - массив максимальных элементов в строке
      for (K=0; K<N-1; K++) 
        { RMAX=0.0; J=K;
        for (I=K; I<N; I++) 
		{ if (fabs(S[L[I]])>epsi) 
				R=fabs(A[L[I]][K])/S[L[I]]; else R=0.0; //делим первые элементы на S[I]
		if (R>RMAX) { J=I; RMAX=R; } } //находим 'элемент с максимальным отношением K-го элемента к максимальному в K-ой строке, записываем его номер
		LK=L[J]; L[J]=L[K]; L[K]=LK; //меняем местами L[K] - текущий и L[J] - с максимальным
		for (I=K+1; I<N; I++)
			{ XMULT = A[L[I]][K]/A[LK][K]; //делим каждый элемент K-ого столбца на максимальный, начиная с (K+1)-го элемента строки 
					if (A[LK][K]==0.0)
						XMULT=0.0;
				for (J=K+1; J<N; J++)
					A[L[I]][J]=A[L[I]][J]-XMULT*A[LK][J]; //вычитаем, чтобы K-ый элемент в K-ом столбце обнулился
				A[L[I]][K] = XMULT; } }
      for (K=0; K<N-1; K++)
        for (I=K+1; I<N; I++)      
          B[L[I]] = B[L[I]] - A[L[I]][K]*B[L[K]];
      if (fabs(A[L[N-1]][N-1])>epsi) 
		  X[N-1] = B[L[N-1]]/A[L[N-1]][N-1]; else X[N-1]=0.0;
      for (I=N-1; I>=0; I--) {
		  SUM = B[L[I]];
		  for (J=I+1; J<N; J++)
          SUM = SUM - A[L[I]][J]*X[J];
        if (fabs(A[L[I]][I])>epsi) 
			X[I] = SUM/A[L[I]][I];
		else X[I]=0.0; }
	  delete []L; delete []S; return X; }
double *MetodGaussa(double **a, double *y, int n, double *xs) 
{ double maxi, temp=0.0, *x=new double[n]; int k=0, index=0, i=0, j=0, q;
if (!x) { cout << "No memory" << endl; i=getchar(); exit(1); } 
for (j=0; j<n; j++) x[j]=xs[j];
k=0; while (k<n) { // Поиск строки с максимальным a[i][k]
    maxi=fabs(a[k][k]); index=k;
    for (i=k+1; i<n; i++) if (fabs(a[i][k])>maxi) { maxi=fabs(a[i][k]); index=i; } 
    if (!(fabs(maxi))) { for (q=0; q<n; q++) x[q]=0.0; return x; //cout << "Resheniye poluchit nevozmozhno iz-za nulevogo stolbtsa " << index << " matritsy A" << endl; k=getchar(); exit(1); 
	} //нет ненулевых диагональных элементов
    for (j=0; j<n; j++) { temp=a[k][j]; a[k][j]=a[index][j]; a[index][j]=temp; } temp=y[k]; y[k]=y[index]; y[index]=temp; //Перестановка строк
    for (i=k; i<n; i++) {
      temp=a[i][k]; if (!(fabs(temp))) continue; // для нулевого коэффициента пропустить
      for (j=0; j<n; j++) if (fabs(temp)>0.0) a[i][j]=a[i][j]/temp; else a[i][j]=0.0;
	  if (fabs(temp)>0.0) y[i]=y[i]/temp; else y[i]=0.0;
      if (i==k) continue; // уравнение не вычитать само из себя
      for (j=0; j<n; j++) a[i][j]=a[i][j]-a[k][j]; y[i]=y[i]-y[k]; } k++; }
for (k=n-1; k>=0; k--) { // обратная подстановка
	x[k]=y[k]; for (i=0; i<k; i++) y[i]=y[i]-a[i][k]*x[k]; } 
return x; }
double SeryeStenkiRasIzlDifPov(double ww, double hh, double ll, double *T, double *EPS, double *HOs, double *RHOs, double *A, int *ns, int chst) //N серых диффузных поверхностей
{ int N=chst, fl=0; double b=0.0, t=0.0, epsi=1e-10; //t=T[ns[0]]; b=T[ns[1]]; if (fabs(t-b)<epsi) return 0.0;
double **Fs=new double*[N], z, r, y; int i=0, j=0, f=3, k=0, d=0; 
double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p=NULL;
for (i=0; i<N; i++) { p=new double[N]; for (j=0; j<N; j++) p[j]=0.0; Fs[i]=p; }
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd=0; for (i=0; i<N; i++) id[i]=1; // Поверхности 1 и 3 (bottom and top)
    	for (i=0; i<N; i++) if (!(id[i])) PIN[i]=q[i]; else PIN[i]=sigma*pow(T[i], 4e0); //Заполнение массива PIN ПТП q и температурой T - перевод температур в ППИ, Fs(1,2)=F1-2, 
	iclsd=1; //угловые коэффициенты; конфигурация замкнута (iclsd=1), диагональные элементы не нужны //for (j=0; j<N; j++) cout << "RHOs ( " << j << " ) = " << RHOs[j] << endl;  for (j=0; j<N; j++) cout << "EPS ( " << j << " ) = " << EPS[j] << endl;
	k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(ww); //h,l,w
	k=39; i=0; j=1; Fs[i][j]=VIEW(k, f, arg, fl); //Fs12 //cout << "Fs(0,1) = " << Fs[0][1] << "\th = " << hh << "\tl = " << ll << "\tw = " << ww << "\tf = " << f << endl; //Fs(1,2)=F1-2, F12
    	k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh);
	k=38; i=0; j=2; Fs[i][j]=VIEW(k, f, arg, fl); //Fs13
    	k=3; d=1; Fs[i][k]=Fs[i][d]; //Fs14
    	k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); 
    k=39; i=0; j=4; Fs[i][j]=VIEW(k, f, arg, fl); //Fs15
	k=5; Fs[i][k]=Fs[i][j]; //Fs16
	   	k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); 
    k=39; i=1; j=2; Fs[i][j]=VIEW(k, f, arg, fl); //Fs23
    	k=0; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ww); 
    k=38; i=1; j=3; Fs[i][j]=VIEW(k, f, arg, fl); //Fs24
    	k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ll);
    k=39; i=1; j=4; Fs[i][j]=VIEW(k, f, arg, fl); //F25 
    i=1; j=5; k=4; Fs[i][j]=Fs[i][k]; //Fs26
	i=2; j=3; k=0; Fs[i][j]=Fs[k][j]; //Fs34
	i=2; j=4; k=0; Fs[i][j]=Fs[k][j]; //Fs35
	i=2; j=5; k=0; Fs[i][j]=Fs[k][j]; //Fs36
	i=3; j=4; k=1; Fs[i][j]=Fs[k][j]; //Fs45
	i=3; j=5; k=1; Fs[i][j]=Fs[k][j]; //Fs46
	k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ll);
    k=38; i=4; j=5; Fs[i][j]=VIEW(k, f, arg, fl); //F56 
	POUT=GRAYDIFF(iclsd, A, EPS, HOs, Fs, id, PIN, POUT, N);
	sumq=0.0; //проверка полной ППИ
    	for (i=0; i<N; i++) //Вывод - перевод в температуры
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; 
			cout << "Surface " << i << "\tT [K] = " << T[i] << "\tq [W/m2] = " << q[i] << "\tQA [W] = " << QA[i]/ll << endl;
		} r=0.0; z=1.0; for (i=0; i<N; i++) if (ns[i]==i) { r=r+z*q[i]; z=-z; } 
		for (i=0; i<N; i++) { p=Fs[i]; delete []p; } delete []Fs; delete []PIN; 
		delete []POUT; delete []arg; delete []q; delete []QA; delete []id; 
		cout << "Summa vsekh teplovykh potokov = " << sumq << endl;
		return r; }