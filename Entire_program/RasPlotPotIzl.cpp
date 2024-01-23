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
//-----
double *GRAYDIFFSPEC(int, double *, double *, double *, double *, double **, int *, double *, double *, int, char *);
double *GRAYDIFF(int, double *, double *, double *, double **, int *, double *, double *, int, char *);
double VIEW(int, int, double *, int);
double perpplates(double, double, double, double, double, double, double);
double parlplates(double, double, double, double, double, double, double);
double *MetodGaussa(double **, double *, int, char *, double *);
double opredKTPTKToch(double *, double *, double, int);
double *SeryeStenkiRasIzl(double, double, double, double *, double *, double *, double *, double *, int *, int, char *, int);
double *SeryeStenkiRasIzlCyl(double, double, double *, double *, double *, double *, double *, int *, int, char *, int);
double *SeryeStenkiRasIzlDifPov(double, double, double, double *, double *, double *, double *, int *, int, char *, int);
double F0_lamT(double, double, char *);
double bbfn(double, double);
double opredPokPreSr(double, int, int, double);
char **napolStrok(int, int);
char **napNazvFile(int, int, char *);
double **dliny_voln(char *, char *, char *);
double usredVelichPlank(double *, double *, double *, double, int, double, char *);
double *Kramers_Kronig(double *, int, double *, char *, char *, char *, int, int, int, char *, double);
double **osvPam(double **, int);
//-----
double *SeryeStenkiRasIzl(double ww, double hh, double ll, double *T, double *EPS, double *HOs, double *RHOs, double *A, int *ns, int chst, 
	char *snm, int v) //N серых диффузно излучающих поверхностей с диффузным и зеркальным коэффициентами отражения
{ 
	int N=chst, fl=0; 
	double b=0.0, t=0.0, epsi=1e-20; //t=T[ns[0]]; b=T[ns[1]]; if (fabs(t-b)<epsi) return 0.0;
	double **Fs=new double*[N], z=0.0, r=0.0, y=0.0, hf=1e0; 
	int i=2, j=0, f=3, k=0, u=0, w=0; 
	double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p=NULL, *vm=new double[i];
for (i=0; i<N; i++) { p=new double[N]; if (!p) { cout << snm << endl; k=getchar(); exit(1); }
for (j=0; j<N; j++) p[j]=0.0; Fs[i]=p; }
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd=0; for (i=0; i<N; i++) id[i]=1; // Поверхности 1 и 3 (bottom and top)
if ((!Fs) || (!id) || (!PIN) || (!POUT) || (!arg) || (!vm) || (!q) || (!QA)) 
{ cout << snm << endl; k=getchar(); exit(1); }
    	for (i=0; i<N; i++) if (!(id[i])) PIN[i]=q[i]; else PIN[i]=sigma*pow(T[i], 4e0); //Заполнение массива PIN ПТП q и температурой T - перевод температур в ППИ, Fs(1,2)=F1-2, 
	iclsd=1; //угловые коэффициенты; конфигурация замкнута (iclsd=1), диагональные элементы не нужны //for (j=0; j<N; j++) cout << "RHOs ( " << j << " ) = " << RHOs[j] << endl;  for (j=0; j<N; j++) cout << "EPS ( " << j << " ) = " << EPS[j] << endl;
	k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(ww); //h,l,w
	k=39; i=0; j=1; Fs[i][j]=VIEW(k, f, arg, fl); //Fs12 //cout << "Fs(1,2) = " << Fs[0][1] << "\th = " << hh << "\tl = " << ll << "\tw = " << ww << "\tf = " << f << endl; //Fs(1,2)=F12
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh);     // a,b,c
	k=38; r=VIEW(k, f, arg, fl); 
	z=0.0; t=parlplates(ww, ww, 2.0*ww, ll, z, ll, hh); k=1; b=RHOs[k];
    i=0; j=2; Fs[i][j]=r+b*t; //Fs13 //Fs13=F13+rhos2*F1(2)3	
	t=perpplates(ww, 2.0*ww, z, hh, ll, z, ll); k=1; b=RHOs[k];
    i=0; j=3; k=1; Fs[i][j]=Fs[i][k]+b*t; //Fs14 //cout << "t = " << t << "\tb = " << b << endl; //Fs14=F14+rhos2*F1(2)4
    i=0; j=4; k=1; u=2; w=3; Fs[i][j]=0.5*(hf-(hf-RHOs[k])*Fs[i][k]-Fs[i][u]-Fs[i][w]); //Fs15 //смещены два прямоугольника на ww по оси X
    i=0; j=5; k=4; Fs[i][j]=Fs[i][k]; //Fs16
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); //h,l,w
    k=39; r=VIEW(k, f, arg, fl); 
	k=0; b=RHOs[k];
	t=perpplates(hh, 2.0*hh, z, ww, ll, z, ll);
	i=1; j=2; Fs[i][j]=r+b*t; //Fs23 //Fs23=F23+rhos1*F1(1)3 //смещены два прямоугольника на hh по оси X, т.к. у них разные размеры
	k=0; arg[k]=fabs(ll); k++; arg[k]=fabs(hh); k++; arg[k]=fabs(ww); //a,b,c
    k=38; r=VIEW(k, f, arg, fl);
	k=0; b=RHOs[k];
	t=parlplates(hh, hh, 2.0*hh, ll, z, ll, ww);
	i=1; j=3; Fs[i][j]=r+b*t;  //Fs24 //Fs24=F24+rhos1*F2(1)4
    i=1; j=4; k=0; r=A[k]*Fs[k][i]/A[i];
	i=1; j=4; u=2; w=3; Fs[i][j]=0.5*(hf-(hf-RHOs[k])*r-Fs[i][u]-Fs[i][w]); //Fs25
    i=1; j=5; k=4; Fs[i][j]=Fs[i][k]; //Fs26
    k=0; arg[k]=fabs(ww); k++; arg[k]=fabs(ll); k++; arg[k]=2.0*fabs(hh);  // a,b,c
	k=38; r=VIEW(k, f, arg, fl);
	k=1; b=RHOs[k];
	t=parlplates(ww, ww, 2.0*ww, ll, z, ll, 2.0*hh);
    i=2; j=0; Fs[i][i]=RHOs[j]*(r+b*t); //Fs33 //Fs33=rhos1*F3(1)3+rhos1*rhos2*F3(12+21)3 //излучение стенки саму на себя
	r=perpplates(z, ww, hh, 2.0*hh, ll, z, ll);
	k=1; b=RHOs[k];
	t=perpplates(ww, 2.0*ww, hh, 2.0*hh, ll, z, ll);
    i=2; j=3; k=0; Fs[i][j]=Fs[k][j]+RHOs[k]*(r+b*t); //Fs34 //Fs34=Fs1-4+rhos1*F3(1)-3+rhos1*rhos2*F3(12+21)-4
    i=2; j=4; k=0; r=A[k]*Fs[k][i]/A[i]; 
	u=1; t=A[u]*Fs[u][i]/A[i];
	i=2; j=4; u=1; w=3; Fs[i][j]=0.5*(hf-(hf-RHOs[k])*r-(hf-RHOs[u])*t-Fs[i][i]-Fs[i][w]); //Fs35
    i=2; j=5; k=4; Fs[i][j]=Fs[i][k]; //Fs36
    k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=2.0*fabs(ww); //a,b,c
	k=38; r=VIEW(k, f, arg, fl);
	k=0; b=RHOs[k];
	t=parlplates(hh, hh, 2.0*hh, ll, z, ll, 2.0*ww);
    i=3; j=1; Fs[i][i]=RHOs[j]*(r+b*t); //Fs44 //Fs44=rhos2*F4(2)-4+rhos1*rhos2*F4(12+21)-4
    i=3; j=4; k=0; r=A[k]*Fs[k][i]; 
	u=1; w=2; t=A[u]*Fs[u][i]; b=A[w]*Fs[w][i];
	i=3; j=4; u=1; Fs[i][j]=0.5*(hf-((hf-RHOs[k])*r+(hf-RHOs[u])*t+b)/A[i]-Fs[i][i]); //Fs45
    i=3; j=5; k=4; Fs[i][j]=Fs[i][k]; //Fs46
    i=4; j=5; k=0; r=A[k]*Fs[k][i];
	u=1; w=4; t=A[u]*Fs[u][w]; 
	u=2; k=3; w=4; b=A[u]*Fs[u][w]+A[k]*Fs[k][w];
	u=1; w=4; Fs[i][j]=0.5*(hf-((hf-RHOs[k])*r+(hf-RHOs[u])*t+b)/A[w]); //Fs56 //for (j=0; j<N; j++) cout << "A ( " << j << " ) = " << A[j] << endl; //for (i=0; i<N; i++) { y=RHOs[i]+EPS[i]; cout << "Edinitsa = " << y << endl; } //for (i=0; i<N; i++) cout << "Temp ( " << i << " ) = " << T[i] << endl; //POUT=GRAYDIFFSPEC(iclsd, A, EPS, RHOs, HOs, Fs, id, PIN, POUT, N); //Решение системы уравнений методом Гаусса //for (i=0; i<N; i++) cout << POUT[i] << endl; //for (j=0; j<N; j++) { r=0.0; for (k=0; k<N; k++) r=r+Fs[j][k]; cout << "s = " << r << "\t"; } k=getchar();
	POUT=GRAYDIFFSPEC(iclsd, A, EPS, RHOs, HOs, Fs, id, PIN, POUT, N, snm);
	sumq=0.0; //проверка полной ППИ
    	for (i=0; i<N; i++) //Вывод - перевод в температуры
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; //cout << "Surface " << i << "\tT [K] = " << T[i] << "\tq [W/m2] = " << q[i] << "\tQA [W] = " << QA[i]/ll << endl; 
		} 
		r=0.0; z=hf; j=0; 
		for (i=0; i<N; i++) 
			if (ns[j]==i) 
			{ t=q[i]; vm[j]=t; j++; r=r+z*t; z=-z; } //cout << "Summa vsekh teplovykh potokov = " << sumq << endl;
		for (i=0; i<N; i++) 
		{ p=Fs[i]; if (p) delete[]p; } 
		if (Fs) delete[]Fs; 
		if (PIN) delete []PIN; 
		if (POUT) delete[]POUT; 
		if (arg) delete[]arg; 
		if (QA) delete[]QA; 
		if (id) delete[]id; 
		if (!v) { if (q) delete[]q; return vm; }
		if (v==1) { if (vm) delete[]vm; return q; } 
		if ((v<0) || (v>1)) { cout << "Nepravilno vybran identifikator!"; k=getchar(); exit(1); }
}
double *GRAYDIFF(int iclsd, double *A, double *EPS, double *HOs, double **Fs, int *ID, double *PIN, double *POUT, int chst, char *snm)
{
int N=chst, i=0, j=0, k=0;
double **qm=new double*[N], **em=new double*[N], **pm=new double*[N];
double *B=new double[N], *p=NULL, *idf=new double[N], hf=1e0, s=hf, ikr=hf, epsi=1e-20;
double sumq=0.0, *T=new double[N], sigma=5.67e-8, *q=new double[N], *QA=new double[N], r=0.0;
if ((!qm) || (!em) || (!pm) || (!B) || (!idf) || (!T) || (!q) || (!QA)) 
{ cout << snm << endl; i=getchar(); exit(1); }
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; if (!p) { cout << snm << endl; i=getchar(); exit(1); }
p=new double[N]; em[i]=p; if (!p) { cout << snm << endl; i=getchar(); exit(1); } 
p=new double[N]; pm[i]=p; if (!p) { cout << snm << endl; i=getchar(); exit(1); } }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<i; j++)
		  if (fabs(A[i])>epsi) 
			  Fs[i][j]=(A[j]/A[i])*Fs[j][i]; else Fs[i][j]=0.0; //Для закрытой конфигурации, необходимо рассчитать диагональные элементы из правила суммирования
	if (iclsd==1) { 
		for (i=0; i<N; i++) {
		Fs[i][i]=hf;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-Fs[i][j]; } } } }
	 for (i=0; i<N; i++) { // Заполнение матриц коэффициентов q и e
		for (j=0; j<N; j++) {
			if (i==j) ikr=hf; else ikr=0.0; //Символ Кронекера delta_ij 
			if (fabs(EPS[j])>epsi) 
				qm[i][j]=ikr/EPS[j]-(hf/EPS[j]-hf)*Fs[i][j];
			em[i][j]=ikr-Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Заполнение матрицы выходных коэффициентов POUT и RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(hf-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(hf-idf[j]))*PIN[j]; } }
	POUT=MetodGaussa(pm, B, N, snm, POUT); 
	for (i=0; i<N; i++) { if (!(ID[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); else q[i]=POUT[i];
    QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; } 
	if (T) delete []T; if (q) delete []q; if (QA) delete []QA; 
	qm=osvPam(qm, N);
	pm=osvPam(pm, N);
	em=osvPam(em, N);
	if (B) delete []B; if (idf) delete []idf; 
	return POUT;
}
double *GRAYDIFFSPEC(int iclsd, double *A, double *EPS, double *RHOs, double *HOs, double **Fs, int *ID, double *PIN, double *POUT, int chst, 
	char *snm)
{ 
	int N=chst, i=0, j=0; 
	double **qm=new double*[N], **em=new double*[N], **pm=new double*[N], *B=new double[N], *p=NULL;
	double *idf=new double[N], hf=1e0, s=0.0, ikr=s, y=s, epsi=1e-20; // Вычисление недостающих угловых коэффициентов - левой нижней части матрицы из свойства взаимности
	if ((!qm) || (!em) || (!pm) || (!B) || (!idf)) { cout << snm << endl; i=getchar(); exit(1); }
for (i=0; i<N; i++) { p=new double[N]; qm[i]=p; if (!p) { cout << snm << endl; i=getchar(); exit(1); } 
p=new double[N]; em[i]=p; if (!p) { cout << snm << endl; i=getchar(); exit(1); } 
p=new double[N]; pm[i]=p; if (!p) { cout << snm << endl; i=getchar(); exit(1); } }
for (i=0; i<N; i++) { s=0.0; for (j=0; j<ID[i]; j++) s=s+hf; idf[i]=s; }
for (i=1; i<N; i++)
	 for (j=0; j<i; j++)
		  if (fabs(A[i])>epsi) 
			  Fs[i][j]=(A[j]/A[i])*Fs[j][i]; else Fs[i][j]=0.0; //Для закрытой конфигурации, необходимо рассчитать диагональные элементы из правила суммирования
	if (iclsd==1) { for (i=0; i<N; i++) {
		Fs[i][i]=hf;
		for (j=0; j<N; j++) {
			if (j!=i) {
				Fs[i][i]=Fs[i][i]-(hf-RHOs[j])*Fs[i][j]; } }
		if (RHOs[i]<hf) Fs[i][i]=Fs[i][i]/(hf-RHOs[i]); else RHOs[i]=0.0; } }
	 for (i=0; i<N; i++) { // Заполнение матриц коэффициентов q и e
		for (j=0; j<N; j++) {
			if (i==j) ikr=hf; else ikr=0.0; // Символ Кронекера delta_ij 
			if (fabs(EPS[j])>epsi) 
				qm[i][j]=ikr/EPS[j]-((hf-RHOs[j])/EPS[j]-hf)*Fs[i][j];
			em[i][j]=ikr-(hf-RHOs[j])*Fs[i][j]; } }
	 for (i=0; i<N; i++) { // Заполнение матрицы выходных коэффициентов POUT и RHS
			B[i]=-HOs[i];
			for (j=0; j<N; j++) {
					pm[i][j]=qm[i][j]*idf[j]-em[i][j]*(hf-idf[j]);
					B[i]=B[i]+(em[i][j]*idf[j]-qm[i][j]*(hf-idf[j]))*PIN[j]; } } 
	 POUT=MetodGaussa(pm, B, N, snm, POUT); 
	 double sumq=0.0, *T=new double[N], sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	for (i=0; i<N; i++) { if (!(ID[i])) T[i]=pow((POUT[i]/sigma),1e0/4e0); else q[i]=POUT[i];
    QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; } 
	if (T) delete[]T; if (q) delete[]q; if (QA) delete[]QA;
	qm=osvPam(qm, N);
	pm=osvPam(pm, N);
	em=osvPam(em, N);
	if (qm) delete[]qm; if (em) delete []em; if (pm) delete []pm; if (B) delete[]B; if (idf) delete[]idf; 
return POUT; } //Инвертирование матрицы коэффициентов POUT и умножение на матрицу RHS, чтобы получить POUT
double VIEW(int NO, int NARG, double *ARG, int fl)
{ 
	int k=0;
	double VIEWF=0.0, epsi=1e-20, e=epsi, hf=1e0, ov=hf/2e0, A=0.0, B=0.0, C=0.0; 
	double RR1=0.0, RR2=0.0, H=0.0;
	double X=0.0, Y=0.0, RTX=0.0, RTY=0.0, RT=0.0;
	double L=0.0, W=0.0, HH=0.0, WW=0.0, W2=0.0, H2=0.0, HW2=0.0, HW=0.0;
	double H12=0.0, W12=0.0, C1=0.0, C2=0.0, C3=0.0;
	double atw=0.0, ath=0.0, athw=0.0, R1=0.0, R2=0.0, R=0.0;
if ((NO<38) || (NO>43)) { cout << "Illegal value for number (NO =" << NO << ") " << endl; k=getchar(); exit(1); }
if (NO==38) { //два противоположных основания прямоугольного параллелепипеда
	if (NARG!=3) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1);}
      k=0; A=ARG[k]; k++; B=ARG[k]; k++; C=ARG[k];
	  if (fabs(C)>e) { X=A/C; Y=B/C; } else { X=0.0; Y=0.0; }
	  RTX=pow(hf+X*X, ov); RTY=pow(hf+Y*Y, ov); RT=pow(hf+X*X+Y*Y, ov);
      if ((fabs(X)>e) && (fabs(Y)>e))     
	  VIEWF=(log(RTX*RTY/RT)+X*RTY*atan(X/RTY)+Y*RTX*atan(Y/RTX)-X*atan(X)-Y*atan(Y))*2.0/(pi*X*Y);
	  else VIEWF=0.0; }
else if (NO==39) { //две перпендикулярные (смежные) грани
	if (NARG!=3) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1); }
      k=0; H=ARG[k]; k++; L=ARG[k]; k++; W=ARG[k];
	  if (fabs(L)>e) { HH=H/L; WW=W/L; } else { HH=0.0; WW=0.0; } 
	  W2=WW*WW; H2=HH*HH; HW2=H2+W2; HW=sqrt(H2+W2); H12=H2+1.0; W12=W2+1.0;
	  if (HW>e) { C2=W2*(H12+W2)/W12/HW2; C3=H2*(H12+W2)/H12/HW2; } else { C2=0.0; C3=0.0; }
      C1=W12*H12/(H12+W2);
	  if (fabs(WW)>e) atw=atan(1e0/WW); else atw=pi/2e0; 
	  if (fabs(HH)>e) ath=atan(1e0/HH); else ath=pi/2e0; 
	  if (fabs(HW)>e) athw=atan(1e0/HW); else athw=pi/2e0;
      if (fabs(WW)>e) VIEWF=(WW*atw+HH*ath-HW*athw+0.25*(log(C1)+W2*log(C2)+H2*log(C3)))/(pi*WW); else VIEWF=0.0; //cout << "H = " << HH << "\tW = " << WW << "\tHW = " << HW << "\tpi = " << pi << "\tC1 = " << C1 << "\tC2 = " << C2 << "\tC3 = " << C3 << "atw = " << atw << "\tath = " << ath << "\tathw = " << athw << "\tV = " << VIEWF << endl; 
} 
else if (NO==40) { //от одного основания цилиндра на другое
	if (NARG!=3) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1); }
		k=0; A=ARG[k]; k++; R1=ARG[k]; k++; R2=ARG[k];
		if (fabs(A)>e) { RR1=R1/A; RR2=R2/A; }
        X=hf+(hf+RR2*RR2)/(RR1*RR1);
		B=pow((R2/R1), 2.0);
		B=pow(X*X-4.0*B, ov);
        VIEWF=ov*(X-B); }
else if (NO==41) { //от внутренней поверхности цилиндра на основание
	if (NARG!=2) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1); }
	k=0; H=ARG[k]; k++; R=ARG[k]; HH=ov*H/R;
	A=hf+pow(HH, 2e0);
	A=pow(A, ov);
	VIEWF=ov*(A-HH);
}
else if (NO==42) { //внутренняя поверхность сама на себя
	if (NARG!=2) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1); }
		k=0; H=ARG[k]; k++; R=ARG[k]; HH=ov*H/R;
		A=hf+pow(HH, 2e0);
		A=pow(A, ov);
        VIEWF=hf+HH-A;
}
else if (NO==43) {
	if (NARG!=2) { cout << "Wrong number of input parameters (NARG =" << NARG << ") for NO = " << NO << endl; k=getchar(); exit(1); }
		k=0; H=ARG[k]; k++; R=ARG[k];
        HH=ov*H/R;
		A=hf+pow(HH, 2e0);
		A=pow(A, ov);
        VIEWF=2.0*HH*(A-HH);
}
return VIEWF; }
double perpplates(double X1, double X2, double Y1, double Y2, double Z1, double Z2, double Z3)
{ int NARG=3, k=3, n=39, fl=0; double perppltf=0.0, A=0.0, F=0.0, *ARG=new double[k], epsi=1e-20, e=epsi;
      k=0; ARG[k]=Y2; k++; ARG[k]=Z3; k++; ARG[k]=X2;
      A=X2*Z3;
      k=39; fl++; if (fabs(A*Y2)>e) F=A*VIEW(k, NARG, ARG, fl); //cout << "F1 per = " << F << "\t";
	  k=0; ARG[k]=Y1; A=X2*Z3; fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F2 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X1; A=X1*Z3; fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F3 per = " << F << "\t";
	  k=0; ARG[k]=Y1; A=X1*Z3; fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F4 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=Z2; A=X1*Z2; fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F5 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X1*Z2; fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F6 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X2; A=X2*Z2; fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F7 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X2*Z2; fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F8 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=(Z3-Z1); A=X2*(Z3-Z1); fl++;
      if (fabs(A*Y2)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F9 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X2*(Z3-Z1); fl++;
      if (fabs(A*Y1)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F10 per = " << F << "\t";
      k=0; ARG[k]=Y2; k=2; ARG[k]=X1; A=X1*(Z3-Z1); fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F11 per = " << F << "\t";
      k=0; ARG[k]=Y1; A=X1*(Z3-Z1); fl++;
      if (fabs(A*Y1)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F12 per = " << F << "\t";
      k=0; ARG[k]=Y2; k++; ARG[k]=(Z2-Z1); k++; ARG[k]=X2; A=X2*(Z2-Z1); fl++;
      if (fabs(A*Y2)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F13 per = " << F << "\t";
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
{ 
	int NARG=3, j=0, k=0, n=38, fl=0; 
	double parlpltf=0.0, A=0.0, F=0.0, *ARG=new double[NARG], epsi=1e-20, e=epsi;
      k=0; ARG[k]=X3; k++; ARG[k]=Y3; k++; ARG[k]=C; k=0; j=1; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=A*VIEW(n, NARG, ARG, fl); //cout << "F1 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F2 par = " << F << "\t";
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F3 par = " << F << "\t";
	  k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F4 par = " << F << "\t";
	  k=0; ARG[k]=X2; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
	  if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F5 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F6 par = " << F << "\t";
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F7 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F8 par = " << F << "\t";
	  k=0; ARG[k]=X3-X1; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F9 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F10 par = " << F << "\t";
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F11 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F12 par = " << F << "\t";
	  k=0; ARG[k]=X2-X1; j=1; ARG[j]=Y3; A=ARG[k]*ARG[j]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F13 par = " << F << "\t";
      k=1; ARG[k]=Y2; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F14 par = " << F << "\t";
	  k=1; ARG[k]=Y3-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F-A*VIEW(n, NARG, ARG, fl); //cout << "F15 par = " << F << "\t";
      k=1; ARG[k]=Y2-Y1; j=0; A=ARG[j]*ARG[k]; fl++;
      if (fabs(A)>e) F=F+A*VIEW(n, NARG, ARG, fl); //cout << "F16 par = " << F << "\n";
      if ((fabs(X1)>e) && (fabs(Y1)>e)) parlpltf=F/(4.0*(X1*Y1));
      if (ARG) delete []ARG; 
	  return parlpltf; }
double *SeryeStenkiRasIzlDifPov(double ww, double hh, double ll, double *T, double *EPS, double *HOs, double *A, int *ns, int chst, 
	char *snm, int v) //N серых диффузных поверхностей
{ 
	int N=chst, fl=0; double b=0.0, t=0.0, epsi=1e-20; //t=T[ns[0]]; b=T[ns[1]]; if (fabs(t-b)<epsi) return 0.0;
	double **Fs=new double*[N], z=0.0, r=0.0, y=0.0, hf=1e0; int i=2, j=0, f=3, k=0, d=0; //f - число аргументов
	double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p=NULL, *vm=new double[i];
for (i=0; i<N; i++) { 
	p=new double[N]; 
	if (!p) { cout << snm << endl; k=getchar(); exit(1); } 
	for (j=0; j<N; j++) p[j]=0.0; 
	Fs[i]=p; }
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd=0; 
	for (i=0; i<N; i++) id[i]=1; // Поверхности 1 и 3 (bottom and top)
	if ((!Fs) || (!id) || (!PIN) || (!POUT) || (!arg) || (!vm) || (!q) || (!QA)) 
	{ cout << snm << endl; k=getchar(); exit(1); }
    	for (i=0; i<N; i++) 
			if (!(id[i])) PIN[i]=q[i]; 
			else PIN[i]=sigma*pow(T[i], 4e0); //Заполнение массива PIN ПТП q и температурой T - перевод температур в ППИ, Fs(1,2)=F1-2, 
	iclsd=1; //угловые коэффициенты; конфигурация замкнута (iclsd=1), диагональные элементы не нужны 
	k=0; arg[k]=fabs(hh); k++; arg[k]=fabs(ll); k++; arg[k]=fabs(ww); //h,l,w
	k=39; i=0; j=1; Fs[i][j]=VIEW(k, f, arg, fl); //Fs12 
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
	POUT=GRAYDIFF(iclsd, A, EPS, HOs, Fs, id, PIN, POUT, N, snm); 
    	for (i=0; i<N; i++) //Вывод - перевод в температуры
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); 
		else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; 
		} 
		r=0.0; z=hf; j=0; 
		for (i=0; i<N; i++) 
			{ if (ns[j]==i) 
				{ t=q[i]; vm[j]=t; j++; r=r+z*t; z=-z; } } 
		Fs=osvPam(Fs, N);
		if (PIN) delete[]PIN; 
		if (POUT) delete[]POUT; 
		if (arg) delete[]arg; 
		if (QA) delete[]QA; 
		if (id) delete[]id; 
		if (!v) { if (q) delete[]q; return vm; }
		if (v==1) { if (vm) delete[]vm; return q; } 
		if ((v<0) || (v>1)) { cout << "Nepravilno vybran identifikator!" << endl; k=getchar(); exit(1); }
}
double F0_lamT(double T0, double lam, char *snm)
{ 
	int k=20, j=0; 
	double *F0=new double[k], *lamT=new double[k], dolya=0.0, hf=1e0, hmi=1e-3, ko=1e-6, lT0=0.0;
if ((!F0) || (!lamT)) { cout << snm << endl; j=getchar(); exit(1); }
if (lam>hmi) lam=lam*ko; lT0=lam*T0;
lamT[j]=0.18889e-2; F0[j]=0.05059; j++; lamT[j]=0.22222e-2; F0[j]=0.10503; j++; lamT[j]=0.24444e-2; F0[j]=0.14953; j++;
lamT[j]=0.27222e-2; F0[j]=0.21033; j++; lamT[j]=0.28889e-2; F0[j]=0.24803; j++; lamT[j]=0.31667e-2; F0[j]=0.31067; j++;
lamT[j]=0.33889e-2; F0[j]=0.35933; j++; lamT[j]=0.36111e-2; F0[j]=0.40585; j++; lamT[j]=0.38889e-2; F0[j]=0.46031; j++;
lamT[j]=0.41111e-2; F0[j]=0.50066; j++; lamT[j]=0.44444e-2; F0[j]=0.55573; j++; lamT[j]=0.47778e-2; F0[j]=0.60449; j++;
lamT[j]=0.51667e-2; F0[j]=0.65402; j++; lamT[j]=0.56111e-2; F0[j]=0.70211; j++; lamT[j]=0.61667e-2; F0[j]=0.75146; j++;
lamT[j]=0.68889e-2; F0[j]=0.80152; j++; lamT[j]=0.78889e-2; F0[j]=0.85171; j++; lamT[j]=0.93889e-2; F0[j]=0.90031; j++;
lamT[j]=1.25556e-2; F0[j]=0.95094; j++; lamT[j]=2.33333e-2; F0[j]=0.99051; j++;
dolya=opredKTPTKToch(F0, lamT, lT0, k); if (dolya<0.0) dolya=0.0; if (dolya>=hf) dolya=hf; 
if (F0) delete[]F0; if (lamT) delete[]lamT; 
return dolya; }
double bbfn(double hp, double T)
{   
	double hmi=1e-3, ko=1e-6; //if (hp>hmi) hp=hp*ko;
	double X=T*hp, CC=1.5e1/pow(pi,4e0), C2=1.4388e4;
	double EPS=1e-16, V=C2/X, EX=exp(V), M=0.0, BBFN=0.0, VM, BM, EM=1e0, hf=EM; 
	int j=0, jk=1000;
do {
    M=M+hf;
    VM=M*V;
    BM=(6e0+VM*(6e0+VM*(3e0+VM)))/pow(M, 4e0);
    EM=EM/EX;
    BBFN=BBFN+BM*EM; j++; }
while (((pow(VM,3e0)*EM)>EPS) && (j<jk));
    BBFN=CC*BBFN;
	if (BBFN<0.0) BBFN=0.0; if (BBFN>=hf) BBFN=hf; 
    return BBFN; }
double opredPokPreSr(double tem, int vybves, int vybmar, double dko)
{ 
	int k=0, j=0, q=j+1, l=q+1, i=0, dlar=0; 
	double hf=1e0, ume=hf;
	char **mauk=napolStrok(vybves, vybmar);
	char *snm=mauk[j], *sfno=mauk[q]; 
	char **mu=napNazvFile(vybves, vybmar, snm);
	char *sppov=mu[j], *skpov=mu[q], *sdvov=mu[l]; 
	double **unau=dliny_voln(snm, sdvov, sfno);
	double *dv=unau[j], *ne=unau[q], e=1e-9, d=e, t=ne[j]; 
	while (d<t) { d=d+hf; i++; } dlar=i; if (ne) delete[]ne; if (unau) delete[]unau;
double *npp=new double[dlar]; if (!npp) { cout << snm << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0;
npp=Kramers_Kronig(dv, dlar, npp, sppov, sfno, snm, vybves, vybmar, dlar, skpov, dko);
t=usredVelichPlank(dv, npp, npp, tem, dlar, ume, snm);
for (k=0; k<=q; k++) { snm=mauk[k]; if (snm) delete[]snm; } if (mauk) delete[]mauk; 
for (k=0; k<=l; k++) { snm=mu[k]; if (snm) delete[]snm; } if (mu) delete[]mu; 
if (dv) delete []dv; if (npp) delete []npp;
return t; }
double *SeryeStenkiRasIzlCyl(double r, double h, double *T, double *EPS, double *HOs, double *RHOs, double *A, int *ns, int chst, 
	char *snm, int v) //N серых диффузно излучающих поверхностей с диффузным и зеркальным коэффициентами отражения
{ 
	int N=chst, fl=0; 
	double b=0.0, t=0.0, epsi=1e-20; 
	double **Fs=new double*[N], z=0.0, y=0.0, hf=1e0, g=0.0; 
	int i=2, j=0, f=3, k=0, u=0, w=0; 
	double *PIN=new double[N], *POUT=new double[N], *arg=new double[f], *p=NULL, *vm=new double[i];
	double sumq=0.0, sigma=5.67e-8, *q=new double[N], *QA=new double[N];
	int *id=new int[N], iclsd=0; 
	if ((!Fs) || (!id) || (!PIN) || (!POUT) || (!arg) || (!vm) || (!q) || (!QA)) 
	{ cout << snm << endl; k=getchar(); exit(1); }
for (i=0; i<N; i++) { p=new double[N]; if (!p) { cout << snm << endl; k=getchar(); exit(1); }
for (j=0; j<N; j++) p[j]=0.0; Fs[i]=p; }
	for (i=0; i<N; i++) id[i]=1; // Поверхности 1 и 3 (bottom and top)
    	for (i=0; i<N; i++) if (!(id[i])) PIN[i]=q[i]; else PIN[i]=sigma*pow(T[i], 4e0); //Заполнение массива PIN ПТП q и температурой T - перевод температур в ППИ
	iclsd=1; //угловые коэффициенты; конфигурация замкнута (iclsd=1)
	k=0; arg[k]=fabs(h); k++; arg[k]=fabs(r); k++; arg[k]=fabs(r); //r, h
	k=43; i=0; j=1; w=2; Fs[i][j]=VIEW(k, w, arg, fl); //F12 //основание на внутреннюю поверхность
    k=40; i=0; j=2; w++; Fs[i][j]=VIEW(k, w, arg, fl); //F13 //от одного основания на другое
    k=42; i=1; j=i; w--; Fs[i][j]=VIEW(k, w, arg, fl); //Fs22 //внутренняя поверхность на себя
	k=41; i=1; j=2; Fs[i][j]=VIEW(k, w, arg, fl); //Fs23 //внутренняя поверхность на основание
	POUT=GRAYDIFF(iclsd, A, EPS, HOs, Fs, id, PIN, POUT, N, snm);
	sumq=0.0; //проверка полной ППИ
    	for (i=0; i<N; i++) //Вывод - перевод в температуры
		{ if (!(id[i])) T[i]=pow((POUT[i]/sigma), 1e0/4e0); else q[i]=POUT[i];
        	QA[i]=q[i]*A[i]; sumq=sumq+QA[i]; 
		} 
		g=0.0; z=hf; j=0; 
		for (i=0; i<N; i++) 
			if (ns[j]==i) 
			{ t=q[i]; vm[j]=t; j++; g=g+z*t; z=-z; } 
		Fs=osvPam(Fs, N);
		if (PIN) delete []PIN; 
		if (POUT) delete[]POUT; 
		if (arg) delete[]arg; 
		if (QA) delete[]QA; 
		if (id) delete[]id; 
		if (!v) { if (q) delete[]q; return vm; }
		if (v==1) { if (vm) delete[]vm; return q; } 
		if ((v<0) || (v>1)) { cout << "Nepravilno vybran identifikator!" << endl; k=getchar(); exit(1); }
}