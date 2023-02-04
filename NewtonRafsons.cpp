#define _CRT_SECURE_NO_WARNINGS 
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <ctime>
#include <iostream>
using namespace std;
const int dsnr=50;
void PrintMatr(double **, int);
double *reshMetKram(double **, double *, int);
double Determinant(double **, int);
double **GetMatr(double **, double **, int, int, int);
double vychopred(double  **, int);
double **polmat(double **, double *, int, int);
double **osvpam(double **, int);
double provresh(double **, double *, double *, double *, int, double);
void vyvodmatr(double **, double *, int);
double *MetodGaussa(double **, double *, int, double *);
double *reshMetGau(double **, double *, int, double *);
void napstrNR();
double *reshMetObrMatr(double **, double *, int, double *);
double **ObratMatr(double **, double **, int);
double vychNevyaz(double **, double **, double *, double *, int);
double NormaSystFunct(double *, double *, int, double **, double **, double);
double *rasMetNewtRafs(double *, double **, double **, double *, double, double, int);
char *snmnr;
//------
void napstrNR()
{ int k=0; snmnr=new char[dsnr]; if (!snmnr) { cout << "No memory!" << endl; k=getchar(); exit(1); }
snmnr[k]='N'; k++; snmnr[k]='o'; k++; snmnr[k]='_'; k++; snmnr[k]='m'; k++; 
snmnr[k]='e'; k++; snmnr[k]='m'; k++; snmnr[k]='o'; k++; snmnr[k]='r'; k++; 
snmnr[k]='y'; k++; snmnr[k]='!'; k++; snmnr[k]='\0'; }
void PrintMatr(double **mas, int m) { // Функция вывода матрицы
  int k, j; for (k = 0; k<m; k++) { for (j = 0; j<m; j++) cout << mas[k][j] << "\t"; cout << endl; } }
double **GetMatr(double **mas, double **p, int i, int j, int m) { //Получение матрицы без i-й строки и j-го столбца
  int ki, kj, di, dj;
  di = 0;
  for (ki = 0; ki<(m-1); ki++) { // проверка индекса строки
    if (ki == i) di = 1;
    dj = 0;
    for (kj = 0; kj<(m-1); kj++) { // проверка индекса столбца
      if (kj == j) dj = 1;
      p[ki][kj] = mas[ki + di][kj + dj];
}  } return p; }
double Determinant(double **mas, int m) { // Рекурсивное вычисление определителя
  int i;
  double **p, k=1., d=0.0, *po;
  p=new double*[m]; if (!p) { cout << snmnr << endl; i=getchar(); exit(1); }
  for (i=0; i<m; i++) { po=new double[m]; if (!po) { cout << snmnr << endl; i=getchar(); exit(1); } p[i]=po; }
  d=0;
  if (m<1) cout << "Определитель вычислить невозможно!";
  if (m==1) { d = mas[0][0]; return (d); } 
  if (m==2) { d = mas[0][0]*mas[1][1]-mas[1][0]*mas[0][1]; return(d); } 
  if (m>2) {
	for (i=0; i<m; i++) {
		p=GetMatr(mas, p, i, 0, m); 
		d=d+k*mas[i][0]*Determinant(p, m-1); //разложение по строкам по первому столбцу
		k=-k; } } //(-1) в степени i
for (i=0; i<m; i++) { po=p[i]; delete []po; } return (d); }
double *reshMetKram(double **mat, double *ssc, int rm)
{ int n=rm, j;
double d, *m=new double[n], mm, **mmat; if (!m) { cout << snmnr << endl; j=getchar(); exit(1); } //PrintMatr(mat, rm);
d=vychopred(mat,n); for (j=0; j<n; j++) { mmat=polmat(mat,ssc,j,n); mm=vychopred(mmat,n); 
m[j]=mm/d; /*cout << "x ( " << j << " ) = " << m[j] << endl;*/
mmat=osvpam(mmat,n); } 
return m; }
double **osvpam(double **a, int n)
{ int j; double *b; for (j=0; j<n; j++) { b=a[j]; delete []b; } return NULL; }
double **polmat(double **mat, double *b, int no, int n)
{ int k, j; double **ma=new double*[n], *m; 
if (!ma) { cout << snmnr << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	m=new double[n]; if (!m) { cout << snmnr << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) {
		if (j==no) m[j]=b[k]; else m[j]=mat[k][j];}
ma[k]=m; } 
return (ma); }
double provresh(double **A, double *x, double *xp, double *b, int n, double toch)
{ int k, j, q=0; double s, *pr=new double[n], m, t, mp;
if (!pr) { cout << snmnr << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { s=0.0; 
for (j=0; j<n; j++) s=s+A[k][j]*x[j]; pr[k]=s-b[k]; //cout << "pr = " << pr[k] << "\t"; 
} //cout << endl; 
m=0.0; for (j=0; j<n; j++) { t=fabs(pr[k]); if (t>m) { m=t; q=k; } } 
mp=0.0; for (k=0; k<n; k++) { t=fabs(xp[k]-x[k]); if (t>mp) { mp=t; } }
delete []pr; //cout << "pr = " << m << endl; //return mp; 
return m; } 
double vychopred(double  **mas, int m) {
  double d=Determinant(mas, m); //Вычисление определителя
  return d; }
void vyvodmatr(double **a, double *b, int n) // Вывод системы уравнений
{ int i, j; for (i=0; i<n; i++) { for (j=0; j<n; j++) { cout << a[i][j] << "\t";} 
cout << " --- " << b[i] << endl; } }
double *MetodGaussa(double **a, double *y, int n, double *xs) 
{ double maxi, temp=0.0, *x=new double[n]; int k=0, index=0, i=0, j=0, q;
if (!x) { cout << snmnr << endl; i=getchar(); exit(1); } 
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
double *reshMetGau(double **a, double *y, int n, double *x) 
{ int i, j; double **mt=new double*[n], *b, *xs;
if (!mt) { cout << snmnr << endl; i=getchar(); exit(1); } 
for (j=0; j<n; j++) y[j]=-y[j];
for (i=0; i<n; i++) { b=new double[n]; if (!b) { cout << snmnr << endl; i=getchar(); exit(1); } 
for (j=0; j<n; j++) b[j]=a[i][j]; mt[i]=b; }  
b=new double[n]; if (!b) { cout << snmnr << endl; i=getchar(); exit(1); } 
for (i=0; i<n; i++) b[i]=y[i]; //PrintMatr(a,n); 
xs=MetodGaussa(a,y,n,x); 
for (j=0; j<n; j++) x[j]=xs[j];
for (i=0; i<n; i++) { for (j=0; j<n; j++) a[i][j]=mt[i][j]; y[i]=b[i]; }
delete []b; delete []xs; for (i=0; i<n; i++) { b=mt[i]; delete []b; } //for (i=0; i<n; i++) cout << "x ( " << i << " ) = " << x[i] << endl; //i=getchar();
return x; }
double *reshMetObrMatr(double **fx, double *f, int n, double *dx)
{ int k, j;
double **md2=new double*[n], *u, s, **md3=new double*[n];
if ((!md2) || (!md3)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	u=new double[n]; if (!u) { cout << "No memory!" << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) if (k==j) u[j]=1.0; else u[j]=0.0; md2[k]=u;  }
md3=ObratMatr(fx, md2, n);
for (k=0; k<n; k++) { s=0.0;
for (j=0; j<n; j++) { s=s+md3[k][j]*f[j]; } dx[k]=s; }
for (j=0; j<n; j++) { u=md2[j]; delete []u; u=md3[j]; delete []u; }
delete []md2; delete []md3;
return dx; }
double **ObratMatr(double **mdp, double **mdv, int n)
{ int k, j, p; double perv, vtor, *u, **md1=new double*[n], **md2=new double*[n];
for (j=0; j<n; j++) { u=new double[n]; if (!u) { cout << "No memory!" << endl; k=getchar(); exit(1); } md1[j]=u; 
u=new double[n]; if (!u) { cout << "No memory!" << endl; k=getchar(); exit(1); } md2[j]=u; }
for (k=0; k<n; k++) for (j=0; j<n; j++) { md1[k][j]=mdp[k][j]; md2[k][j]=mdv[k][j]; }
for (k=0; k<n-1; k++) { perv=md1[k][k];
	if (!perv) {}
    else { for (j=0; j<n; j++) { md1[k][j]=md1[k][j]/perv; md2[k][j]=md2[k][j]/perv; } } perv=1.0;
	for (j=k+1; j<n; j++) {
        vtor=md1[j][k];
		if (!vtor) {}
        else { for (p=0; p<n; p++) { md1[j][p]=md1[j][p]*perv/vtor-md1[k][p]; 
		md2[j][p]=md2[j][p]*perv/vtor-md2[k][p]; } } } }
k=n-1; while (k>0) {
    perv=md1[k][k];
	if (!perv) {}
    else { j=n-1; while (j>=0) {
		md1[k][j]=md1[k][j]/perv; md2[k][j]=md2[k][j]/perv; j--; }	}
    j=k-1; while (j>=0) {
        vtor=md1[j][k];
		if (!vtor) {}
        else { p=n-1; while (p>=0) {
            md1[j][p]=md1[j][p]-md1[k][p]*vtor; 
			md2[j][p]=md2[j][p]-md2[k][p]*vtor; p--; } } j--; } k--; }
for (j=0; j<n; j++) { u=md1[j]; delete []u; }
return md2; }
double *rasMetNewtRafs(double *tet, double **ms1, double **ms4, double *ssc, double ksuf, double tora, int ksu)
{ napstrNR(); long q=0, qk=10000, k, j, w; //for (j=0; j<ksu; j++) cout << "te n = ( " << j << " ) = " << tet[j] << "\t"; cout << endl; //for (j=0; j<ksu; j++) cout << "ssc = ( " << j << " ) = " << ssc[j] << "\t"; cout << endl; //cout << "ms1" << endl; for (k=0; k<ksu; k++) { for (j=0; j<ksu; j++) cout << ms1[k][j] << "\t"; cout << endl; } //cout << "ms4" << endl; for (k=0; k<ksu; k++) { for (j=0; j<ksu; j++) cout << ms4[k][j] << "\t"; cout << endl; }
double *fst, *ts1, *ts4, **fdm=new double*[ksu], m=1e2, s, *cp=new double[ksu], *ct=new double[ksu], *fs=new double[ksu], *tett=new double[ksu], t;
if ((!fdm) || (!cp) || (!fs) || (!ct) || (!tett)) { cout << snmnr << endl; k=getchar(); exit(1); }
for (j=0; j<ksu; j++) { cp[j]=0.0; ct[j]=0.0; tett[j]=tet[j]; }
q=0; while ((q<qk) && (m>tora)) {
for (k=0; k<ksu; k++) { ts1=ms1[k]; ts4=ms4[k]; m=0.0; fst=new double[ksu]; if (!fst) { cout << snmnr << endl; j=getchar(); exit(1); }
	for (j=0; j<ksu; j++) { fst[j]=ts1[j]+4.0*ts4[j]*pow(tet[j],3.0); m=m+ts1[j]*tet[j]+ts4[j]*pow(tet[j],4.0); }
fs[k]=m-ssc[k]; fdm[k]=fst; }
ct=reshMetObrMatr(fdm, fs, ksu, ct); 
s=0.0; for (j=0; j<ksu; j++) { tet[j]=tet[j]+ct[j]; s=s+tet[j]; cp[j]=ct[j]; } s=s/ksuf;
m=vychNevyaz(ms1, ms4, tet, ssc, ksu);
q++; } w=q; t=s; //cout << "tsr MOM = " << s << "\tq = " << q << "\t"; //for (j=0; j<ksu; j++) cout << "te k = ( " << j << " ) = " << tet[j] << "\t"; cout << endl; 
//if (w==qk) { 
for (j=0; j<ksu; j++) { fst=fdm[j]; if (fst) delete []fst; tet[j]=tett[j]; cp[j]=0.0; ct[j]=0.0; }
m=1e3; q=0; while ((q<qk) && (m>tora)) {
for (k=0; k<ksu; k++) { ts1=ms1[k]; ts4=ms4[k]; m=0.0; fst=new double[ksu]; if (!fst) { cout << snmnr << endl; j=getchar(); exit(1); }
	for (j=0; j<ksu; j++) { fst[j]=ts1[j]+4.0*ts4[j]*pow(tet[j],3.0); m=m+ts1[j]*tet[j]+ts4[j]*pow(tet[j],4.0); }
fs[k]=m-ssc[k]; fdm[k]=fst; }
ct=reshMetGau(fdm, fs, ksu, ct);
s=0.0; for (j=0; j<ksu; j++) { tet[j]=tet[j]+ct[j]; s=s+tet[j]; cp[j]=ct[j]; } s=s/ksuf;
m=provresh(fdm,tet,ct,fs,ksu,tora); 
q++; } //cout << "tsNR = " << (s+t)/2e0 << "\t"; //cout << "tsr MG = " << s << "\tq = " << q << endl; //} 
for (j=0; j<ksu; j++) { fst=fdm[j]; if (fst) delete []fst; } delete []fdm; delete []cp; delete []ct; delete []fs; delete []tett;
return tet; }
double vychNevyaz(double **ms1, double **ms4, double *tet, double *ssc, int n)
{ int k, j; double *ts1, *ts4, m, t=0.0;
	for (k=0; k<n; k++) { ts1=ms1[k]; ts4=ms4[k]; m=0.0;
	for (j=0; j<n; j++) m=m+ts1[j]*tet[j]+ts4[j]*pow(tet[j],4.0);
	t=t+(m-ssc[k]); } return t; }
double NormaSystFunct(double *tet, double *pdt, int n, double **ms1, double **ms4, double s)
{ double *fst, *ts1, *ts4, m1, m2, m=0.0, t1, m0; int k, j;
	for (k=0; k<n; k++) { ts1=ms1[k]; ts4=ms4[k]; m0=0.0;
	for (j=0; j<n; j++) { 
		t1=tet[j]+s*pdt[j];
		m1=ts1[j]*t1;
		m2=ts4[j]*pow(t1,4.0);
		m0=m0+fabs(m1)+fabs(m2); } m=m+fabs(m0); }
return m; }