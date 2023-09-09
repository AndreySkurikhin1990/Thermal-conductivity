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
void PrintMatr(double **, int);
double *reshMetKram(double **, double *, int, char *);
double Determinant(double **, int, char *);
double **GetMatr(double **, double **, int, int, int);
double **polmat(double **, double *, int, int, char *);
double **osvpam(double **, int);
double provresh(double **, double *, double *, double *, int, double);
double *MetodGaussa(double **, double *, int, char *);
double *reshMetObrMatr(double **, double *, int, char *);
double **ObratMatr(double **, double **, int, char *);
double vychNevyaz(double **, double **, double *, double *, int);
double NormaSystFunct(double *, double *, int, double **, double **, double);
double *opredTempStenFragmMasTem(int, double *, double *, double *, double *, int, int, double, double, double, double, double, char *, int, int);
double *reshnewrafsokon(double *, double **, double **, double *,  int, double, double *, int, double *, int, double, double, double, 
	double, double, double *, double *, double *, int, char *, double);
double *reshnewtrafs(double ***, double ***, double *, double *, double *, double *, int, double *, double *, int, double *, double *, int, double, 
	double, double *, double, double *, int, double, double, double, char *, double, double *, double, double, int, double, double *, 
	double, double, double, double, double, int, int, double *, int, double *, int, int, int);
double opredKTPTKToch(double *, double *, double, int);
double *usrednen(double *, double *, int, int);
double *koefPrib(double *, double *, int, char *);
double *GAUSS(double **, double *, int, double *);
void zapisvfile(double *, int, char *);
double *koefPribl(double *, double *, int, char *);
double *RasIzlSerStenNacSok(double *, double *, double *, double ***, double ***, double, double, double, double, double, double, double, int, 
	double *, double *, int, double *, double *, char *, double *, int, double, double, double, double, double, double, double, int, int, int,
	double *, int, double, int, int, int, int);
double *oprProizPTPTepl(double *, double *, int, char *);
double **poiskMatrYakobi(double *, int, double *, double *, double *, double ***, double ***, double, double, int, double, double, double, int, double *, 
	double *, int, double *, double *, char *, double *, int, double, double, double, double, double, double, int, int, double *, int, 
	double *, double *, int, int, int);
double *opredMatrF(double *, int, double *, double *, double *, double ***, double ***, double, double, int, double, double, double, int, double *, 
	double *, int, double *, double *, char *, double *, int, double, double, double, double, double, double, int, int, double *, int, 
	double *, double *, int, int, int);
double opredTempStenFragm(double, int, double *, double *, double *, double *, int, int, double, double, double, double);
int provDeterNul(double **, int, char *);
double *perPlotTeplPot(double *, double *, int, char *, double ***, int, int);
//------
void PrintMatr(double **mas, int m) { // Функция вывода матрицы
  int k=0, j=0; for (k = 0; k<m; k++) { for (j = 0; j<m; j++) cout << "A ( " << k << " , " << j << " ) = " << mas[k][j] << "\t"; cout << endl; } }
double **GetMatr(double **mas, double **p, int i, int j, int m) 
{ //Получение матрицы без i-й строки и j-го столбца
  int ki=0, kj=0, di=0, dj=0;
  di = 0;
  for (ki = 0; ki<(m-1); ki++) { // проверка индекса строки
    if (ki == i) di = 1;
    dj = 0;
    for (kj = 0; kj<(m-1); kj++) { // проверка индекса столбца
      if (kj == j) dj = 1;
      p[ki][kj] = mas[ki + di][kj + dj];
} } 
  return p; }
double Determinant(double **mas, int m, char *snmnr) { // Рекурсивное вычисление определителя
  int i=0, w=0, v=1;
  double **p=NULL, k=1e0, d=0.0, *po=NULL;
  p=new double*[m]; if (!p) { cout << snmnr << endl; i=getchar(); exit(1); }
  for (i=0; i<m; i++) { po=new double[m]; if (!po) { cout << snmnr << endl; i=getchar(); exit(1); } p[i]=po; }
  d=0;
  if (m<1) cout << "Определитель вычислить невозможно!";
  if (m==1) { d = mas[w][w]; return (d); } 
  if (m==2) { d = mas[w][w]*mas[v][v]-mas[v][w]*mas[w][v]; 
  return(d); } 
  if (m>2) {
	for (i=0; i<m; i++) {
		p=GetMatr(mas, p, i, 0, m); 
		d=d+k*mas[i][w]*Determinant(p, m-1, snmnr); //разложение по строкам по первому столбцу
		k=-k; } } //(-1) в степени i
for (i=0; i<m; i++) { po=p[i]; if (po) delete []po; } 
return (d); }
double *reshMetKram(double **mat, double *ssc, int rm, char *snmnr)
{ 
	int n=rm, j=0;
	double d=0.0, *m=new double[n], mm=d, **mmat=NULL; 
if (!m) { cout << snmnr << endl; j=getchar(); exit(1); } //PrintMatr(mat, rm);
d=Determinant(mat, n, snmnr); 
for (j=0; j<n; j++) { mmat=polmat(mat, ssc, j, n, snmnr); mm=Determinant(mmat, n, snmnr); 
if (fabs(d)>0.0) m[j]=mm/d; 
mmat=osvpam(mmat, n); } 
return m; }
double **osvpam(double **a, int n)
{ int j=0; double *b=NULL; for (j=0; j<n; j++) { b=a[j]; if (b) delete []b; } return NULL; }
double **polmat(double **mat, double *b, int no, int n, char *snmnr)
{ int k=0, j=0; double **ma=new double*[n], *m=NULL; 
if (!ma) { cout << snmnr << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	m=new double[n]; if (!m) { cout << snmnr << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) {
		if (j==no) m[j]=b[k]; else m[j]=mat[k][j];}
ma[k]=m; } 
return (ma); }
double provresh(double **A, double *x, double *xp, double *b, int n, double toch, char *snmnr)
{ 
	int k=0, j=0, q=0; 
	double s=0.0, *pr=new double[n], m=0.0, t=0.0, mp=0.0;
if (!pr) { cout << snmnr << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { s=0.0; 
for (j=0; j<n; j++) s=s+A[k][j]*x[j]; pr[k]=s-b[k]; //cout << "pr = " << pr[k] << "\t"; 
} //cout << endl; 
m=0.0; for (j=0; j<n; j++) { t=fabs(pr[k]); if (t>m) { m=t; q=k; } } 
mp=0.0; for (k=0; k<n; k++) { t=fabs(xp[k]-x[k]); if (t>mp) { mp=t; } }
if (pr) delete[]pr; //cout << "pr = " << m << endl; //return mp; 
return m; } 
double *MetodGaussa(double **a, double *y, int n, char *snm) 
{
	double maxi=0.0, temp=0.0, *x=new double[n], e=1e-15;
	int k=0, index=0, i=0, j=0, q=0;
if (!x) { cout << snm << endl; i=getchar(); exit(1); } 
k=0; while (k<n) { // Поиск строки с максимальным a[i][k]
    maxi=fabs(a[k][k]); index=k;
    for (i=k+1; i<n; i++) if (fabs(a[i][k])>maxi) { maxi=fabs(a[i][k]); index=i; } 
    if (fabs(maxi)<e) { for (q=0; q<n; q++) x[q]=0.0; return x; //cout << "Resheniye poluchit nevozmozhno iz-za nulevogo stolbtsa " << index << " matritsy A" << endl; k=getchar(); exit(1); 
	} //нет ненулевых диагональных элементов
    for (j=0; j<n; j++) { temp=a[k][j]; a[k][j]=a[index][j]; a[index][j]=temp; } 
	temp=y[k]; y[k]=y[index]; y[index]=temp; //Перестановка строк
    for (i=k; i<n; i++) {
      temp=a[i][k]; if (fabs(temp)<e) continue; // для нулевого коэффициента пропустить
      for (j=0; j<n; j++) if (fabs(temp)>0.0) a[i][j]=a[i][j]/temp; else a[i][j]=0.0;
	  if (fabs(temp)>e) y[i]=y[i]/temp; else y[i]=0.0;
      if (i==k) continue; // уравнение не вычитать само из себя
      for (j=0; j<n; j++) a[i][j]=a[i][j]-a[k][j]; y[i]=y[i]-y[k]; } k++; }
for (k=n-1; k>=0; k--) { // обратная подстановка
	x[k]=y[k]; for (i=0; i<k; i++) y[i]=y[i]-a[i][k]*x[k]; } 
return x; }
double *reshMetObrMatr(double **fx, double *f, int n, char *snmnr)
{ 
	int k=0, j=0;
	double **md2=new double*[n], *u=NULL, s=0.0, **md3=new double*[n], hf=1e0, *dx=new double[n];
if ((!md2) || (!md3) || (!dx)) { cout << snmnr << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { //создаем единичную матрицу
	u=new double[n]; if (!u) { cout << snmnr << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) if (k==j) u[j]=hf; else u[j]=0.0; md2[k]=u;  }
md3=ObratMatr(fx, md2, n, snmnr);
for (k=0; k<n; k++) { s=0.0;
for (j=0; j<n; j++) s=s+md3[k][j]*f[j]; dx[k]=s; }
for (j=0; j<n; j++) { u=md2[j]; if (u) delete[]u; u=md3[j]; if (u) delete []u; }
if (md2) delete[]md2; if (md3) delete[]md3;
return dx; }
double **ObratMatr(double **mdp, double **mdv, int n, char *snmnr)
{ 
	int k=0, j=0, p=0; 
	double perv=0.0, vtor=0.0, *u=NULL, **md1=new double*[n], **md2=new double*[n], e=1e-9;
for (j=0; j<n; j++) { u=new double[n]; if (!u) { cout << snmnr << endl; k=getchar(); exit(1); } md1[j]=u; 
u=new double[n]; if (!u) { cout << snmnr << endl; k=getchar(); exit(1); } md2[j]=u; }
for (k=0; k<n; k++) for (j=0; j<n; j++) { md1[k][j]=mdp[k][j]; md2[k][j]=mdv[k][j]; }
for (k=0; k<n-1; k++) { perv=md1[k][k];
	if (fabs(perv)>e)
    { for (j=0; j<n; j++) { md1[k][j]=md1[k][j]/perv; md2[k][j]=md2[k][j]/perv; } } perv=1.0;
	for (j=k+1; j<n; j++) {
        vtor=md1[j][k];
		if (fabs(vtor)>e)
        { for (p=0; p<n; p++) { md1[j][p]=md1[j][p]*perv/vtor-md1[k][p]; 
		md2[j][p]=md2[j][p]*perv/vtor-md2[k][p]; } } } }
k=n-1; while (k>0) {
    perv=md1[k][k];
	if (fabs(perv)>0.0)
    { j=n-1; while (j>=0) {
		md1[k][j]=md1[k][j]/perv; md2[k][j]=md2[k][j]/perv; j--; }	}
    j=k-1; while (j>=0) {
        vtor=md1[j][k];
		if (fabs(vtor)>0.0)
        { p=n-1; while (p>=0) {
            md1[j][p]=md1[j][p]-md1[k][p]*vtor; 
			md2[j][p]=md2[j][p]-md2[k][p]*vtor; p--; } } j--; } k--; }
for (j=0; j<n; j++) { u=md1[j]; if (u) delete[]u; }
return md2; }
double vychNevyaz(double **ms1, double **ms4, double *tet, double *ssc, int n)
{ 
	int k=0, j=0; double *ts1=NULL, *ts4=NULL, m=0.0, t=0.0;
	for (k=0; k<n; k++) { ts1=ms1[k]; ts4=ms4[k]; m=0.0;
	for (j=0; j<n; j++) m=m+ts1[j]*tet[j]+ts4[j]*pow(tet[j],4.0);
	t=t+(m-ssc[k]); } 
return t; }
double NormaSystFunct(double *tet, double *pdt, int n, double **ms1, double **ms4, double s)
{ 
	double *fst=NULL, *ts1=NULL, *ts4=NULL;
	double m1=0.0, m2=0.0, m=0.0, t1=0.0, m0=0.0; 
	int k=0, j=0;
	for (k=0; k<n; k++) { 
		ts1=ms1[k]; ts4=ms4[k]; m0=0.0;
	for (j=0; j<n; j++) { 
		t1=tet[j]+s*pdt[j];
		m1=ts1[j]*t1;
		m2=ts4[j]*pow(t1,4.0);
		m0=m0+fabs(m1)+fabs(m2); } 
	m=m+fabs(m0); }
return m; }
double *reshnewtrafs(double ***protv, double ***prots, double *hlr, double *hrl, double *Tm, double *Am, int kost, double *ktptk, double *ete, int cem, 
	double *ktpvo, double *tvo, int dmkvoz, double htk, double hvp, double *te, double tora, double *laefm, int koelvyma, double tmax, 
	double tmin, double tolob, char *snm, double qobsh, double *Rm, double hlr0, double hrl0, int vybves, double por, double *stch, 
	double hk, double p1, double p2, double dkosp, double dkoal, int vybmar, int vp, double *ktr, int dmktr, double *dpctp, int vfqi, 
	int kP, int kolPer)
{
	int k=0, j=0, m=200, n=2*kost, w=0, u=4; kolPer=0;
	double **mJ=NULL, *y=NULL, *po=NULL, e=1e-1, hf=1e0, ra=hf, *temk=new double[n], pr=0.0, tk=0.0;
	double *temn=new double[n], *delta=NULL, s=0.0, t=0.0, *temp=new double[n], tmaxi=5e3;
	if ((!temn) || (!temk) || (!temp)) { cout << snm << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) { temk[j]=te[j]; temn[j]=temk[j]; temp[j]=temk[j]; } //for (j=0; j<n; j++) cout << "temk ( " << j << " ) = " << temk[j] << "\t"; cout << endl; //for (j=0; j<kost; j++) cout << "R ( " << j << " ) = " << Rm[j] << "\tT = " << Tm[j] << "\tA = " << Am[j] << "\t"; cout << endl; j=getchar(); //for (w=0; w<u; w++) { cout << "Matritsa " << w << endl; PrintMatr(protv[w], kost); } j=getchar(); //for (w=0; w<u; w++) { cout << "Matritsa " << w << endl; PrintMatr(prots[w], kost); } j=getchar(); 
	k=0; ra=hf;
	while ((ra>e) && (k<m)) {
	for (j=0; j<n; j++) temn[j]=temk[j]; pr=tk;
	mJ=poiskMatrYakobi(temn, n, Rm, Tm, Am, protv, prots, hlr0, hrl0, vybves, por, htk, hvp, kost, ete, ktptk, cem, ktpvo, tvo, snm, 
		stch, dmkvoz, qobsh, hk, p1, p2, dkosp, dkoal, vybmar, vp, ktr, dmktr, dpctp, laefm, vfqi, kP, kolPer); 
	y=opredMatrF(temn, n, Rm, Tm, Am, protv, prots, hlr0, hrl0, vybves, por, htk, hvp, kost, ete, ktptk, cem, ktpvo, tvo, snm, 
		stch, dmkvoz, qobsh, hk, p1, p2, dkosp, dkoal, vybmar, vp, ktr, dmktr, dpctp, laefm, vfqi, kP, kolPer); //for (j=0; j<n; j++) cout << "y ( " << j << " ) = " << y[j] << "\t"; cout << endl; j=getchar();
	j=provDeterNul(mJ, n, snm); if (!j) { cout << "Determinant = " << j << endl; PrintMatr(mJ, n); j=getchar(); k++; continue; } 
	delta=MetodGaussa(mJ, y, n, snm); //delta=reshMetObrMatr(mJ, y, n, snm); delta=reshMetKram(mJ, y, n, snm); //s=0.0; t=0.0; for (j=0; j<n; j++) { s=s+temp[j]; t=t+hf; } s=s/t; if (fabs(s)>tmaxi) { cout << "k = " << k << "\t"; k++; continue; }
	for (j=0; j<n; j++) { temk[j]=temn[j]-delta[j]; temp[j]=temn[j]+delta[j]; } //for (j=0; j<n; j++) cout << "temp ( " << j << " ) = " << temp[j] << "\t"; cout << endl; j=getchar();
	if (y) delete[]y;
	y=opredMatrF(temp, n, Rm, Tm, Am, protv, prots, hlr0, hrl0, vybves, por, htk, hvp, kost, ete, ktptk, cem, ktpvo, tvo, snm, 
		stch, dmkvoz, qobsh, hk, p1, p2, dkosp, dkoal, vybmar, vp, ktr, dmktr, dpctp, laefm, vfqi, kP, kolPer);
	ra=0.0;
	for (j=0; j<n; j++) {
		ra=ra+fabs(y[j]);
		po=mJ[j]; if (po) delete[]po; }
	s=0.0; t=0.0; for (j=0; j<n; j++) { s=s+temp[j]; t=t+hf; } s=s/t; tk=s; //cout << "tsr = " << s << "\tra = " << tk-pr << "\t"; 
	if ((fabs(s)>tmaxi) || (fabs(tk-pr)<e)) break;
	if (mJ) delete[]mJ; if (delta) delete[]delta; if (y) delete[]y;
	k++; if (!(k%10)) kolPer++; } 
	if (fabs(s)<tmaxi) { //for (j=0; j<n; j++) cout << "temn ( " << j << " ) = " << temn[j] << "\ttemk = " << temk[j] << "\t"; cout << endl;
		printf("shagi = %d\tra = %0.4le\t", k, ra); } 
	for (j=0; j<n; j++) { s=temn[j]; if (fabs(s)>tmaxi) s=0.0; temn[j]=s; 
	s=temk[j]; if (fabs(s)>tmaxi) s=0.0; temk[j]=s; }
	if (temn) delete[]temn; if (temp) delete[]temp;
	return temk;
}
double *koefPrib(double *ktp, double *te, int le, char *snm)
{
	int k=0, kem=3, w=0, v=1, u=2, r=0, j=0; 
	double **A=new double*[kem], *AA=NULL, *b=new double[kem], *t=new double[kem];
	double e=1e-15, yx2=0.0, yx=0.0, p=0.0, hf=1e0, *ko=NULL;
	double x4=0.0, x3=x4, x2=x3, x=x2, y=x, de=x, de1=x, de2=x, de3=x;
	if ((!A) || (!b) || (!t)) { cout << snm << endl; k=getchar(); exit(1); }
	for (k = 0; k < kem; k++) { AA = new double[kem]; 
	if (AA) A[k] = AA; else { cout << snm << endl; k=getchar(); exit(1); } } 
	for (k = 0; k < le; k++) {
		yx2 = yx2 + ktp[k] * pow(te[k], 2e0); 
		yx = yx + ktp[k] * te[k]; 
		y = y + ktp[k];
		x4 = x4 + pow(te[k], 4e0); 
		x3 = x3 + pow(te[k], 3e0); 
		x2 = x2 + pow(te[k], 2e0); 
		x = x + te[k]; }
	k=0; b[k] = yx2; k++; b[k] = yx; k++; b[k] = y; p = 0.0; for (k = 0; k < le; k++) p = p + hf;
	A[w][w]=x4; A[w][v]=x3; A[w][u]=x2; A[v][w]=x3; A[v][v]=x2; A[v][u]=x; A[u][w]=x2; A[u][v]=x; A[u][u]=p;
	ko=reshMetKram(A, b, kem, snm); k=0; j=kem-1; p=ko[k]; ko[k]=ko[j]; ko[j]=p; //for (k=0; k<kem; k++) cout << "ko ( " << k << " ) = " << ko[k] << "\t"; cout << endl; //if (ko) delete[]ko; //ko=MetodGaussa(A, b, kem, snm); k=0; j=kem-1; p=ko[k]; ko[k]=ko[j]; ko[j]=p; //for (k=0; k<kem; k++) cout << "ko ( " << k << " ) = " << ko[k] << "\t"; cout << endl; 
	if (b) delete[]b; if (t) delete[]t; 
	for (k=0; k<kem; k++) { AA=A[k]; if (AA) delete[]AA; } if (A) delete[]A;
	return ko;
}
double *usrednen(double *ao, double *usr, int sao, int kost)
{
	double s=0.0, k=0.0, minz=1e-7, hf = 1e0; int j=0;
	for (j=0; j<kost; j++) if (fabs(usr[j])>minz) { s=s+usr[j]; k=k+hf; }
	if (fabs(k)>minz) s=s/k; else s=0.0; ao[sao]=s; //средняя величина
	return ao;
}
double opredKTPTKToch(double *ktptks, double *te, double temp, int ce)
{
	int n=ce, f=1, p=0, k=0, nn=0, nk=n-1; 
	double e=1e-6, ktp=0.0, ko=0.0, x1=0.0, x2=0.0, y1=0.0, y2=0.0, b=0.0, dt=0.0; 
	if ((temp>=te[nn]) && (temp<te[nk-1])) {
	for (k=0; k<(n-1); k++)
		if ((te[k+1]>temp) && (f>0) && (te[k]<=temp)) { p=k+1; f=nn; break; } }
		else if ((temp>=te[nk-1]) && (temp<te[nk])) { p=nk; f=nn; }
		else if (temp<te[nn]) { p=1; f=nn; }
		else if (temp>=te[nk]) { p=nk; f=nn; }
	if ((!f) && (p>nn)) {
		x2=te[p];
		x1=te[p - 1];
		dt = x2 - x1;
		if (fabs(dt) > e) {
			y2=ktptks[p];
			y1=ktptks[p - 1];
			b=y1;
			if ((nk==p) && (temp>te[p])) b=y2;
			else if ((nk==p) && (temp<=te[p])) b=y1;
			else if ((p==1) && (temp<te[nn])) b=y1;
			ko = (y2 - y1) / dt;
			ktp = b + ko*(temp - x1); }
		else ktp=0.0; } else ktp=0.0;
	return ktp;
}
double *GAUSS(double **A, double *B, int chst, double *X)
{ int N=chst, I=0, *L=new int[N], K=0, J=0, LK=0;
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
	  if (L) delete []L; if (S) delete []S; 
	  return X; }
void zapisvfile(double *vyvod, int dlina, char *nazf)
{ 	if (dlina>0) { int k; FILE *fo = fopen(nazf, "a+"); 
	if (!fo) { cout << "File is not open!" << endl; k=getchar(); exit(1); } cout << "dl_zap" << dlina << endl;
	for (k=0; k<dlina; k++) fprintf(fo,"%0.15le\n",vyvod[k]); //fprintf(fo,"%\n\n\n"); 
	fclose(fo); } }
double *koefPribl(double *ktp, double *te, int le, char *snm)
{
	int k=0, kem=3, w=0, v=1, u=2; 
	double **A=new double*[kem], *AA=NULL, *b=new double[kem];
	double *ko=new double[kem], yx2=0.0, yx=0.0, p=0.0, hf=1e0;
	if ((!A) || (!b) || (!ko)) { cout << snm << endl; k=getchar(); exit(1); }
	double x4=0.0, x3=x4, x2=x3, x=x2, y=x, de=x, de1=x, de2=x, de3=x;
	if (A) { for (k=0; k<kem; k++) { AA=new double[kem]; 
	if (AA) A[k]=AA; else { cout << snm << endl; k=getchar(); exit(1); } } }
	else { cout << snm << endl; k=getchar(); exit(1); } //for (k=0; k<le; k++) cout << "y ( " << k << " ) = " << ktp[k] << "\tx = " << te[k] << "\t"; cout << endl; 
	for (k=0; k<le; k++) {
		yx2=yx2+ktp[k]*pow(te[k], 2e0); 
		yx=yx+ktp[k]*te[k]; 
		y=y+ktp[k];
		x4=x4+pow(te[k], 4e0); 
		x3=x3+pow(te[k], 3e0); 
		x2=x2+pow(te[k], 2e0); 
		x=x+te[k]; }
	k=0; b[k]=yx2; k++; b[k]=yx; k++; b[k]=y; p=0.0; for (k=0; k<le; k++) p=p+hf;
	A[w][w]=x4; A[w][v]=x3; A[w][u]=x2; A[v][w]=x3; A[v][v]=x2; A[v][u]=x; A[u][w]=x2; A[u][v]=x; A[u][u]=p;
	//for (k=0; k<kem; k++) cout << "b ( " << k << " ) = " << b[k] << "\t"; cout << endl; 
	de=A[w][w]*(A[u][u]*A[v][v]-A[u][v]*A[v][u]);
	de=de-A[w][v]*(A[v][w]*A[u][u]-A[u][w]*A[v][u]);
	de=de+A[w][u]*(A[v][w]*A[u][v]-A[u][w]*A[v][v]);
	de1=b[w]*(A[v][v]*A[u][u]-A[u][v]*A[v][u]);
	de1=de1-A[w][v]*(b[v]*A[u][u]-b[u]*A[v][u]);
	de1=de1+A[w][u]*(b[v]*A[u][v]-b[u]*A[v][v]);
	de2=A[w][w]*(b[v]*A[u][u]-b[u]*A[v][u]);
	de2=de2-b[w]*(A[v][w]*A[u][u]-A[u][w]*A[v][u]);
	de2=de2+A[w][u]*(A[v][w]*b[u]-A[u][w]*b[v]);
	de3=A[w][w]*(A[v][v]*b[u]-A[u][v]*b[v]);
	de3=de3-A[w][v]*(A[v][w]*b[u]-A[u][w]*b[v]);
	de3=de3+b[w]*(A[v][w]*A[u][v]-A[u][w]*A[v][v]);
	if (fabs(de)>0.0) { k=2; ko[k]=de1/de; k--; ko[k]=de2/de; k--; ko[k]=de3/de; } else for (k=0; k<kem; k++) ko[k]=0.0;
	if (b) delete[]b; for (k = 0; k < kem; k++) { AA = A[k]; if (AA) delete[]AA; } //for (k=0; k<kem; k++) cout << "ko_Prbl ( " << k << " ) = " << ko[k] << "\t"; cout << endl;
	return ko;
}
double **poiskMatrYakobi(double *te, int dlmt, double *Ref, double *Tra, double *Ab, double ***protv, double ***prots, double hlr0, double hrl0, 
	int vybves, double por, double htk, double hvp, int kost, double *ete, double *ktptk, int cem, double *ktpvo, double *vte, char *snm, 
	double *stch, int dmkvoz, double qobsh, double hk, double p1, double p2, double dkosp, double dkoal, int vybmar, int vp, double *ktr, 
	int dmktr, double *dpctp, double *laefm, int vfqi, int kolPok, int kolPer)
{
	int k=0, n=dlmt, j=0, q=0, w=0, x=2, f=0;
	double r=0.0, e=1e-4, tx=r, txs=r, *hlr=NULL, *hrl=NULL, *hlrs=NULL, *dqdt=NULL;
	double *hrls=NULL, **mJ=new double*[n], dpctpe=r, t=r, s=r, **mu=NULL;
	double hf=1e0, latk=r, *dqtdt=NULL, lav=r, lae=r, d=r, y=r;
	double *hes=new double[n], *he=new double[n]; //cout << "n = " << n << "\tks = " << kost << "\tkolPer = " << kolPer << "\t"; //for (k=0; k<n; k++) cout << "tem ( " << k << " ) = " << te[k] << "\t";
	if ((!mJ) || (!hes) || (!he)) { cout << snm << endl; k=getchar(); exit(1); }
	d=0.0; for (k=0; k<kost; k++) d=d+hf;
	for (k=0; k<n; k++) {
		tx=te[k];
		latk=opredKTPTKToch(ktptk, ete, tx, cem);
		lav=opredKTPTKToch(ktpvo, vte, tx, dmkvoz);
		dpctpe=opredKTPTKToch(dpctp, ete, tx, cem);
		dqdt=new double[n];
		if (!dqdt) { cout << snm << endl; k=getchar(); exit(1); }
		for (j=0; j<n; j++) dqdt[j]=0.0;
		mJ[k]=dqdt;
		if (k<(n-1)) { if (!(k%2)) { r=latk/htk; dqdt[k]=r; dqdt[k+1]=-r; }
		else { s=(hf-dpctpe)*lav/hvp; r=dpctpe*latk/htk/2e0; r=r+s; dqdt[k]=r; dqdt[k+1]=-r; } } 
		else { lae=opredKTPTKToch(laefm, ete, tx, cem); j=0; r=lae/(d*htk+(d-hf)*hvp); dqdt[j]=r; dqdt[k]=-r; } 
	hes[k]=0.0; he[k]=0.0; }
	for (k=0; k<(n-1); k++) {
	dqdt=mJ[k]; f=0;
		for (q=0; q<n; q++) {
	txs=te[q]+e;
	dpctpe=opredKTPTKToch(dpctp, ete, txs, cem);
	vp=2; hlrs=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, txs, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, 
		dmkvoz, qobsh, hk, p1, p2, txs, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, e, f, q%2, vfqi, kolPok); //if (!q) { for (j=0; j<kost; j++) cout << "hlrs ( " << j << " ) = " << hlrs[j] << "\t"; cout << endl << endl; }
	vp=3; hrls=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, txs, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, 
		dmkvoz, qobsh, hk, p1, p2, txs, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, e, f, q%2, vfqi, kolPok); //if (!q) { for (j=0; j<kost; j++) cout << "hrls ( " << j << " ) = " << hrls[j] << "\t"; cout << endl << endl; }
	hes=perPlotTeplPot(hlrs, hrls, kost, snm, protv, n, kolPer); //if (!q) { for (j=0; j<n; j++) cout << "hes ( " << j << " ) = " << hes[j] << "\t"; cout << endl << endl; }
	if (hlrs) delete[]hlrs; if (hrls) delete[]hrls;
	r=0.0; s=r; t=r; tx=te[q]; 
	dpctpe=opredKTPTKToch(dpctp, ete, tx, cem);
	vp=2; hlr=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tx, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, 
		dmkvoz, qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, r, f, q%2, vfqi, kolPok); //if (!q) { for (j=0; j<kost; j++) cout << "hlr ( " << j << " ) = " << hlr[j] << "\t"; cout << endl << endl; }
	vp=3; hrl=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tx, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, 
		dmkvoz, qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, r, f, q%2, vfqi, kolPok); //for (j=0; j<kost; j++) cout << "hrl ( " << j << " ) = " << hrl[j] << "\t"; cout << endl << endl; //if (!q) { for (j=0; j<kost; j++) cout << "hlr ( " << j << " ) = " << hlr[j] << "\thrl = " << hrl[j] << "\t"; cout << endl << endl; } j=getchar();
	he=perPlotTeplPot(hlr, hrl, kost, snm, protv, n, kolPer); //if (!q) { for (j=0; j<n; j++) cout << "he ( " << j << " ) = " << he[j] << "\t"; cout << endl << endl; } j=getchar();
	if (hlr) delete[]hlr; if (hrl) delete[]hrl; //0 - hrl, 1 - hlr
	if ((!k) && ((q+x+1)<n)) { r=hes[q+1]-he[q+1]; s=hes[q]-he[q]; t=hes[q+x]-he[q+x]; y=r; r=r-(s+t); } //r=dhlr(k)/dt(+); s=dhrl(k)/dt(-); t=dhrl(k+1)/dt(-)
	else { if ((k==(n-2)) && ((q+1)<n) && ((q-x+1)>0)) { r=hes[q+1]-he[q+1]; s=hes[q-x+1]-he[q-x+1]; t=hes[q]-he[q]; y=r; r=r+s-t; } //r=dhlr(k)/dt(+); s=hrl(k-1)/dt(+); t=hrl(k)/dt(-);
	else { if ((q%2) && (k>0) && (k<(n-2)) && ((q+x+1)<n)) { r=hes[q+1]-he[q+1]; s=hes[q+x]-he[q+x]; y=r; r=r-s; } //r=dhlr(k)/dt(+); s=hrl(k+1)/dt(-); 
	else if ((!(q%2)) && (k>0) && (k<(n-2)) && ((q+x+1)<n) && ((q-x+1)>0)) { r=hes[q-x+1]-he[q-x+1]; s=hes[q+x]-he[q+x]; y=r; r=r-s; } } } //r=dhlr(k-1)/dt(+); s=hrl(k+1)/dt(-); 
	r=r/(txs-tx); //if (!q) cout << "r = " << r << "\tq = " << q << "\tk = " << k << "\ts = " << s << "\tt = " << t << "\ty = " << y << "\ttxs = " << txs << "\ttx = " << tx << endl;
	dqdt[q]=dqdt[q]+r; if (q%2) f++; } 
		mJ[k]=dqdt; }
	if (hes) delete[]hes; if (he) delete[]he; //PrintMatr(mJ, n); k=getchar();
	return mJ;
}
double *opredMatrF(double *te, int dlmt, double *Ref, double *Tra, double *Ab, double ***protv, double ***prots, double hlr0, double hrl0, int vybves, 
	double por, double htk, double hvp, int kost, double *ete, double *ktptk, int cem, double *ktpvo, double *vte, char *snm, double *stch, 
	int dmkvoz, double qobsh, double hk, double p1, double p2, double dkosp, double dkoal, int vybmar, int vp, double *ktr, int dmktr, 
	double *dpctp, double *laefm, int vfqi, int kolPok, int kolPer)
{
	int k=0, n=dlmt, j=0, q=0;
	double r=0.0, t=r, s=r, hf=1e0, latk=r, lav=r, lae=r, d=r, dpctpe=r, mo=-hf, *he=NULL;
	double tx=r, *hlr=NULL, *hrl=NULL, *vf=new double[n], dd=r, ts=r, *Tssk=new double[kost];
	double *Ts=new double[kost], *Tss=new double[kost], Tn=te[k], *Tsk=new double[kost];
	if ((!vf) || (!Ts) || (!Tss) || (!Tssk) || (!Tsk)) { cout << snm << endl; k=getchar(); exit(1); }
	d=0.0; for (k=0; k<kost; k++) { d=d+hf; Ts[k]=0.0; Tss[k]=0.0; Tsk[k]=0.0; Tssk[k]=0.0; }
	q=0; dd=0.0; s=dd; 
	for (k=0; k<n; k++) { if (!(k%2)) Ts[q]=te[k]; else { Tss[q]=te[k]; q++; } 
	dd=dd+hf; s=s+te[k]; }
	s=s/dd; ts=s;
	k=0; Tsk[k]=Ts[k]; q=0;
	for (k=1; k<(2*kost); k++) {
	t=opredTempStenFragm(Tn, k, ktpvo, vte, ete, ktptk, cem, dmkvoz, htk, hvp, qobsh, mo); 
	t=Tn-t;
	if (k%2) { Tssk[q]=t; q++; } else Tsk[q]=t; }
	for (k=0; k<n; k++) vf[k]=-qobsh;
	q=0;
	for (k=0; k<n; k++) {
		tx=te[k];
		latk=opredKTPTKToch(ktptk, ete, tx, cem);
		lav=opredKTPTKToch(ktpvo, vte, tx, dmkvoz);
		dpctpe=opredKTPTKToch(dpctp, ete, tx, cem);
		if (k<(n-1)) { if (!(k%2)) { r=latk/htk; t=Ts[q]-Tss[q]; vf[k]=vf[k]+r*t; }
		else { s=(hf-dpctpe)*lav/hvp; t=Tss[q]-Ts[q+1]; vf[k]=vf[k]+s*t; r=dpctpe*latk/htk; t=Tssk[q]-Tsk[q+1]; vf[k]=vf[k]+r*t; } } 
		else { lae=opredKTPTKToch(laefm, ete, tx, cem); j=0; r=lae/(d*htk+(d-hf)*hvp); t=Ts[j]-Tss[k]; vf[k]=vf[k]+r*t; } 
	if (k%2) q++; }
	dpctpe=opredKTPTKToch(dpctp, ete, ts, cem); r=0.0;
	vp=2; hlr=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tx, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, 
		dmkvoz, qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, r, q, -1, vfqi, kolPok); //for (j=0; j<kost; j++) cout << "hlr ( " << j << " ) = " << hlr[j] << "\t"; cout << endl << endl;
	vp=3; hrl=RasIzlSerStenNacSok(Ref, Tra, Ab, protv, prots, tx, hlr0, hrl0, por, dpctpe, htk, hvp, kost, ete, ktptk, cem, ktpvo, vte, snm, stch, 
		dmkvoz, qobsh, hk, p1, p2, tx, dkosp, dkoal, vybves, vybmar, vp, ktr, dmktr, r, q, -1, vfqi, kolPok); //for (j=0; j<kost; j++) cout << "hrl ( " << j << " ) = " << hrl[j] << "\t"; cout << endl << endl;
	he=perPlotTeplPot(hlr, hrl, kost, snm, protv, n, kolPer);
	for (j=0; j<n; j++) if (j<kost) hrl[j]=he[j]; else hlr[j-kost]=he[j];
	q=0;
	for (k=0; k<(n-1); k++) {
		tx=te[k];
	if (!(k%2)) r=hlr[q]-hrl[q]; else r=hlr[q]-hrl[q+1];
		vf[k]=vf[k]+r; //cout << "r = " << r << "\t";
	if (k%2) q++;
	}
	for (k=0; k<n; k++) vf[k]=-vf[k];
	if (hlr) delete[]hlr; if (hrl) delete[]hrl; 
	if (Ts) delete[]Ts; if (Tss) delete[]Tss; 
	if (Tsk) delete[]Tsk; if (Tssk) delete[]Tssk;
	return vf;
}
int provDeterNul(double **A, int n, char *snm)
{ 
	int k=0, i=0, j=0, w=0, f=0;
	double **kA=new double*[n], *p=NULL, d=0.0, c=0.0, hf=1e0, r=0.0, e=1e-10;
	if (!kA) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) { p=new double[n]; if (!p) { cout << snm << endl; k=getchar(); exit(1); } 
	for (j=0; j<n; j++) p[j]=A[k][j]; kA[k]=p; }
	for (k=0; k<(n-1); k++) {
		d=kA[k][k];
		//-----
		if (fabs(d)<e) { f=0;
			for (i=k+1; i<n; i++) {
				if (fabs(kA[i][k])>e) { f=1; w=i; break; } } 
				if (f) { 
		for (i=0; i<n; i++) { r=kA[k][i]; kA[k][i]=kA[w][i]; kA[w][i]=r; } 
		d=kA[k][k]; } else continue; }
		//-----
		for (j=k+1; j<n; j++) { c=kA[j][k]; 
		//-----
		if (fabs(c)>e) { c=d/c;
		for (i=k; i<n; i++) kA[j][i]=kA[j][i]*c-kA[k][i]; } else continue; } }
	f=1; for (k=0; k<n; k++) { r=0.0; 
	for (j=0; j<n; j++) r=r+fabs(kA[k][j]); 
	if (fabs(r)<e) { f=0; break; } }
	for (k=0; k<n; k++) { p=kA[k]; if (p) delete[]p; } if (kA) delete[]kA;
	return f;
}
double *perPlotTeplPot(double *hlr, double *hrl, int kost, char *snm, double ***mpov, int n, int kolPer)
{
	int i=0, w=0, x=w+1, y=x+1, z=y+1, j=2, v=0;
	double **mu0=mpov[w], **mu1=mpov[x], **mu2=mpov[y], **mu3=mpov[z], *he=new double[n];
	double w0=0.0, w1=0.0, w2=0.0, w3=0.0, wo=0.0, e=1e-9, srl=0.0, slr=0.0, l=0.0, r=0.0;
	if (!he) { cout << snm << endl; j=getchar(); exit(1); }
	if (kolPer>0) {
	double *hlrn=new double[kost], *hrln=new double[kost];
	double *hlrs=new double[kost], *hrls=new double[kost];
	if ((!hlrn) || (!hrln) || (!hrls) || (!hlrs)) { cout << snm << endl; j=getchar(); exit(1); }
	for (i=0; i<kost; i++) { hlrn[i]=hlr[i]; hrln[i]=hrl[i]; hlrs[i]=hlr[i]; hrls[i]=hrl[i]; }
	for (v=0; v<kolPer; v++) {
	for (i=0; i<kost; i++) {
	slr=0.0; srl=0.0; w=0; x=w+1; y=x+1; z=y+1;
	for (j=0; j<kost; j++) { 
		l=hrl[j]; //то, что идет влево
		r=hlr[j]; //то, что идет вправо
		w0=mu0[i][j];
		w2=mu2[i][j];
		w3=mu3[i][j];
		wo=w2+w3;
		if (fabs(wo)>e) slr=slr+w0*(w2*l+r*w3)/wo;
		w1=mu1[i][j];
		w2=mu2[i][j]; //влево
		w3=mu3[i][j]; //вправо
		wo=w2+w3;
		if (fabs(wo)>e) srl=srl+w1*(w2*l+r*w3)/wo; } 
	hlrn[i]=slr; hrln[i]=srl; } 
	for (i=0; i<kost; i++) { hlr[i]=hlrn[i]; hrl[i]=hrln[i]; } }
	for (i=0; i<kost; i++) { hlr[i]=hlrs[i]; hrl[i]=hrls[i]; }
	if (hlrn) delete[]hlrn; if (hrln) delete[]hrln;	if (hlrs) delete[]hlrs; if (hrls) delete[]hrls; }
	w=0; for (j=0; j<n; j++) { if (j%2) { he[j]=hlr[w]; w++; } else he[j]=hrl[w]; }
	return he;
}