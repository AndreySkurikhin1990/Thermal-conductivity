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
const int dmkoosck = 14, dkosckl = 6, vybkvi=4;
const double tek0 = 273.15, tnosck = 3e2, dtosck = 1e2, tnack = 2e2, detek = 1e2;
//-------------
double opredKTPTKTochSha(double *, double *, double, int);
void oprsodoxkvi();
double *koefoslab(double, double, double, double *, int);
double *oslaintestepcher(int, int, double *);
double *RasKorOV(double, double, double, double, double, double);
double UsredMasOV(double **, int, int);
double vychopred(double  **, int);
double **osvpam(double **, int);
double Determinant(double **, int);
double vychopred(double  **, int);
double *reshMetKram(double **, double *, int);
double **polmat(double **, double *, int, int);
double **GetMatr(double **, double **, int, int, int);
//-------------
double *dkosckm = NULL, *dkosckt = NULL, *tkusck = NULL;
double *kusck = NULL, T, *sodoxkvi, *etek;
int cemk=11;
//-------------
void main()
{	double wmg, wsi, wal, dkusct, dkosce, ume; int k; 
	dkosckm = new double[dkosckl]; dkosckt = new double[dkosckl]; etek=new double[cemk]; sodoxkvi=new double[6];
	if ((!dkosckm) || (!dkosckt) || (!sodoxkvi)) { cout << "No memory" << endl; k = getchar(); exit(1); }
	oprsodoxkvi(); wal=sodoxkvi[0]; wsi=sodoxkvi[1]; wmg=sodoxkvi[2]; delete []sodoxkvi; 
	tkusck = new double[dmkoosck]; kusck = new double[dmkoosck];
	if ((!tkusck) || (!kusck)) { cout << "No memory" << endl; k = getchar(); exit(1); }
	tkusck[0] = tnosck; for (k = 1; k < dmkoosck; k++) tkusck[k] = tkusck[k - 1] + dtosck; //cout << "wsi = " << wsi << "\twal = " << wal << "\twmg = " << wmg << endl;
	kusck = koefoslab(wmg, wsi, wal, tkusck, dmkoosck); 
	etek[0]=tnack+tek0; for (k=1; k<cemk; k++) etek[k] = etek[k-1]+detek; //for (k=0; k<dmkoosck; k++) cout << "dkusct = " << kusck[k] << "\ttkusck = " << tkusck[k] << "\t"; cout << endl; //for (k=0; k<dkosckl; k++) cout << "dkoscem = " << dkosckm[k] << "\tdkoscet = " << dkosckt[k] << "\t"; cout << endl;
for (k=0; k<cemk; k++) { T=etek[k]; 
dkusct=opredKTPTKTochSha(kusck, tkusck, T, dmkoosck); 
if (dkusct>1e0) dkusct=1e0; if (dkusct<0.0) dkusct=0.0; cout << "dkusct = " << dkusct << endl;
dkosce=opredKTPTKTochSha(dkosckm, dkosckt, T, dkosckl); 
if (dkosce>1e0) dkosce=1e0; if (dkosce<0.0) dkosce=0.0;  cout << "dkosce = " << dkosce << endl;
ume=dkosce*dkusct; cout << "ume = " << ume << endl; 
}
k=getchar();
}
double opredKTPTKTochSha(double *ktptks, double *te, double temp, int ce)
{
	int n = ce, f = 1, p = 0, k; 
	double e=1e-4, ktp = 0.0, ko = 0.0, x1=0.0, x2=0.0, y1=0.0, y2=0.0, b=0.0, dt=0.0; 
	if ((temp>=te[0]) && (temp<=te[n-1])) {
	for (k = 0; k<n; k++)
		if ((te[k] >= temp) && (f>0)) { p = k; f = 0; } 
	if (!p) {
		if (!f) p = 1; else {
			p = n - 1; f = 0;
		}
	} }
	else if (temp<te[0]) { p=1; f=0; }
	else if (temp>te[n-1]) { p=n-1; f=0; } //cout << "p = " << p << "\tf = " << f << endl;
	if ((!f) && (p>0)) {
		x2=te[p];
		x1=te[p - 1];
		dt = x2 - x1;
		if (fabs(dt) > e) {
			y2=ktptks[p];
			y1=ktptks[p - 1];
			b=y1;
			if ((n-1)==p) b=y2;
			ko = (y2 - y1) / dt;
			ktp = b + ko*(temp - x1);
		}
	}
	return ktp;
}
void oprsodoxkvi()
{ int k=0; double wal, wsi, wmg, salok, smgok, ssiok, ko=1e-2;
double tnd = 6e2, dtd = 2e2, tm; dkosckt[0] = tnd; 
	for (k = 1; k < dkosckl; k++) dkosckt[k] = dkosckt[k - 1] + dtd;
	if (vybkvi == 4) {
		salok = 33e0; smgok = 15e0; ssiok = 52e0; 
		wal = 25e0; wsi = 11e0; wmg = 4e1;
	} //КВИ-400
	else if (vybkvi == 5) {
		salok = 34e0; smgok = 11e0; ssiok = 54e0;
		wal = 28e0; wmg = 8e0; wsi = 44e0; 
	} //КВИ-500
	else if (vybkvi == 6) {
		salok = 36e0; smgok = 9e0; ssiok = 55e0;
		wal = 3e1; wsi = 7e0; wmg = 45e0; 
	} //КВИ-600
	else if (vybkvi == 7) {
		salok = 37e0; smgok = 8e0; ssiok = 55e0; 
		wal = 31e0; wmg = 6e0; wsi = 45e0; 
	} //КВИ-700
	else if (vybkvi == 8) {
		salok = 38e0; smgok = 7e0; ssiok = 55e0; 
		wal = 3e1; wmg = 5e0; wsi = 45e0; 
	} //КВИ-800
	else if (vybkvi == 9) {
		salok = 39e0; smgok = 6e0; ssiok = 55e0; 
		wal = 32e0; wmg = 5e0; wsi = 45e0;
	} //КВИ-900
	else if (vybkvi == 10) {
		salok = 39e0; smgok = 6e0; ssiok = 55e0; 
		wal = 32e0; wmg = 4e0; wsi = 45e0;
	} //КВИ-1000
	else { cout << "Net takoy marki KVI!" << endl; k = getchar(); exit(1); }
	salok = salok*ko; smgok = smgok*ko; ssiok = ssiok*ko;
	wal = wal*ko; wsi = wsi*ko; wmg = wmg*ko;
		k = 0; dkosckm[k] = 3e0; k++; dkosckm[k] = 6e0; k++; dkosckm[k] = 7e0; k++; 
		dkosckm[k] = 15e0; k++; dkosckm[k] = 2e1; k++; dkosckm[k] = 26e0;
		if (vybkvi>6) { k=1; dkosckm[k]=7e0; } if (vybkvi>9) { k=2; dkosckm[k]=8e0; }
	for (k = 0; k < dkosckl; k++) {
		tm = dkosckm[k]*ko; dkosckm[k] = 1e0 - tm; //cout << "dko = " <<dkosckm[k] << "\t";
	}
	sodoxkvi[0]=wal; sodoxkvi[1]=wsi; sodoxkvi[2]=wmg; 
	sodoxkvi[3]=salok; sodoxkvi[4]=ssiok; sodoxkvi[5]=smgok; 
}
double *koefoslab(double wmg, double wsi, double wal, double *tere, int n)
{ int lt=n, k; 
double wo=wmg+wsi+wal, kmg=0.0, kal=0.0, ksi=0.0, *kuo=new double[lt], ht=1e0, eps=1e-6;
if (!kuo) { cout << "No memory" << endl; k=getchar(); exit(1); }
double *mgo=new double[lt], *alo=new double[lt], *sio=new double[lt];
if ((!mgo) || (!alo) || (!sio)) { cout << "No memory" << endl; k=getchar(); exit(1); }
sio=oslaintestepcher(lt,0,sio);
alo=oslaintestepcher(lt,1,alo);
mgo=oslaintestepcher(lt,2,mgo); //for (k=0; k<n; k++) cout << "sio = " << sio[k] << "\talo = " << alo[k] << "\tmgo = " << mgo[k] << endl; //cout << "wmg = " << wmg << "\twsio = " << wsi << "\twal = " << wal << "\two = " << wo << endl;
for (k=0; k<n; k++) {
kmg=mgo[k]; kal=alo[k]; ksi=sio[k];
if (kmg<0.0) kmg=0.0; if (kmg>ht) kmg=ht; 
if (kal<0.0) kal=0.0; if (kal>ht) kal=ht; 
if (ksi<0.0) ksi=0.0; if (ksi>ht) ksi=ht;
if (fabs(wo)>eps) 
kuo[k]=(kmg*wmg+kal*wal+ksi*wsi)/wo; 
else kuo[k]=0.0; } //for (k=0; k<n; k++) cout << "kuo = " << kuo[k] << "\ttere = " << tere[k] << endl;
delete []mgo; delete []alo; delete []sio; return kuo; }
double *oslaintestepcher(int dm, int vy, double *sc) //рассчитывает ослабление интегральной степени черноты для SiO2, Al2O3, MgO
{	int n=2, k, j, w=0, p=6, q=10, m, l; double **scv=new double*[p], **hsv=new double*[p];
	double *hs=new double[n], **kor=new double*[q], t0, t1, *vscs=new double[dm], *vscm=new double[dm], *vsca=new double[dm], *po=NULL;
if ((!scv) || (!sc) || (!hsv) || (!hs) || (!kor)) { cout << "No memory" << endl; k=getchar(); exit(1);}
	k=0;
	sc[k]=1e0;   k++; sc[k]=0.94;  k++; sc[k]=0.87;  k++; sc[k]=0.801; k++; sc[k]=0.736; k++; 
	sc[k]=0.676; k++; sc[k]=0.635; k++; sc[k]=0.590; k++; sc[k]=0.567; k++; sc[k]=0.543; k++; 
	sc[k]=0.53;  k++; sc[k]=0.525; k++; sc[k]=0.515; k++; sc[k]=0.507; //степень черноты магнезита
	w=0; scv[w]=sc; sc=new double[dm];
	k=0;
	hs[k]=98e-2; k++; hs[k]=2e-2;
	k=0;
	t0=hs[k]; k++; t1=hs[k];
	k=0;
	hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1);
	hsv[w]=hs; w++; hs=new double[n]; //Магнезит: 98 % - MgO, 2 % - SiO2
	if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1); }
k=0;
sc[k]=1e0;   k++; sc[k]=0.976; k++; sc[k]=0.949; k++; sc[k]=0.905; k++; sc[k]=0.859; k++; 
sc[k]=0.812; k++; sc[k]=0.774; k++; sc[k]=0.737; k++; sc[k]=0.709; k++; sc[k]=0.681; k++; 
sc[k]=0.661; k++; sc[k]=0.639; k++; sc[k]=0.626; k++; sc[k]=0.620; //степень черноты шамота
scv[w]=sc; sc=new double[dm]; 
k=0; 
hs[k]=56e-2; k++; hs[k]=396e-3; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Шамот: 56 % - SiO2, 39,6 % - Al2O3 - шамот, считаем, что магния нет в шамоте
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1e0;    k++; sc[k]=98e-2;  k++; sc[k]=951e-3; k++; sc[k]=92e-2;  k++; sc[k]=883e-3; k++; 
sc[k]=853e-3; k++; sc[k]=821e-3; k++; sc[k]=79e-2;  k++; sc[k]=767e-3; k++; sc[k]=746e-3; k++;
sc[k]=73e-2;  k++; sc[k]=715e-3; k++; sc[k]=705e-3; k++; sc[k]=692e-3; //степень черноты корундошамота
scv[w]=sc; sc=new double[dm]; 
k=0; 
hs[k]=28e-2; k++; hs[k]=7e-1; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Корундошамот: 28 % - SiO2, 70 % - Al2O3 - корундошамот
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1.0000; k++; sc[k]=983e-3; k++; sc[k]=936e-3; k++; sc[k]=867e-3; k++; sc[k]=819e-3; k++; 
sc[k]=721e-3; k++; sc[k]=659e-3; k++; sc[k]=593e-3; k++; sc[k]=541e-3; k++; sc[k]=49e-2;  k++; 
sc[k]=453e-3; k++; sc[k]=429e-3; k++; sc[k]=403e-3; k++; sc[k]=384e-3; //степень черноты каолинового теплоизоляционного кирпича (КТК)
scv[w]=sc; sc=new double[dm];
k=0; 
hs[k]=57e-2; k++; hs[k]=4e-1; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //КТК: 57 % - SiO2, 40 % - Al2O3
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1); }
k=0;
sc[k]=1.000;  k++; sc[k]=984e-3; k++; sc[k]=941e-3; k++; sc[k]=882e-3; k++; sc[k]=813e-3; k++; 
sc[k]=751e-3; k++; sc[k]=695e-3; k++; sc[k]=641e-3; k++; sc[k]=594e-3; k++; sc[k]=558e-3; k++; 
sc[k]=53e-2;  k++; sc[k]=499e-3; k++; sc[k]=479e-3; k++; sc[k]=462e-3; //степень черноты муллита
scv[w]=sc; sc=new double[dm];
k=0; 
hs[k]=28e-2; k++; hs[k]=72e-2; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Муллит: 28 % - SiO2, 72 % - Al2O3
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0;
sc[k]=1e0;    k++; sc[k]=984e-3; k++; sc[k]=953e-3; k++; sc[k]=917e-3; k++; sc[k]=854e-3; k++; 
sc[k]=808e-3; k++; sc[k]=756e-3; k++; sc[k]=711e-3; k++; sc[k]=578e-3; k++; sc[k]=523e-3; k++; 
sc[k]=495e-3; k++; sc[k]=468e-3; k++; sc[k]=448e-3; k++; sc[k]=429e-3; //степень черноты кремнезема
scv[w]=sc; 
k=0; 
hs[k]=985e-3; k++; hs[k]=1e-2; 
k=0; 
t0=hs[k]; k++; t1=hs[k]; 
k=0; 
hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Кремнезем: 98,5 % - SiO2, 1 % - Al2O3
hsv[w]=hs; w++;
for (j=0; j<dm; j++) { k=0; 
for (m=1; m<p; m++)
	for (l=m+1; l<p; l++) {
		kor[k]=RasKorOV(hsv[m][0],hsv[m][1],hsv[l][0],hsv[l][1],scv[m][j],scv[l][j]); k++; }
		vscs[j]=UsredMasOV(kor,q,0); vsca[j]=UsredMasOV(kor,q,1); } 
if (vy==2) { hs=hsv[0]; sc=scv[0]; 
for (k=0; k<dm; k++) vscm[k]=(sc[k]-vscs[k]*hs[1])/hs[0]; } 
for (k=0; k<q; k++) { sc=kor[k]; delete []sc; }
for (k=0; k<p; k++) { hs=hsv[k]; delete []hs; } delete []hsv; //0 - SiO2, 1 - Al2O3, 2 - MgO
for (k=0; k<p; k++) { po=scv[k]; delete []po; } delete []scv;
if (!vy) { delete []vsca; delete []vscm; po=vscs; }
if (vy==1) { delete []vscm; delete []vscs; po=vsca; }
if (vy==2) { delete []vsca; delete []vscs; po=vscm; } 
return po;
}
double *RasKorOV(double a11, double a12, double a21, double a22, double b1, double b2)
{ int l=2, k; double **mas=new double*[l], *st=new double[l], *bs=new double[l], *x, *kor=new double[l];
if ((!mas) || (!st) || (!bs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
st[0]=a11; st[1]=a12; mas[0]=st; st=new double[l];
if (st) { st[0]=a21; st[1]=a22; mas[1]=st; } else { cout << "No memory" << endl; k=getchar(); exit(1);}
bs[0]=b1; bs[1]=b2;
x=reshMetKram(mas,bs,l);
if ((x[0]>=0.0) && (x[1]>=0.0) && (x[0]<=1.0) && (x[1]<=1.0)) {
	kor[0]=x[0]; kor[1]=x[1]; }
else { kor[0]=0.0; kor[1]=0.0; }
for (k=0; k<l; k++) { st=mas[k]; delete []st; } delete []bs;
return kor; }
double UsredMasOV(double **kor, int q, int vy)
{ double s1=0.0, s2=0.0, p=0.0, *s; int j;
for (j=0; j<q; j++) { s=kor[j];
    if ((s[0]>0.0) && (s[1]>0.0)) {
        s1=s1+s[0]; s2=s2+s[1]; p=p+1e0; } }
if (!vy) return s1/p; else return s2/p; }
double *reshMetKram(double **mat, double *ssc, int rm)
{ int n=rm, j;
double d, *m=new double[n], mm, **mmat; if (!m) { cout << "No memory" << endl; j=getchar(); exit(1); } //PrintMatr(mat, rm);
d=vychopred(mat,n); for (j=0; j<n; j++) { mmat=polmat(mat,ssc,j,n); mm=vychopred(mmat,n); 
m[j]=mm/d; /*cout << "x ( " << j << " ) = " << m[j] << endl;*/
mmat=osvpam(mmat,n); } 
return m; }
double **osvpam(double **a, int n)
{ int j; double *b; for (j=0; j<n; j++) { b=a[j]; delete []b; } return NULL; }
double **polmat(double **mat, double *b, int no, int n)
{ int k, j; double **ma=new double*[n], *m; 
if (!ma) { cout << "No memory" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	m=new double[n]; if (!m) { cout << "No memory" << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) {
		if (j==no) m[j]=b[k]; else m[j]=mat[k][j];}
ma[k]=m; } 
return (ma); }
double vychopred(double  **mas, int m) {
  double d=Determinant(mas, m); //Вычисление определителя
  return d; }
double Determinant(double **mas, int m) { // Рекурсивное вычисление определителя
  int i;
  double **p, k=1., d=0.0, *po;
  p=new double*[m]; if (!p) { cout << "No memory" << endl; i=getchar(); exit(1); }
  for (i=0; i<m; i++) { po=new double[m]; if (!po) { cout << "No memory" << endl; i=getchar(); exit(1); } p[i]=po; }
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