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
#define vmivmf 1 //выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный
#define vyfv 1 //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм
#define vyuv 1 //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
#define vysv 0 //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град С
using namespace std;
const double pi = acos(-1e0), enu=1e-3, tocrasov=1e-4, tnoscv = 3e2, tnacv = 2e2, dtoscv = 1e2;
const double tev0=273.15, tna=2e2+tev0, dtv=1e2, pksvv = 0.0;
const double ssiv84 = 467.0*1e-3, salv84 = 129.0*1e-3, smgv84 = 282.0*1e-3;
const double ssiv207 = 458.0*1e-3, salv207 = 14.0*1e-2, smgv207 = 29.0*1e-2;
const int dsov=60, dmkooscv=14, cemv=12, nxtv=30, nnxtv=0;
const int vtvn = 0, isrp=0, vpkf=0, vpmf=0; //0 - старые, 1 - новые
const double por207 = 66.35*1e-2, poris84 = 55.75*1e-2, porin84=81.53*1e-2;
const double poro84=86.61*1e-2, poro16035=84.36*1e-2, por16035=83.97*1e-2;
const char snmv[]="No memory";
//-----
int dlar=1;
char *sndov=NULL, *ssdov=NULL, *szfovk=NULL, *szfov=NULL, *snmov=NULL, *sfnoov=NULL;
char *sppovv=NULL, *sppovvk=NULL, *sdvovv=NULL, *sdvovvk=NULL, *skpovv=NULL, *skpovvk=NULL;
double porver = 0.0;
//-----
double opredKTPTKTochSha(double *, double *, double, int);
void osvpamov(int);
double *Kramers_Kronig_ver(double *, int, double *);
double *epsilnu(double *, int, double *);
double *Koef_Pogl_ver(double *);
double *dliny_voln_ver(double *, int);
double epsisredver(double, double *, double *, int, double *, double *, int);
double usredVelichPlank(double *, double *, double *, double, int, double);
double *izmMasChast(double *, double, int, double *);
double *PokazPrelomAl(double *, double *, double *, double, int, double *);
double *PokazPrelomPribl(double *, double *, double *, int, double *);
double trapz(double *, double *, int);
double *oslaintestepcher(int, int, double *);
double *RasKorOV(double, double, double, double, double, double);
double UsredMasOV(double **, int, int);
double vychopred(double  **, int);
double Determinant(double **, int);
double *reshMetKram(double **, double *, int);
double **polmat(double **, double *, int, int);
double **osvpam(double **, int);
double *koefoslab(double, double, double, double *, int);
double **GetMatr(double **, double **, int, int, int);
void initarrver(int, double, double, double);
double **rasPorpoRazVer(double, int, int, int, int, int);
double **napMasEKTPVerNac();
int napMasEKTPitomNac();
int napMasEKTPkviNac();
int napMasEKTPshaNac();
//-----
int main()
{ int k=0; //int j, jk = nxtv; double smgv=smgv84, ssiv=ssiv84, salv=salv84; initarrver(jk, smgv, ssiv - pksvv, salv + pksvv); 
double **mu=napMasEKTPVerNac(); //k=napMasEKTPshaNac(); //k=napMasEKTPkviNac(); //k=napMasEKTPitomNac();
k=getchar();//int dmi=0; double **mu=rasPorpoRazVer(porver, vyfv, k, vysv, isrp, vpkf), *po; //0 - старые, 1 - новые значения //for (k=0; k<6; k++) { po=mu[k]; if (po) delete[]po; } if (mu) delete[]mu; 
return 0;
}
void novNapMas()
{
	int j, k;
	if ((!vysv) || (vysv==1)) {  //исходный или после повторных измерений
	if (!vyfv) { if (!vpmf) porver = por207; else if (vpmf==1) porver=por16035; }
	else if (vyfv==1) { if (!vpkf) porver = poris84; else if (vpkf==1) porver=porin84; }
	else if (vyfv==2) porver = poro16035; }
	if (vysv==2) { //после обжига
	if ((!vyfv) || (vyfv==2)) porver = poro16035;
	else if (vyfv==1) porver = poro84; }
}
void napstrdir(int vyve, int vm)
{	int k; 
sndov=new char[dsov]; ssdov=new char[dsov]; szfovk=new char[dsov]; 
snmov=new char[dsov]; sfnoov=new char[dsov]; szfov=new char[2*dsov];
for (k=0; k<(2*dsov); k++) szfov[k]='\0';
for (k=0; k<dsov; k++) { sndov[k]='\0'; ssdov[k]='\0'; szfovk[k]='\0'; snmov[k]='\0'; sfnoov[k]='\0'; }
k=0; 
sndov[k]='D'; k++; sndov[k]=':'; k++; sndov[k]='\\';k++; sndov[k]='\\'; k++; sndov[k]='_';  k++; sndov[k]='А';  k++;
sndov[k]='с'; k++; sndov[k]='п'; k++; sndov[k]='и'; k++; sndov[k]='р';  k++; sndov[k]='а';  k++; sndov[k]='н';  k++;
sndov[k]='т'; k++; sndov[k]='у'; k++; sndov[k]='р'; k++; sndov[k]='а';  k++; sndov[k]='\\'; k++; sndov[k]='\\'; k++;
sndov[k]='t'; k++; sndov[k]='m'; k++; sndov[k]='p'; k++; sndov[k]='\\'; k++; sndov[k]='\\'; k++; sndov[k]='\0';
k=0;
ssdov[k]='C';  k++; ssdov[k]=':'; k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++; ssdov[k]='U';  k++; ssdov[k]='s';  k++;
ssdov[k]='e';  k++; ssdov[k]='r'; k++; ssdov[k]='s';  k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++; ssdov[k]='А';  k++;
ssdov[k]='н';  k++; ssdov[k]='д'; k++; ssdov[k]='р';  k++; ssdov[k]='е';  k++; ssdov[k]='й';  k++; ssdov[k]='\\'; k++;
ssdov[k]='\\'; k++; ssdov[k]='D'; k++; ssdov[k]='o';  k++; ssdov[k]='c';  k++; ssdov[k]='u';  k++; ssdov[k]='m';  k++;
ssdov[k]='e';  k++; ssdov[k]='n'; k++; ssdov[k]='t';  k++; ssdov[k]='s';  k++; ssdov[k]='\\'; k++; ssdov[k]='\\'; k++;
ssdov[k]='_';  k++; ssdov[k]='А'; k++; ssdov[k]='с';  k++; ssdov[k]='п';  k++; ssdov[k]='и';  k++; ssdov[k]='р';  k++;
ssdov[k]='а';  k++; ssdov[k]='н'; k++; ssdov[k]='т';  k++; ssdov[k]='у';  k++; ssdov[k]='р';  k++; ssdov[k]='а';  k++;
ssdov[k]='\\'; k++; ssdov[k]='\\';k++; ssdov[k]='t';  k++; ssdov[k]='m';  k++; ssdov[k]='p';  k++; ssdov[k]='\\'; k++;
ssdov[k]='\\'; k++; ssdov[k]='\0';
k=0;
szfovk[k]='V'; k++; szfovk[k]='y'; k++; szfovk[k]='v'; k++; szfovk[k]='o'; k++; szfovk[k]='d'; k++; szfovk[k]='v'; k++; 
szfovk[k]='F'; k++; szfovk[k]='i'; k++; szfovk[k]='l'; k++; szfovk[k]='e'; k++; szfovk[k]='.'; k++; szfovk[k]='t'; k++;
szfovk[k]='x'; k++; szfovk[k]='t'; k++; szfovk[k]='\0';
k=0;
snmov[k]='N'; k++; snmov[k]='o'; k++; snmov[k]='_'; k++; snmov[k]='m'; k++; snmov[k]='e'; k++; snmov[k]='m'; k++; 
snmov[k]='o'; k++; snmov[k]='r'; k++; snmov[k]='y'; k++; snmov[k]='!'; k++; snmov[k]='\0';
k=0;
sfnoov[k]='F';  k++; sfnoov[k]='i';  k++; sfnoov[k]='l';  k++; sfnoov[k]='e';  k++; sfnoov[k]='_';  k++; sfnoov[k]='i'; k++;
sfnoov[k]='s';  k++; sfnoov[k]='_';  k++; sfnoov[k]='n';  k++; sfnoov[k]='o';  k++; sfnoov[k]='t';  k++; sfnoov[k]='_'; k++;
sfnoov[k]='o';  k++; sfnoov[k]='p';  k++; sfnoov[k]='e';  k++; sfnoov[k]='n';  k++; sfnoov[k]='!';  k++; sfnoov[k]='\0'; 
strcpy(szfov,ssdov); strcat(szfov,szfovk); k=strlen(szfov)+1; szfov[k]='\0';
if (vyve==1) { skpovvk=new char[dsov]; sdvovvk=new char[dsov]; sppovvk=new char[dsov]; 
skpovv=new char[2*dsov]; sdvovv=new char[2*dsov]; sppovv=new char[2*dsov]; 
if ((!skpovvk) || (!sdvovvk) || (!sppovvk) || (!skpovv) || (!sdvovv) || (!sppovv)) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<(2*dsov); k++) { skpovv[k]='\0'; sdvovv[k]='\0'; sppovv[k]='\0'; } 
for (k=0; k<dsov; k++) { skpovvk[k]='\0'; sdvovvk[k]='\0'; sppovvk[k]='\0'; }
k=0;
skpovvk[k]='K'; k++; skpovvk[k]='o'; k++; skpovvk[k]='e'; k++; skpovvk[k]='f'; k++; skpovvk[k]='f'; k++; skpovvk[k]='i'; k++;
skpovvk[k]='c'; k++; skpovvk[k]='i'; k++; skpovvk[k]='e'; k++; skpovvk[k]='n'; k++; skpovvk[k]='t'; k++; skpovvk[k]='_'; k++;
skpovvk[k]='p'; k++; skpovvk[k]='o'; k++; skpovvk[k]='g'; k++; skpovvk[k]='l'; k++; skpovvk[k]='o'; k++; skpovvk[k]='s'; k++;
skpovvk[k]='c'; k++; skpovvk[k]='h'; k++; skpovvk[k]='e'; k++; skpovvk[k]='n'; k++; skpovvk[k]='i'; k++; skpovvk[k]='y'; k++;
skpovvk[k]='a'; k++; skpovvk[k]='_'; k++; skpovvk[k]='v'; k++; skpovvk[k]='e'; k++; skpovvk[k]='r'; k++; skpovvk[k]='.'; k++;
skpovvk[k]='t'; k++; skpovvk[k]='x'; k++; skpovvk[k]='t'; k++; skpovvk[k]='\0';
k=0;
sdvovvk[k]='D'; k++; sdvovvk[k]='l'; k++; sdvovvk[k]='i'; k++; sdvovvk[k]='n'; k++; sdvovvk[k]='y'; k++; sdvovvk[k]='_'; k++; 
sdvovvk[k]='v'; k++; sdvovvk[k]='o'; k++; sdvovvk[k]='l'; k++; sdvovvk[k]='n'; k++; sdvovvk[k]='_'; k++; sdvovvk[k]='v'; k++; 
sdvovvk[k]='e'; k++; sdvovvk[k]='r'; k++; sdvovvk[k]='.'; k++; sdvovvk[k]='t'; k++; sdvovvk[k]='x'; k++; sdvovvk[k]='t'; k++;
sdvovvk[k]='\0'; 
k=0;
sppovvk[k]='P'; k++; sppovvk[k]='o'; k++; sppovvk[k]='k'; k++; sppovvk[k]='a'; k++; sppovvk[k]='z'; k++; sppovvk[k]='a'; k++;
sppovvk[k]='t'; k++; sppovvk[k]='_'; k++; sppovvk[k]='p'; k++; sppovvk[k]='r'; k++; sppovvk[k]='e'; k++; sppovvk[k]='l'; k++;
sppovvk[k]='o'; k++; sppovvk[k]='m'; k++; sppovvk[k]='l'; k++; sppovvk[k]='e'; k++; sppovvk[k]='n'; k++; sppovvk[k]='_'; k++;
sppovvk[k]='v'; k++; sppovvk[k]='e'; k++; sppovvk[k]='r'; k++; sppovvk[k]='.'; k++; sppovvk[k]='t'; k++; sppovvk[k]='x'; k++;
sppovvk[k]='t'; k++; sppovvk[k]='\0';
strcpy(sppovv,ssdov); strcat(sppovv,sppovvk); k=strlen(sppovv)+1; sppovv[k]='\0'; //cout << sppovv << endl;
strcpy(skpovv,ssdov); strcat(skpovv,skpovvk); k=strlen(skpovv)+1; skpovv[k]='\0'; //cout << skpovv << endl;
strcpy(sdvovv,ssdov); strcat(sdvovv,sdvovvk); k=strlen(sdvovv)+1; sdvovv[k]='\0'; //cout << sdvovv << endl;
}
}
double epsisredver(double T, double *tdkusctm, double *dkusctm, int dkusctl, double *dkoscet, double *dkoscem, int dkoscel)
{ int k=1; napstrdir(k,k); 
double *dv=NULL; dv=dliny_voln_ver(dv, 0);
dv=new double[dlar]; 
if (!dv) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<dlar; k++) dv[k]=0.0;
dv=dliny_voln_ver(dv, 1);
double *npp=new double[dlar];
if (!npp) { cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0; cout << "ce = " << dlar << endl;
npp=Kramers_Kronig_ver(dv, dlar, npp); cout << "k = " << k << endl;
double dkusct=opredKTPTKTochSha(dkusctm, tdkusctm, T, dkusctl); 
if (dkusct>1e0) dkusct=1e0; if (dkusct<0.0) dkusct=0.0; 
double dkosce=opredKTPTKTochSha(dkoscem, dkoscet, T, dkoscel); 
if (dkosce>1e0) dkosce=1e0; if (dkosce<0.0) dkosce=0.0; 
double ume=dkusct*dkosce; 
double *epsil=new double[dlar];
if (!epsil) { cout << snmov << endl; k=getchar(); exit(1); } 
epsil=epsilnu(npp, dlar, epsil);
double epssr=usredVelichPlank(dv, epsil, npp, T, dlar, ume); cout << "T = " << T << "\teps_s = " << epssr << "\tume = " << ume << "\t";
if (dv) delete []dv; if (npp) delete []npp; if (epsil) delete []epsil;
k=1; osvpamov(k); return epssr; }
double *dliny_voln_ver(double *dv, int ide)
{ ifstream fin; double p; int k, q=100, lear; char *s=new char[q]; 
if (!s) { cout << snmov << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
fin.open(sdvovv,ios_base::in); 
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); k++; } 
fin.clear(); fin.seekg(0);
lear=k; dlar=lear; 
if (!ide) { delete []s; fin.close(); return NULL; } else {
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); dv[k]=p;  k++; }
fin.close();
for (k=0; k<lear; k++) dv[k]=(1e-2)/dv[k];
delete []s; return dv; } }
double *Kramers_Kronig_ver(double *dv, int lear, double *nnu)
{	int k=0; ofstream fo;
fo.open(sppovv,ios_base::out | ios_base::trunc); 
if (!fo.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
	double c0=299792458.0, *nu=new double[lear], *nus=new double[lear], *nnus=new double[lear];
	double *alsr=new double[lear], e=1e-10; 
	if ((!nu) || (!nnus) || (!nus) || (!alsr)) { cout << snmov << endl; k=getchar(); exit(1); } 
	alsr=Koef_Pogl_ver(alsr); 
	for (k=0; k<lear; k++) if (dv[k]>e) nu[k]=2e0*pi*c0/dv[k]; else nu[k]=0.0; 
nus=izmMasChast(nu, enu, lear, nus); 
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear, nnus); 
nnu=PokazPrelomPribl(nu, nus, nnus, lear, nnu); 
if (nu) delete[]nu; if (nnus) delete[]nnus; if (nus) delete[]nus; if (alsr) delete[]alsr;
for (k=0; k<lear; k++) { nnu[k]=nnu[k]-1e0+1.543; fo << nnu[k] << endl; } fo.close(); 
return nnu; }
double opredKTPTKTochSha(double *ktptks, double *te, double temp, int ce)
{
	int n = ce, f = 1, p = 0, k; 
	double e=1e-4, ktp = 0.0, ko = 0.0, x1=0.0, x2=0.0, y1=0.0, y2=0.0, b=0.0, dt=0.0; 
	if ((temp>=te[0]) && (temp<=te[n-1])) {
	for (k = 0; k<n; k++)
		if ((te[k] >= temp) && (f>0)) { p = k; f = 0; break; } 
	if (!p) { if (!f) p = 1; else { p = n - 1; f = 0; }	} }
	else if (temp<te[0]) { p=1; f=0; }
	else if (temp>te[n-1]) { p=n-1; f=0; }
	if ((!f) && (p>0)) {
		x2=te[p]; x1=te[p - 1]; dt = x2 - x1;
		if (fabs(dt) > e) {
			y2=ktptks[p]; y1=ktptks[p - 1]; b=y1;
			if ((n-1)==p) b=y2;
			ko = (y2 - y1) / dt;
			ktp = b + ko*(temp - x1); } }
return ktp; }
void osvpamov(int vyve)
{ delete []sndov; delete []ssdov; delete []szfov; 
delete []szfovk; delete []snmov; delete []sfnoov; 
if (vyve==1) { 
delete []sppovv; delete []sppovvk; delete []sdvovv; 
delete []sdvovvk; delete []skpovv; delete []skpovvk; }
}
double *izmMasChast(double *nu, double epnu, int chel, double *nus)
{ int k, p=chel-1;
for (k=0; k<p; k++) nus[k]=epnu*(nu[k+1]-nu[k])+nu[k];
nus[p]=epnu*(nu[p-1]-nu[p-2])+nu[p-2];
return nus; }
double *PokazPrelomPribl(double *x, double *xs, double *ys, int lear, double *y)
{ double ko=0.0; int k, j, p=lear-1;
for (k=0; k<p; k++) {	
	ko=(ys[k+1]-ys[k])/(xs[k+1]-xs[k]); 
	y[k]=ko*(x[k]-xs[k])+ys[k]; }
k=p-1;
ko=(ys[k]-ys[k-1])/(xs[k]-xs[k-1]);
y[k]=(x[k]-xs[k-1])*ko+ys[k];
ko=(ys[p]-ys[p-1])/(xs[p]-xs[p-1]);
y[p]=(x[p]-xs[p-1])*ko+ys[p]; //for (k=0; k<lear; k++) { cout << "k = " << k << "\tnus = " << xs[k] << "\tns = " << ys[k] << "\tnu = " << x[k] << "\tn = " << y[k] << endl; } k=getchar();
return y; }
double *Koef_Pogl_ver(double *kp)
{	ifstream fin; double p; int k=0, q=100; char *s=new char[q]; 
	if (!s) { cout << snmov << endl; k=getchar(); exit(1); }
fin.open(skpovv,ios_base::in); for (k=0; k<q; k++) s[k]='\0'; //cout << "k = " << k << endl;
if (!fin.is_open()) { cout << sfnoov << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p; k++; }
fin.close(); delete []s; return kp; }
double usredVelichPlank(double *dv, double *uv, double *npp, double tem, int n, double umensh)
{ double PP=6.6260755e-34, PB=1.380658e-23, c0=299792458.0, vl, c1, c2, lambda; 
int k; double *Ib, *Ibn, *dl, nc, nz, e=1e-20;
Ib=new double[n]; Ibn=new double[n]; dl=new double[n];
if ((!Ib) || (!Ibn) || (!dl)) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<n; k++) {
    vl=c0/npp[k];
    c1=PP*pow(vl,2e0);
    c2=PP*vl/PB;
    lambda=dv[k]/npp[k];
    Ib[k]=2e0*pi*c1/(pow(lambda,5e0)*(exp(c2/lambda/tem)-1e0));
Ibn[k]=uv[k]*Ib[k];
dl[k]=lambda; }
nc=trapz(dl,Ibn,n);
nz=trapz(dl,Ib,n);
if (fabs(nz)>e) nc=umensh*nc/nz; else nc=0.0;
delete []Ibn; delete []Ib; delete []dl;
return nc; }
void initarrver(int koel, double wmg, double wsi, double wal)
{
	int dkoscvl = 6, cemv = 11, k=0, j=0;
	double *dkoscvm = new double[dkoscvl], *dkoscvt = new double[dkoscvl];
	if ((!dkoscvm) || (!dkoscvt)) { cout << snmv << endl; k = getchar(); exit(1); }
	double tnd = 6e2, dtd = 2e2, tm; 
	k=0; dkoscvt[k] = tnd; 
	for (k = 1; k < dkoscvl; k++) dkoscvt[k] = dkoscvt[k - 1] + dtd;
	k = 0; 
	dkoscvm[k] = 4.68; k++; dkoscvm[k] = 4.69; k++; 
	dkoscvm[k] = 5.65; k++; dkoscvm[k] = 13.17; k++; 
	dkoscvm[k] = 20.2; k++; dkoscvm[k] = 27.81;
	for (k = 0; k < dkoscvl; k++) {
		tm = dkoscvm[k] / 1e2; dkoscvm[k] = 1e0 - tm;
	} 
	double *tkuscv = new double[dmkooscv], *kuscv = new double[dmkooscv], g;
	double *etev=new double[cemv];
	if ((!tkuscv) || (!kuscv) || (!etev)) { cout << snmv << endl; k = getchar(); exit(1); }
	k=0; tkuscv[k] = tnoscv; 
	for (k = 1; k < dmkooscv; k++) tkuscv[k] = tkuscv[k - 1] + dtoscv;
	kuscv = koefoslab(wmg, wsi, wal, tkuscv, dmkooscv); 
	for (k=0; k<dmkooscv; k++) cout << "te = " << tkuscv[k] << "\tku = " << kuscv[k] << endl;
	double *stchsrver = new double[cemv];
	if (!stchsrver) { cout << snmv << endl; k = getchar(); exit(1); }
	etev[0] = tna; for (k = 1; k<cemv; k++) etev[k] = etev[k - 1] + dtv;
	for (k = 0; k<cemv; k++) { g = epsisredver(etev[k], tkuscv, kuscv, dmkooscv, dkoscvt, dkoscvm, dkoscvl); stchsrver[k] = g; }
}
double *koefoslab(double wmg, double wsi, double wal, double *tere, int n)
{ int lt=n, k; 
double wo=wmg+wsi+wal, kmg=0.0, kal=0.0, ksi=0.0, *kuo=new double[lt], ht=1e0, eps=1e-6;
double *mgo=new double[lt], *alo=new double[lt], *sio=new double[lt];
if ((!mgo) || (!alo) || (!sio) || (!kuo)) { cout << "No memory" << endl; k=getchar(); exit(1); }
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
if ((!mas) || (!st) || (!bs)) { cout << snmov << endl; k=getchar(); exit(1);}
st[0]=a11; st[1]=a12; mas[0]=st; st=new double[l];
if (st) { st[0]=a21; st[1]=a22; mas[1]=st; } else { cout << snmov << endl; k=getchar(); exit(1);}
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
double d, *m=new double[n], mm, **mmat; if (!m) { cout << snmv << endl; j=getchar(); exit(1); } //PrintMatr(mat, rm);
d=vychopred(mat,n); for (j=0; j<n; j++) { mmat=polmat(mat,ssc,j,n); mm=vychopred(mmat,n); 
m[j]=mm/d; /*cout << "x ( " << j << " ) = " << m[j] << endl;*/
mmat=osvpam(mmat,n); } 
return m; }
double **osvpam(double **a, int n)
{ int j; double *b; for (j=0; j<n; j++) { b=a[j]; delete []b; } return NULL; }
double **polmat(double **mat, double *b, int no, int n)
{ int k, j; double **ma=new double*[n], *m; 
if (!ma) { cout << snmv << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) {
	m=new double[n]; if (!m) { cout << snmv << endl; k=getchar(); exit(1); }
	for (j=0; j<n; j++) {
		if (j==no) m[j]=b[k]; else m[j]=mat[k][j];}
ma[k]=m; } 
return (ma); }
double provresh(double **A, double *x, double *xp, double *b, int n, double toch)
{ int k, j, q=0; double s, *pr=new double[n], m, t, mp;
if (!pr) { cout << snmv << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { s=0.0; 
for (j=0; j<n; j++) s=s+A[k][j]*x[j]; pr[k]=s-b[k]; 
} //cout << endl; 
m=0.0; for (j=0; j<n; j++) { t=fabs(pr[k]); if (t>m) { m=t; q=k; } } 
mp=0.0; for (k=0; k<n; k++) { t=fabs(xp[k]-x[k]); if (t>mp) { mp=t; } }
delete []pr; 
return m; } 
double vychopred(double  **mas, int m) {
  double d=Determinant(mas, m); //Вычисление определителя
  return d; }
double Determinant(double **mas, int m) { // Рекурсивное вычисление определителя
  int i;
  double **p, k=1., d=0.0, *po;
  p=new double*[m]; if (!p) { cout << snmv << endl; i=getchar(); exit(1); }
  for (i=0; i<m; i++) { po=new double[m]; if (!po) { cout << snmv << endl; i=getchar(); exit(1); } p[i]=po; }
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
double trapz(double *arx, double *ary, int lenarr)
{ int k; double s=0.0;
for (k=1; k<lenarr; k++) s=s+(ary[k]+ary[k-1])*(arx[k]-arx[k-1])/2e0;
return s; }
double *PokazPrelomAl(double *nu, double *nus, double *ka, double vl, int arle, double *np)
{ int k, j, p=arle-1, q; double dkpdo, podln, *fn=new double[arle]; 
if (!fn) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<arle; k++) {
    q=0; for (j=0; j<p; j++) {
        dkpdo=(ka[j+1]-ka[j])/(nu[j+1]-nu[j]);
        podln=((nu[j]+nus[k])/(nu[j]-nus[k]));
        podln=fabs(podln);
        fn[q]=dkpdo*log(podln);
        q++; }
    fn[p]=fn[p-1];
    np[k]=1.0+(vl/pi)*trapz(nu,fn,arle)/2.0/nu[k]; }
delete []fn; return np; }
double *epsilnu(double *npp, int lear, double *epsarr)
{ int k; double eps, n;
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
	eps=(4.0*n+2.0)/3.0/pow((n+1.0),2.0);
	eps=eps+2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-1.0)/(pow(n,2.0)+1.0)/(pow(n,4.0)-1.0);
	eps=eps-8.0*pow(n,4.0)*(pow(n,4.0)+1.0)*log(n)/(pow(n,2.0)+1.0)/pow((pow(n,4.0)-1.0),2.0);
	eps=eps-pow(n,2.0)*log((n-1.0)/(n+1.0))*pow((pow(n,2.0)-1.0),2.0)/pow((pow(n,2.0)+1.0),3.0); 
	epsarr[k]=eps; }
return epsarr; }