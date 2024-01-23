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
double epsisred(double, double *, double *, int, double *, double *, int, int, int, double);
double *izmMasChast(double *, double, int, double *);
double trapz(double *, double *, int);
double *Kramers_Kronig(double *, int, double *, char *, char *, char *, int, int, int, char *, double);
double **dliny_voln(char *, char *, char *);
double *epsilnu(double *, int, double *);
double *Koef_Pogl(double *, char *, char *, char *, int, double);
double *PokazPrelomAl(double *, double *, double *, double, int, double *, char *);
double *PokazPrelomPribl(double *, double *, double *, int, double *, double);
double nsredPlank2(double *, double *, double, int, char *);
double usredVelichPlank(double *, double *, double *, double, int, double, char *);
double opredKTPTKToch(double *, double *, double, int);
char **napolStrok(int, int);
char **napNazvFile(int, int, char *);
double koefPoglSred(double, int, int, double);
//-------------------
char **napolStrok(int vyve, int vm)
{	
	int k=0, dsov=60, n=2; 
	char *snmov=new char[dsov], *sfnoov=new char[dsov], **ms=new char*[n];
if ((!snmov) || (!sfnoov) || (!ms)) { cout << "No memory!" << endl; k=getchar(); exit(1); } 
for (k=0; k<dsov; k++) { snmov[k]='\0'; sfnoov[k]='\0'; } 
k=0;
snmov[k]='N'; k++; snmov[k]='o'; k++; snmov[k]='_'; k++; snmov[k]='m'; k++; snmov[k]='e'; k++; snmov[k]='m'; k++; 
snmov[k]='o'; k++; snmov[k]='r'; k++; snmov[k]='y'; k++; snmov[k]='!'; k++; snmov[k]='\0';
k=0;
sfnoov[k]='F';  k++; sfnoov[k]='i';  k++; sfnoov[k]='l';  k++; sfnoov[k]='e';  k++; sfnoov[k]='_';  k++; sfnoov[k]='i'; k++;
sfnoov[k]='s';  k++; sfnoov[k]='_';  k++; sfnoov[k]='n';  k++; sfnoov[k]='o';  k++; sfnoov[k]='t';  k++; sfnoov[k]='_'; k++;
sfnoov[k]='o';  k++; sfnoov[k]='p';  k++; sfnoov[k]='e';  k++; sfnoov[k]='n';  k++; sfnoov[k]='!';  k++; sfnoov[k]='\0'; 
k=0; ms[k]=snmov; k++; ms[k]=sfnoov; 
return ms;
}
char **napNazvFile(int vyve, int vm, char *snm)
{
	int k=0, dsov=60, j=0, q=0, f=4, m=0; 
	char *skpovk=new char[dsov], *sdvovk=new char[dsov], *sppovk=new char[dsov], *skpovt=new char[dsov];
if ((!skpovk) || (!sdvovk) || (!sppovk) || (!skpovt)) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<dsov; k++) { skpovk[k]='\0'; sdvovk[k]='\0'; sppovk[k]='\0'; }
k=0;
skpovk[k]='K'; k++; skpovk[k]='o'; k++; skpovk[k]='e'; k++; skpovk[k]='f'; k++; skpovk[k]='f'; k++; skpovk[k]='i'; k++;
skpovk[k]='c'; k++; skpovk[k]='i'; k++; skpovk[k]='e'; k++; skpovk[k]='n'; k++; skpovk[k]='t'; k++; skpovk[k]='_'; k++;
skpovk[k]='p'; k++; skpovk[k]='o'; k++; skpovk[k]='g'; k++; skpovk[k]='l'; k++; skpovk[k]='o'; k++; skpovk[k]='s'; k++;
skpovk[k]='c'; k++; skpovk[k]='h'; k++; skpovk[k]='e'; k++; skpovk[k]='n'; k++; skpovk[k]='i'; k++; skpovk[k]='y'; k++;
skpovk[k]='a'; k++; skpovk[k]='_'; k++; 
if (!vyve) { skpovk[k]='s'; k++; skpovk[k]='h'; k++; skpovk[k]='a'; k++; }
if (vyve==1) { skpovk[k]='v'; k++; skpovk[k]='e'; k++; skpovk[k]='r'; k++; }
if (vyve==2) { skpovk[k]='i'; k++; skpovk[k]='t'; k++; skpovk[k]='o'; k++; skpovk[k]='m'; k++; 
	if (!vm) { skpovk[k]='4'; k++; skpovk[k]='4'; k++; }
	if (vm==1) { skpovk[k]='6'; k++; skpovk[k]='2'; k++; }
	if (vm==2) { skpovk[k]='8'; k++; skpovk[k]='6'; k++; }
	if (vm==3) { skpovk[k]='1'; k++; skpovk[k]='0'; k++; skpovk[k]='0'; k++; }
skpovk[k]='0'; k++; } 
if (vyve==3) { skpovk[k]='k'; k++; skpovk[k]='v'; k++; skpovk[k]='i'; k++; 
	 if (vm==4) { skpovk[k]='4'; k++; }
	 if (vm==5) { skpovk[k]='5'; k++; }
	 if (vm==6) { skpovk[k]='6'; k++; }
	 if (vm==7) { skpovk[k]='7'; k++; }
	 if (vm==8) { skpovk[k]='8'; k++; }
	 if (vm==9) { skpovk[k]='9'; k++; }
	 if (vm==10) { skpovk[k]='1'; k++; skpovk[k]='0'; k++; }
	 if ((vm<=2) || (vm>10)) { cout << "Net takoy marki KVI!" << endl; k=getchar(); exit(1); }
	 skpovk[k]='0'; k++; skpovk[k]='0'; k++; }
for (j=0; j<k; j++) skpovt[j]=skpovk[j];
j=k; q=j; skpovt[j]='_'; j++; skpovt[j]='T'; j++; m=j-k;
skpovk[k]='.'; k++; skpovk[k]='t'; k++; skpovk[k]='x'; k++; skpovk[k]='t'; k++; skpovk[k]='\0'; k++; q=k-q+j;
for (k=j; k<q; k++) skpovt[k]=skpovk[k-m]; //cout << "skpovt = " << skpovt << "\tskpovk = " << skpovk << "\t"; k=getchar();
k=0;
sdvovk[k]='D'; k++; sdvovk[k]='l'; k++; sdvovk[k]='i'; k++; sdvovk[k]='n'; k++; sdvovk[k]='a'; k++; sdvovk[k]='_'; k++; 
sdvovk[k]='v'; k++; sdvovk[k]='o'; k++; sdvovk[k]='l'; k++; sdvovk[k]='n'; k++; sdvovk[k]='y'; k++; sdvovk[k]='_'; k++; 
if (!vyve) { sdvovk[k]='s'; k++; sdvovk[k]='h'; k++; sdvovk[k]='a'; k++; }
else if (vyve==1) { sdvovk[k]='v'; k++; sdvovk[k]='e'; k++; sdvovk[k]='r'; k++; }
else if (vyve==2) { sdvovk[k]='i'; k++; sdvovk[k]='t'; k++; sdvovk[k]='o'; k++; sdvovk[k]='m'; k++; }
else if (vyve==3) { sdvovk[k]='k'; k++; sdvovk[k]='v'; k++; sdvovk[k]='i'; k++; }
sdvovk[k]='.'; k++; sdvovk[k]='t'; k++; sdvovk[k]='x'; k++; sdvovk[k]='t'; k++; sdvovk[k]='\0';
k=0;
sppovk[k]='P'; k++; sppovk[k]='o'; k++; sppovk[k]='k'; k++; sppovk[k]='a'; k++; sppovk[k]='z'; k++; sppovk[k]='a'; k++; 
sppovk[k]='t'; k++; sppovk[k]='e'; k++; sppovk[k]='l'; k++; sppovk[k]='_'; k++; sppovk[k]='p'; k++; sppovk[k]='r'; k++; 
sppovk[k]='e'; k++; sppovk[k]='l'; k++; sppovk[k]='o'; k++; sppovk[k]='m'; k++; sppovk[k]='l'; k++; sppovk[k]='e'; k++; 
sppovk[k]='n'; k++; sppovk[k]='i'; k++; sppovk[k]='y'; k++; sppovk[k]='a'; k++; sppovk[k]='_'; k++;
if (!vyve) { sppovk[k]='s'; k++; sppovk[k]='h'; k++; sppovk[k]='a'; k++; }
else if (vyve==1) { sppovk[k]='v'; k++; sppovk[k]='e'; k++; sppovk[k]='r'; k++; }
else if (vyve==2) { sppovk[k]='i'; k++; sppovk[k]='t'; k++; sppovk[k]='o'; k++; sppovk[k]='m'; k++; }
else if (vyve==3) { sppovk[k]='k'; k++; sppovk[k]='v'; k++; sppovk[k]='i'; k++; }
sppovk[k]='.'; k++; sppovk[k]='t'; k++; sppovk[k]='x'; k++; sppovk[k]='t'; k++; sppovk[k]='\0'; //0 - шамот, 1 - вермикулит, 2 - ИТОМ, 3 - КВИ
char **unau=new char*[f];
if (!unau) { cout << snm << endl; k=getchar(); exit(1); } 
k=0; unau[k]=sppovk; k++; unau[k]=skpovk; k++; unau[k]=sdvovk; k++; unau[k]=skpovt;
return unau;
}
double epsisred(double T, double *tdkusctm, double *dkusctm, int dkusctl, double *dkoscet, double *dkoscem, int dkoscel, 
	int vybves, int vybmar, double dko)
{ 
	int k=0, j=0, q=j+1, l=q+1, m=l+1, i=0, dlar=0; 
	char **mauk=napolStrok(vybves, vybmar);
	char *snm=mauk[j], *sfno=mauk[q]; 
	char **mu=napNazvFile(vybves, vybmar, snm);
	char *sppov=mu[j], *skpov=mu[q], *sdvov=mu[l], *skpovt=mu[m]; 
	double **unau=dliny_voln(snm, sdvov, sfno), hf=1e0;
	double *dv=unau[j], *ne=unau[q], e=1e-9, d=e, t=ne[j]; 
while (d<t) { d=d+hf; i++; } dlar=i; if (ne) delete[]ne; if (unau) delete[]unau;
double *npp=new double[dlar]; if (!npp) { cout << snm << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0; 
npp=Kramers_Kronig(dv, dlar, npp, sppov, sfno, snm, vybves, vybmar, dlar, skpov, dko); 
double dkusct=opredKTPTKToch(dkusctm, tdkusctm, T, dkusctl); 
if (dkusct>hf) dkusct=hf; if (dkusct<0.0) dkusct=0.0; 
double dkosce=opredKTPTKToch(dkoscem, dkoscet, T, dkoscel); 
if (dkosce>hf) dkosce=hf; if (dkosce<0.0) dkosce=0.0; 
double ume=dkusct*dkosce; 
double *epsil=new double[dlar];
if (!epsil) { cout << snm << endl; k=getchar(); exit(1); } 
epsil=epsilnu(npp, dlar, epsil); 
double epssr=usredVelichPlank(dv, epsil, npp, T, dlar, ume, snm); 
for (k=0; k<=q; k++) { snm=mauk[k]; if (snm) delete[]snm; } if (mauk) delete[]mauk; 
for (k=0; k<=m; k++) { snm=mu[k]; if (snm) delete[]snm; } if (mu) delete[]mu; 
if (dv) delete []dv; if (npp) delete []npp; if (epsil) delete []epsil; 
return epssr; }
double *Kramers_Kronig(double *dv, int lear, double *nnu, char *sppov, char *sfno, char *snm, int vv, int vm, int dlar, char *skpov, double dko)
{	
	int k=0, z=0; 
	double c0=299792458.0, *nu=new double[lear], *nus=new double[lear], *nnus=new double[lear], c=0.0;
	double *alsr=new double[lear], e=1e-15, ht=1e0, ko=1e-2, enu=1e-3; 
	if ((!nu) || (!nnus) || (!nus) || (!alsr)) { cout << snm << endl; k=getchar(); exit(1); } 
	alsr=Koef_Pogl(alsr, snm, skpov, sfno, dlar, dko); 
	for (k=0; k<lear; k++) if (dv[k]>e) nu[k]=2e0*pi*c0/dv[k]; else nu[k]=0.0; 
nus=izmMasChast(nu, enu, lear, nus); 
nnus=PokazPrelomAl(nu, nus, alsr, c0, lear, nnus, snm); 
nnu=PokazPrelomPribl(nu, nus, nnus, lear, nnu, enu); 
if (nu) delete[]nu; if (nnus) delete[]nnus; if (nus) delete[]nus; if (alsr) delete[]alsr;
double nsh=1.56, nver=1.543;
	if (!vv) c=nsh-ht;
	else if (vv==1) c=nver-ht;
	else if (vv==2) { 
		double w1=0.0, w2=0.0, rosha=219.0*1e1, rover=25.0*1e1, phi1=0.0, phi2=0.0, v1=0.0, v2=0.0;
			if (!vm) { w1=6e1*ko; w2=ht-w1; }
			if (vm==1) { w1=5e1*ko; w2=ht-w1; }
			if (vm==2) { w1=25.0*ko; w2=ht-w1; }
			if (vm==3) { w1=2e1*ko; w2=ht-w1; }
			if ((vm<0) || (vm>3)) { cout << "Oshibka v vybore marki ITOM!"; k=getchar(); exit(1); }
v1=w1/rover; v2=w2/rosha; 
phi1=v1/(v1+v2); phi2=v2/(v1+v2); 
c=(nver-ht)*phi1+(nsh-ht)*phi2; }
	else if (vv==3) { 
		double wv=0.0; 
			if (vm==4) wv=49e-2; 
			if (vm==5) wv=37e-2; 
			if (vm==6) wv=31e-2; 
			if (vm==7) wv=27e-2;
			if (vm==8) wv=23e-2; 
			if (vm==9) wv=22e-2; 
			if (vm==10) wv=2e-1; 
			if ((vm<4) || (vm>10)) { cout << "Oshibka v vybore marki KVI!"; k=getchar(); exit(1); }
				c=(nver-ht)*wv+(ht-wv)*(nsh-ht); }
	else { cout << "Oshibka v vybore veschestva!"; k=getchar(); exit(1); }
	for (k=0; k<lear; k++) nnu[k]=nnu[k]+c; 
	if (z) { 
ofstream fo;
fo.open(sppov, ios_base::out | ios_base::trunc); 
if (!fo.is_open()) { cout << sfno << endl; k=getchar(); exit(1); }
for (k=0; k<lear; k++) fo << nnu[k] << endl; 
fo.close(); }
return nnu; }
double *epsilnu(double *npp, int lear, double *epsarr)
{ int k=0; double eps=0.0, n=0.0, hf=1e0;
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
	eps=(4.0*n+2.0)/3.0/pow((n+hf),2.0);
	eps=eps+2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-hf)/(pow(n,2.0)+hf)/(pow(n,4.0)-hf);
	eps=eps-8.0*pow(n,4.0)*(pow(n,4.0)+hf)*log(n)/(pow(n,2.0)+hf)/pow((pow(n,4.0)-hf),2.0);
	eps=eps-pow(n,2.0)*log((n-hf)/(n+hf))*pow((pow(n,2.0)-hf),2.0)/pow((pow(n,2.0)+hf),3.0); 
	epsarr[k]=eps; }
return epsarr; }
double *PokazPrelomAl(double *nu, double *nus, double *ka, double vl, int arle, double *np, char *snm)
{ 
	int k=0, j=0, p=arle-1, q=0; 
double dkpdo=0.0, podln=0.0, *fn=new double[arle], hf=1e0; 
if (!fn) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<arle; k++) {
    q=0; for (j=0; j<p; j++) {
        dkpdo=(ka[j+1]-ka[j])/(nu[j+1]-nu[j]);
        podln=((nu[j]+nus[k])/(nu[j]-nus[k]));
        podln=fabs(podln);
        fn[q]=dkpdo*log(podln);
        q++; }
    fn[p]=fn[p-1];
    np[k]=hf+(vl/pi)*trapz(nu,fn,arle)/2.0/nu[k]; }
if (fn) delete []fn; 
return np; }
double *izmMasChast(double *nu, double epnu, int chel, double *nus)
{ 
	int k, p=chel-1;
for (k=0; k<p; k++) nus[k]=epnu*(nu[k+1]-nu[k])+nu[k];
nus[p]=epnu*(nu[p-1]-nu[p-2])+nu[p-1];
return nus; }
double *PokazPrelomPribl(double *x, double *xs, double *ys, int arle, double *y, double e)
{ int k, p=arle-1; double ko=0.0;
for (k=0; k<arle; k++) y[k]=0.0;
for (k=0; k<(p-1); k++) {
	ko=(ys[k+1]-ys[k])/(xs[k+1]-xs[k]); 
	y[k]=ko*(x[k]-xs[k])+ys[k]; }
k=p-1;
ko=(ys[k]-ys[k-1])/(xs[k]-xs[k-1]);
y[k]=(x[k]-xs[k-1])*ko+ys[k];
y[p]=(ys[p]-ys[k])*e+ys[k];
return y; }
double *Koef_Pogl(double *kp, char *snm, char *skpov, char *sfno, int dlar, double dko)
{
	ifstream fin; 
	double p=0.0; int k=0, q=100; char *s=new char[q]; 
if (!s) { cout << snm << endl; k=getchar(); exit(1); }
fin.open(skpov, ios_base::in); for (k=0; k<q; k++) s[k]='\0'; 
if (!fin.is_open()) { cout << sfno << endl; k=getchar(); exit(1); }
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); kp[k]=p*dko; k++; }
fin.close(); if (s) delete []s; 
return kp; }
double **dliny_voln(char *snm, char *sdvov, char *sfno)
{
	ifstream fin; 
int k=1, q=100, lear=2, dlar=lear; char *s=new char[q]; 
double p=0.0, ht=1e-2, d=0.0, *dv=NULL, **mu=new double*[lear], hf=1e0, *nf=new double[k];
if ((!s) || (!mu) || (!nf)) { cout << snm << endl; k=getchar(); exit(1); } for (k=0; k<q; k++) s[k]='\0';
fin.open(sdvov, ios_base::in); 
if (!fin.is_open()) { cout << sfno << endl; k=getchar(); exit(1); }
k=0; while (!fin.eof()) { fin.getline(s,q,'\n'); p=atof(s); k++; } 
fin.clear(); fin.seekg(0); //cout << "sdvov = " << sdvov << "\t";
lear=k; dlar=lear; dv=new double[dlar]; 
if (!dv) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<dlar; k++) dv[k]=0.0; 
k=0; while ((!fin.eof()) && (k<dlar)) { fin.getline(s,q,'\n'); p=atof(s); dv[k]=p; k++; }
fin.close(); //cout << "dlar = " << dlar << "\t";
for (k=0; k<lear; k++) { dv[k]=ht/dv[k]; d=d+hf; }
if (s) delete[]s; 
k=0; mu[k]=dv; nf[k]=d; k++; mu[k]=nf;
return mu; } 
double trapz(double *arx, double *ary, int lenarr)
{ int k; double s=0.0;
for (k=1; k<lenarr; k++) s=s+(ary[k]+ary[k-1])*(arx[k]-arx[k-1])/2e0;
return s; }
double nsredPlank2(double *dv, double *npp, double tem, int n, char *snmov)
{ int k; double *npp2=new double[n], t, ume=1e0; 
if (!npp2) {cout << snmov << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) npp2[k]=pow(npp[k],2e0); 
t=usredVelichPlank(dv, npp2, npp, tem, n, ume, snmov);
if (npp2) delete []npp2; return t; }
double usredVelichPlank(double *dv, double *uv, double *npp, double tem, int n, double umensh, char *snmov)
{ 
	double PP=6.6260755e-34, PB=1.380658e-23, c0=299792458.0, vl=0.0, c1=0.0, c2=0.0, lambda=0.0; 
	int k=0; double *Ib=new double[n], *Ibn=new double[n], *dl=new double[n], nc=0.0, nz=0.0, e=1e-20;
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
if (Ibn) delete []Ibn; if (Ib) delete []Ib; if (dl) delete []dl;
return nc; }
double koefPoglSred(double T, int vybves, int vybmar, double dko)
{ 
	int k=0, j=0, q=j+1, l=q+1, i=0, dlar=0; 
	char **mauk=napolStrok(vybves, vybmar);
	char *snm=mauk[j], *sfno=mauk[q]; 
	char **mu=napNazvFile(vybves, vybmar, snm);
	char *sppov=mu[j], *skpov=mu[q], *sdvov=mu[l]; 
	double **unau=dliny_voln(snm, sdvov, sfno), hf=1e0, ume=hf;
	double *dv=unau[j], *ne=unau[q], e=1e-9, d=e, t=ne[j]; 
while (d<t) { d=d+hf; i++; } dlar=i; if (ne) delete[]ne; if (unau) delete[]unau;
double *alpha=new double[dlar]; if (!alpha) { cout << snm << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) alpha[k]=0.0;
alpha=Koef_Pogl(alpha, snm, skpov, sfno, dlar, dko); 
double *npp=new double[dlar]; if (!npp) { cout << snm << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0; 
npp=Kramers_Kronig(dv, dlar, npp, sppov, sfno, snm, vybves, vybmar, dlar, skpov, dko); 
double alsr=usredVelichPlank(dv, alpha, npp, T, dlar, ume, snm); 
for (k=0; k<=q; k++) { snm=mauk[k]; if (snm) delete[]snm; } if (mauk) delete[]mauk; 
for (k=0; k<=l; k++) { snm=mu[k]; if (snm) delete[]snm; } if (mu) delete[]mu; 
if (dv) delete []dv; if (alpha) delete []alpha; if (npp) delete []npp;
return alsr; }