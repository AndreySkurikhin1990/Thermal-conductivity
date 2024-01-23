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
double trapz(double *, double *, int);
double *Kramers_Kronig(double *, int, double *, char *, char *, char *, int, int, int, char *, double);
double **dliny_voln(char *, char *, char *);
double *Koef_Pogl(double *, char *, char *, char *, int, double);
double KoefPoglRossel(double, double *, double *, double *, int, char *);
double *KoefPoglRosselNac(double *, int, double, double, int, int);
double nsredPlank2(double *, double *, double, int, char *);
double usredVelichPlank(double *, double *, double *, double, int, double, char *);
double LuchKTPChudnovsky(double *, double, int, double);
double opredKTPTKToch(double *, double *, double, int);
char **napolStrok(int, int);
char **napNazvFile(int, int, char *);
//-----
double KoefPoglRossel(double tem, double *dl, double *alfs, double *npp, int n, char *snmov)
{ 
	double PP=6.6260755e-34, PB=1.380658e-23, c0=299792458.0, C1=PP*pow(c0,2.0), C2=PP*c0/PB;
	double sig=2e0*C1*pow(pi,5e0)/15.0/pow(C2,4e0), np2=nsredPlank2(dl, npp, tem, n, snmov); 
	double t=0.0, chi=0.0, zna=0.0, chasc=0.0, chasz=0.0, rs=0.0, hf=1e0;
	double *Ibc=new double[n], *Ibz=new double[n], *dv=new double[n], dlv=0.0, c1m=0.0, c2m=0.0; int k=0;
if ((!Ibc) || (!Ibz) || (!dv)) { cout << snmov << endl; k=getchar(); exit(1); } 
for (k=0; k<n; k++) {
c1m=C1/pow(npp[k],2e0);
c2m=C2/npp[k];
t=(pi/2e0)*(c1m*c2m)/sig;
dlv=dl[k]/npp[k];
chi=exp(c2m/tem/dlv);
zna=pow((chi-hf),2e0);
Ibz[k]=t*chi/zna/pow(dlv,6e0)/pow(tem,5e0)/np2; 
Ibc[k]=Ibz[k]/alfs[k];
dv[k]=dlv; }
chasc=trapz(dv,Ibc,n);
chasz=trapz(dv,Ibz,n); 
rs=chasc/chasz; 
if (fabs(rs)>0.0) rs=hf/rs;
if (Ibc) delete []Ibc; if (Ibz) delete []Ibz; if (dv) delete []dv;
return rs; }
double *KoefPoglRosselNac(double *tem, int arle, double dkosups, double dkokpiko, int vybves, int vybmar)
{ //расчет лучистого КТП через КП по Росселанду
	int vm=vybmar, vyve=vybves; 
	int k=0, j=0, q=j+1, l=q+1, i=0, dlar=0;
	char **mauk=napolStrok(vybves, vybmar);
	char *snm=mauk[j], *sfno=mauk[q]; if (mauk) delete[]mauk; 
	char **mu=napNazvFile(vybves, vybmar, snm);
	char *sppov=mu[j], *skpov=mu[q], *sdvov=mu[l]; if (mu) delete[]mu;
	//-----
	double *alfs=NULL, *npp=NULL, *kpr=new double[arle], krpk=0.0, temk=0.0, hf=1e0;
	double *ktr=new double[arle], dkusct=hf, dkosce=hf; //КТП по Росселанду
	double dko2=dkosups, dko4=dkokpiko, dko=dko2*dko4; //учет пористой структуры и соотношения КП/КО 
	double **unau=dliny_voln(snm, sdvov, sfno);
	double *dl=unau[j], *ne=unau[q], e=1e-9, d=e, t=ne[j]; 
	while (d<t) { d=d+hf; i++; } dlar=i; if (ne) delete[]ne; if (unau) delete[]unau;
	if (dko<0.0) dko=0.0; if (dko>hf) dko=hf; 
	//-----
	alfs=new double[dlar]; if (!alfs) { cout << snm << endl; k=getchar(); exit(1); } 
	for (k=0; k<dlar; k++) alfs[k]=0.0; 
	alfs=Koef_Pogl(alfs, snm, skpov, sfno, dlar, dko); 
	npp=new double[dlar]; if (!npp) { cout << snm << endl; k=getchar(); exit(1); } 
	for (k=0; k<dlar; k++) npp[k]=0.0; 
	npp=Kramers_Kronig(dl, dlar, npp, sppov, sfno, snm, vybves, vybmar, dlar, skpov, dko); 
	double *GnpT=new double[arle], sigma=5.67e-8, *npT=new double[arle]; //производная dn/dT, функция <n>(T)
	if ((!kpr) || (!GnpT) || (!ktr) || (!npT)) { cout << snm << endl; k=getchar(); exit(1); }
	//-----
	for (k=0; k<arle; k++) {
	temk=tem[k]; GnpT[k]=0.0; ktr[k]=0.0; npT[k]=0.0; kpr[k]=0.0;
	npT[k]=usredVelichPlank(dl, npp, npp, temk, dlar, hf, snm); 
	krpk=KoefPoglRossel(temk, dl, alfs, npp, dlar, snm); kpr[k]=krpk; } 
	if (arle>1) { for (k=1; k<arle; k++) GnpT[k-1]=(npT[k]-npT[k-1])/(tem[k]-tem[k-1]); 
	k=arle-1; GnpT[k]=2e0*GnpT[k-1]-GnpT[k-2]; } 
	if (arle<0) { cout << "Oshibka pri vvode chisla elementov!" << endl; k=getchar(); exit(1); }
	for (k=0; k<arle; k++) {
	temk=tem[k]; krpk=kpr[k];
    krpk=fabs(8e0*npT[k]*sigma*pow(temk,3e0)/(3e0*krpk));
    krpk=krpk*fabs(2e0*npT[k]+temk*GnpT[k]); 
	ktr[k]=krpk; } 
	if (dl) delete []dl; if (npp) delete []npp; if (alfs) delete []alfs; 
	if (kpr) delete []kpr; if (GnpT) delete []GnpT; if (npT) delete []npT; 
	if (snm) delete[]snm; if (sfno) delete[]sfno; 
	if (sppov) delete[]sppov; if (skpov) delete[]skpov; if (sdvov) delete[]sdvov;
	return ktr; }
double LuchKTPChudnovsky(double *Ab, double tem, int kost, double razm)
{ 
	double s=0.0, t=s, hf=1e0, sig=5.67e-8, r=s; 
	int k=0; 
	for (k=0; k<kost; k++) { 
		s=s+Ab[k]; t=t+hf; } 
	if (fabs(t)>0.0) s=s/t;
	r=(3e0/4e0)*2e0*sig*pow(s,2e0)*pow(tem,3e0)*razm;
return r; }