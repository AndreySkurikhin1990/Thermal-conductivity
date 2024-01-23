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
double *Kramers_Kronig(double *, int, double *, char *, char *, char *, int, int, int, char *, double);
double **dliny_voln(char *, char *, char *);
double *ronuVes(double *, int, char *);
double usredVelichPlank(double *, double *, double *, double, int, double, char *);
double PoiskReflPhi(double, double, double, double);
double ReflSred(double, int, int, double, int);
double PoiskTetaShtrNach(double, double, double, double);
double PraCha10(double, double, double, double);
double *PoiSpeKoeOtr(double *, int);
double provUgla(double);
char **napolStrok(int, int);
char **napNazvFile(int, int, char *);
double *opredKoefOtrEdin(int, int, char *, double, double *, int, int);
//-----
double *opredKoefOtrEdin(int vybves, int vybmar, char *snm, double dko, double *tem, int dmko, int vyb)
{
	int k=0; 
	double *ko=new double[dmko], hf=1e0;
	if (!ko) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<dmko; k++) { ko[k]=ReflSred(tem[k], vybmar, vybves, dko, vyb); //cout << "ko ( " << k << " ) = " << ko[k] << "\t"; 
	} return ko;
}
double *ronuVes(double *npp, int lear, char *snm)
{ int k=0; double *ronu=new double[lear], n=0.0, hf=1e0; if (!ronu) { cout << snm << endl; k=getchar(); exit(1); } 
for (k=0; k<lear; k++) {
	n=fabs(npp[k]);
    ronu[k]=0.5+((n-hf)*(3.0*n+hf))/(6.0*pow(n+hf,2.0))-(2.0*pow(n,3.0)*(pow(n,2.0)+2.0*n-hf))/((pow(n,2.0)+hf)*(pow(n,4.0)-hf)); 
    ronu[k]=ronu[k]+(8.0*pow(n,4.0)*(pow(n,4.0)+hf)*log(n))/((pow(n,2.0)+hf)*pow((pow(n,4.0)-hf),2.0));
    ronu[k]=ronu[k]+pow(n,2.0)*pow((pow(n,2.0)-hf),2.0)*log(fabs((n-hf)/(n+hf)))/pow((pow(n,2.0)+hf),3.0); }
return ronu; }
double PoiskReflPhi(double phi, double phis, double n1, double n2)
{ 
	double rpa=(n2*cos(phi)-n1*cos(phis))/(n2*cos(phi)+n1*cos(phis)); //коэффициент отражени€ параллельный плоскости падени€
	double rpe=(n1*cos(phi)-n2*cos(phis))/(n1*cos(phi)+n2*cos(phis)); //коэффициент отражени€ перпендикул€рный к плоскости падени€
return (fabs(pow(rpa,2e0))+fabs(pow(rpe,2e0)))/2e0; }
double ReflSred(double T, int vybmar, int vybves, double dko, int vyb)
{ 
	int k=0, j=0, q=j+1, l=q+1, m=l+1, i=0, dlar=0; 
	char **mauk=napolStrok(vybves, vybmar);
	char *snm=mauk[j], *sfno=mauk[q]; 
	char **mu=napNazvFile(vybves, vybmar, snm);
	char *sppov=mu[j], *skpov=mu[q], *sdvov=mu[l], *skpovt=mu[m]; 
	double **unau=dliny_voln(snm, sdvov, sfno), hf=1e0;
	double *dv=unau[j], *ne=unau[q], e=1e-9, d=e, t=ne[j]; 
	i=0; while (d<t) { d=d+hf; i++; } dlar=i; if (ne) delete[]ne; if (unau) delete[]unau; 
	double *npp=new double[dlar]; if (!npp) { cout << snm << endl; k=getchar(); exit(1); }
for (k=0; k<dlar; k++) npp[k]=0.0;
npp=Kramers_Kronig(dv, dlar, npp, sppov, sfno, snm, vybves, vybmar, dlar, skpov, dko);
double *ronu=NULL, ume=1e0;
if (!vyb) ronu=PoiSpeKoeOtr(npp, dlar); else ronu=ronuVes(npp, dlar, snm);
t=usredVelichPlank(dv, ronu, npp, T, dlar, ume, snm); 
for (k=0; k<=q; k++) { snm=mauk[k]; if (snm) delete[]snm; } if (mauk) delete[]mauk; 
for (k=0; k<=m; k++) { snm=mu[k]; if (snm) delete[]snm; } if (mu) delete[]mu; 
if (dv) delete []dv; if (npp) delete []npp; if (ronu) delete []ronu; 
return t; }
double *PoiSpeKoeOtr(double *npp, int lear)
{ 
	double hf=1e0, maxphi=pi/2e0, minphi=0.0, phi=0.0, enu=1e-3, hphi=enu;
	double n2=hf, n1=hf, *rs=new double[lear], p=hf, r=hf, psi=hf; int k=0;
for (k=0; k<lear; k++) { 
	phi=minphi; n2=npp[k]; p=0.0; r=0.0;
while (phi<maxphi) { 
	psi=PoiskTetaShtrNach(n2-n1, phi, n1, n2); 
	phi=phi+hphi; p=p+hf; 
	r=r+PoiskReflPhi(phi, psi, n1, n2); }
if (fabs(p)>0.0) rs[k]=r/p; } 
return rs; }
double PoiskTetaShtrNach(double dnr, double phipad, double n1, double n2)
{ double phipre=0.0, a=0.0, b=pi/2e0, c=a, ep=1e-8, ra=fabs(a-b);
double fa=0.0, fb=fa, fc=fa, phiprel=fa;
int Nit=10000, k=0;
if (dnr<0.0) { phipre=asin(n2/n1); phipre=provUgla(phipre); }
if (((phipad<phipre) && (dnr<0.0)) || (dnr>0.0)) {
while ((ra>ep) && (k<Nit)) {
    c=(a+b)/2e0;
    fa=PraCha10(phipad,a,n1,n2);
    fb=PraCha10(phipad,b,n1,n2);
    fc=PraCha10(phipad,c,n1,n2);
    if ((fc*fb>0.0) && (fa*fc<0.0)) b=c; 
	if ((fc*fa>0.0) && (fb*fc<0.0)) a=c; 
	k++; ra=fabs(a-b); }
	phiprel=provUgla(c); }
else phiprel=pi/2e0;
return phiprel; }
double PraCha10(double phi, double phis, double n1, double n2)
{ return n1*sin(phi)-n2*sin(phis); }
double provUgla(double ugol)
{ double ugo=ugol;
if (ugo<0.0) ugo=fabs(ugo);
if (fabs(ugo)>(pi/2e0)) ugo=pi-ugo;
return ugo; }