#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <cmath>
#include <windows.h>
#include <iostream>
using namespace std;
const double pi = acos(-1e0), stpoSha33=33.0*1e-2, stpoSha24=24.0*1e-2, stpoSha20=2e1*1e-2;
const double stpoSha16=16.0*1e-2, stpoSha30=3e1*1e-2, stpoSha10=1e1*1e-2, tem0=273.15, prk=3e2, pri=16e1;
const int dsnftk=50, dmvtk=28, dmsrps=23;
const int dmsrpv=18, dmsrpi=17, dmrprs3=31, dmrprs4=32, dmrprs5=25, dmrprs6=33;
const int kv4=23, kv5=46, kv6=41, kv7=45, kv8=54, kv9=39, kv10=17, cvym=6;
struct raspporporazmm { double rprm; struct raspporporazmm *slel; struct raspporporazmm *prel; };
double *NapMasRaspSha(double *);
double *NapMasPrGrSha(double *, int);
double **poisrasprpor1(double, double, int, double, double, int, double *, double *, double *, double **);
double **poisrasprpor2(double, double, int, double, double, int, double *, double *, double *, double **);
double **rasppor2(double, double, double *, double *, double *, int, double, double, int, double **);
double *rasPorpoRazmVerI207(double *, int);
double *rasPorpoRazmVerIs84(double *, int);
double *rasPorpoRazmVerIn84(double *, int);
double *rasPorpoRazmVerI16035(double *, int);
double *rasPorpoRazmVerO84(double *, int);
double *rasPorpoRazmVerO16035(double *, int);
double **rasPorpoRazkvi(int);
double *NapMasRaspVer(int, int, int, int, double *);
double *NapMasRasLePritom(int, int, int, int);
double **NapMasRasLePrSha(int, double *, double *, double *, double *, double *, double **);
double **NapMasRasLePrSha2(int, double *, double *, double *, double *, double *, double **);
double **NapMasRasLePrSha3(int, double *, double *, double *, double *, double *, double **);
double **NapMasRasLePrSha4(int, double *, double *, double *, double *, double *, double **);
double **NapMasRasLePrSha5(int, double *, double *, double *, double *, double *, double **);
double **NapMasRasLePrSha6(int, double *, double *, double *, double *, double *, double **);
double **NapMasRasLePrVer(int, int, int, int, double **);
double **rasPorpoRazVer(double, int, int, int, int, int);
double **rasrasporporaz1(double *, double, int, int, double *, double *, double *, double, double, double, int, double **);
double **rasrasporporaz0(double *, double *, double *, double, int, int, double, double **);
double *rasPorpoRazmitom440(double *, int);
double *rasPorpoRazmitom620(double *, int);
double *rasPorpoRazmitom860(double *, int);
double *rasPorpoRazmitom1000(double *, int);
double *NapMasRaspitom(int, int, double *);
double *PravGranPoromitom(double *, int);
double *LevGranPoromitom(double *, double *, int);
double **rasPorpoRazitom(int);
double *NapMasRaspSha2(double *, int);
double *NapMasPrGrSha2(double *, int);
double **rasPorpoRazSha2(double, int);
double *NapMasRaspSha3(double *, int);
double *NapMasPrGrSha3(double *, int);
double **rasPorpoRazSha3(double, int);
double *NapMasRaspSha4(double *, int);
double *NapMasPrGrSha4(double *, int);
double **rasPorpoRazSha4(double, int);
double *NapMasRaspSha5(double *, int);
double *NapMasPrGrSha5(double *, int);
double **rasPorpoRazSha5(double, int);
double *NapMasRaspSha6(double *, int);
double *NapMasPrGrSha6(double *, int);
double **rasPorpoRazSha6(double, int);
double *rasprPorpoRazmAbskvi400();
double *LegrRasPorpoRazmkvi400();
double *rasprPorpoRazmAbskvi500();
double *LegrRasPorpoRazmkvi500();
double *rasprPorpoRazmAbskvi600(); 
double *LegrRasPorpoRazmkvi600();
double *rasprPorpoRazmAbskvi700(); 
double *LegrRasPorpoRazmkvi700();
double *rasprPorpoRazmAbskvi800(); 
double *LegrRasPorpoRazmkvi800();
double *rasprPorpoRazmAbskvi900(); 
double *LegrRasPorpoRazmkvi900();
double *rasprPorpoRazmAbskvi1000(); 
double *LegrRasPorpoRazmkvi1000();
double *PravGranPoromVer(double *, int);
double **obrabMaskvi(double *,double *, int, double **);
double *LevGranPoromVer(double *, double *, int);
double **poisrasprporkvi(int, double **);
double **rasrasporporazitom(double *, double *, double *, double, int, double, double **);
double **vybRasPorRazVer(double *, int, int, int, int);
double **vybRasPorRazSha(double, double, int, int, double *, double*, double *, double **);
double **vybFunRasPorpoRazSha(double, int, int);
//----------------//Задание распределения пор по размерам - ИТОМ
double *rasPorpoRazmitom440(double *rapo, int n)
{ int k=0; double s=0.0;
rapo[k]=0.57; k++; rapo[k]=0.3;   k++; rapo[k]=0.6;   k++; rapo[k]=0.87;  k++; rapo[k]=1.35; k++;
rapo[k]=2.07; k++; rapo[k]=3.72;  k++; rapo[k]=3.81;  k++; rapo[k]=5.38;  k++; rapo[k]=7.6;  k++;
rapo[k]=9.67; k++; rapo[k]=10.87; k++; rapo[k]=34.68; k++; rapo[k]=10.78; k++; rapo[k]=6.25; k++;
rapo[k]=1.47; k++; rapo[k]=0.0; k++;
for (k=0; k<n; k++) s=s+rapo[k]; for (k=0; k<n; k++) rapo[k]=rapo[k]*1e2/s;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *rasPorpoRazmitom620(double *rapo, int n)
{ int k=0; double s=0.0;
rapo[k]=0.26; k++; rapo[k]=0.15; k++; rapo[k]=0.26; k++; rapo[k]=0.22;  k++; rapo[k]=0.81; k++; 
rapo[k]=0.81; k++; rapo[k]=1.88; k++; rapo[k]=3.95; k++; rapo[k]=5.54;  k++; rapo[k]=7.35; k++;
rapo[k]=7.09; k++; rapo[k]=9.01; k++; rapo[k]=34.9; k++; rapo[k]=13.59; k++; rapo[k]=7.5;  k++;
rapo[k]=2.14; k++; rapo[k]=4.54;
for (k=0; k<n; k++) s=s+rapo[k]; for (k=0; k<n; k++) rapo[k]=rapo[k]*1e2/s;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *rasPorpoRazmitom860(double *rapo, int n)
{ int k=0; double s=0.0;
rapo[k]=0.4;   k++; rapo[k]=0.09;  k++; rapo[k]=0.44;  k++; rapo[k]=0.22;  k++;
rapo[k]=0.66;  k++; rapo[k]=1.02;  k++; rapo[k]=1.33;  k++; rapo[k]=2.66;  k++;
rapo[k]=4.07;  k++; rapo[k]=10.71; k++; rapo[k]=12.17; k++; rapo[k]=11.29; k++;
rapo[k]=35.06; k++; rapo[k]=11.24; k++; rapo[k]=7.13;  k++; rapo[k]=1.51;  k++;
rapo[k]=0.0;
for (k=0; k<n; k++) s=s+rapo[k]; for (k=0; k<n; k++) rapo[k]=rapo[k]*1e2/s;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *rasPorpoRazmitom1000(double *rapo, int n)
{ int k=0; double s=0.0;
rapo[k]=0.23;  k++; rapo[k]=0.19;  k++; rapo[k]=0.04;  k++; rapo[k]=0.61; k++;
rapo[k]=0.23;  k++; rapo[k]=1.03;  k++; rapo[k]=0.8;   k++; rapo[k]=2.47; k++;
rapo[k]=5.66;  k++; rapo[k]=10.87; k++; rapo[k]=14.18; k++; rapo[k]=12.5; k++;
rapo[k]=32.61; k++; rapo[k]=11.59; k++; rapo[k]=5.25;  k++; rapo[k]=1.75; k++;
rapo[k]=0.0;
for (k=0; k<n; k++) s=s+rapo[k]; for (k=0; k<n; k++) rapo[k]=rapo[k]*1e2/s;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *PravGranPoromitom(double *rapo, int n)
{ int k=0; 
rapo[k]=15e1;  k++; rapo[k]=95.0; k++; rapo[k]=85.0; k++; rapo[k]=75.0; k++; rapo[k]=65.0; k++;
rapo[k]=55.0;  k++; rapo[k]=45.0; k++; rapo[k]=35.0; k++; rapo[k]=25.0; k++; rapo[k]=15.0; k++;
rapo[k]=7.5;   k++; rapo[k]=4.0;  k++; rapo[k]=2.0;  k++; rapo[k]=0.75; k++; rapo[k]=0.3;  k++;
rapo[k]=0.055; k++; rapo[k]=0.0055;
for (k=0; k<n; k++) rapo[k]=(1e-6)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; } 
return rapo; }
double *LevGranPoromitom(double *rapo, double *levgr, int n)
{ int k; levgr[0]=0.0; for (k=1; k<n; k++) levgr[k]=rapo[k-1]; return levgr; }
double **rasPorpoRazitom(int vvi)
{ int n=dmsrpi, k=0, f=cvym, j=0, qg=0; //число выходных массивов
double rpn=1e0, dp=1e0, ht=dp, **mu=new double*[f], *rp=NULL, *mk=NULL, *ms=NULL, *legr0=NULL, srp=0.0;
double s=0.0, l=0.0, h=1e-6, p=h, *uv=NULL, *prgr0=NULL, *rpr0=NULL, *po=NULL, m=0.0, e=1e-1, marp=0.0;
prgr0=new double[n]; rpr0=new double[n]; legr0=new double[n];
if ((!rpr0) || (!prgr0) || (!mu) || (!legr0)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { rpr0[k]=0.0; prgr0[k]=0.0; legr0[k]=0.0; }
prgr0=PravGranPoromitom(prgr0, n); 
for (k=1; k<n; k++) legr0[k]=prgr0[k-1];
rpr0=NapMasRaspitom(vvi, n, rpr0);
mu=rasrasporporazitom(prgr0, legr0, rpr0, rpn, n, dp, mu);
if (prgr0) delete[]prgr0; if (rpr0) delete[]rpr0; if (legr0) delete[]legr0;
return mu; }
double **rasrasporporazitom(double *prgr0, double *legr0, double *raspr0, double rpn, int n, double de, double **mu)
{ double s=0.0, e=1e-6, r=0.0, ht=1e0, koef=1e6, *prgr00=new double[n], srpil, srpir, rprit;
double *legr00=new double[n], m=0.0, t=0.0, e1=ht+e, lf, prl=0.0, prr=0.0, legr01=0.0;
int k=2, q=0, p=0, j=0, jk=0, l=0, b=0, u=0, fl, flg, w=0, pd, ld, pflg=-1, qg;
for (k=0; k<n; k++) prgr00[k]=prgr0[k]*koef;
legr00[0]=0.0; for (k=1; k<n; k++) legr00[k]=prgr00[k-1]; //for (k=0; k<n; k++) cout << "le ( " << k << " ) = " << legr00[k] << "\tpr = " << prgr00[k] << endl;
for (k=0; k<n; k++) if (legr00[k]<e1) p=k; else break;
for (k=0; k<p; k++) s=s+raspr0[k];
if (legr00[p]<=ht) m=raspr0[p]*(ht-legr00[p])/(prgr00[p]-legr00[p]); //размер пор до 1 мкм
s=s+m; q=1; //cout << "p = " << p << "\ts = " << s << endl;
t=0.0; j=0; while (t<pri) { j++; t=t+ht; } jk=j; qg=jk;
double *ras01=new double[qg], *prgrm=new double[qg], *legrm=new double[qg];
double *srra=new double[qg], prgr01=rpn, srp, marp; 
if ((!ras01) || (!prgrm) || (!legrm)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
fl=1; t=0.0; for (j=0; j<qg; j++) { for (u=0; u<n; u++) if (fabs(legr00[u]-t)>e) { fl=-1; break; } if (fl<0) break; t=t+ht; } //fl - между целыми - факт существования
for (k=0; k<qg; k++) { prgrm[k]=0.0; legrm[k]=0.0; ras01[k]=0.0; } 
w=0; ras01[w]=s; prgrm[w]=prgr01; legrm[w]=legr01; legr01=prgr01; prgr01=prgr01+de; w++;
for (k=p; k<n; k++) { 
    srpil=legr00[k]; srpir=prgr00[k]; r=srpir-srpil; rprit=raspr0[k];
	ld=0; t=0.0; while (t<(srpil-e)) { ld++; t=t+ht; } pd=ld; while (t<(srpir-e)) { pd++; t=t+ht; } 
	t=srpil; j=0; while (t<srpir) { j++; t=t+ht; } jk=j; l=jk; flg=-1; 
	if (fl<0) { t=0.0; lf=t+ht; for (j=0; j<qg; j++) { if ((srpir>t) && (srpir<lf) && (fabs(t-srpir)>e)) 
	{ flg=1; pd--; jk--; l--; break; } else { t=lf; lf=lf+ht; } } } //cout << "ld = " << ld << "\tpd = " << pd << endl;
	m=rprit/r; //cout << "jk = " << jk << "\tk = " << k << "\tleft = " << srpil << "\tright = " << srpir << "\trprit = " << rprit << "\tprgr01 = " << prgr01 << endl;
	for (b=ld; b<pd; b++) { //до 2 мкм
		q++; s=s+m; ras01[w]=m; legrm[w]=legr01; prgrm[w]=prgr01; legr01=prgr01; prgr01=prgr01+de; w++; } 
	if ((flg>0) && ((k+1)<n)) {
		m=rprit*(srpir-t)/r; 
		rprit=raspr0[k+1]; srpil=legr00[k+1]; srpir=prgr00[k+1]; r=srpir-srpil;
		m=m+rprit*(lf-srpil)/r; 
		q++; s=s+m; prgrm[w]=prgr01; ras01[w]=m; legrm[w]=legr01; legr01=prgr01; prgr01=prgr01+de; w++; //cout << "l = " << l << "\tm = " << m << endl; 
	} pflg=flg; }
jk=15;
delete []prgr00; delete []legr00; //int *mm=new int[jk]; double *mf=new double[jk], *kl=new double[jk]; if ((!mm) || (!mf) || (!kl)) { cout << "No memory!" << endl; k=getchar(); exit(1); } for (k=0; k<jk; k++) { mm[k]=0; mf[k]=0.0; kl[k]=0.0; } u=1; q=u; j=0; r=ras01[j]; prl=legrm[j]; for (k=1; k<qg; k++) { m=ras01[k]; s=fabs(m-r); srpil=legrm[k]; if ((m>e) && (s>e)) { mm[j]=q; mf[j]=m; kl[j]=prl; j++; q=u; prl=srpil; } else q++; r=m; } //for (k=0; k<jk; k++) cout << "mm ( " << k << " ) = " << mm[k] << "\tmf = " << mf[k] << "\tle = " << kl[k] << endl; delete[]mf; delete[]mm; delete[]kl; cout << "qg = " << qg << "\ts = " << s << endl; //t=prgr01; u=k; cout << "u = " << u << endl; for (j=u; j<qg; j++) { prgrm[k]=prgr01; legrm[k]=prgrm[k-1]; prgr01=prgr01+de; } lf=0.0; for (k=0; k<u; k++) lf=lf+ras01[k]; cout << "s ras01 = " << lf << endl; 
s=0.0; for (k=0; k<n; k++) s=s+raspr0[k]; //cout << "s raspr0 = " << s << endl; //s=s-lf; s=s/(t-prgr01); for (j=u; j<qg; j++) ras01[j]=s; 
s=0.0; for (k=0; k<qg; k++) { //cout << "rp " << k << " = " << ras01[k] << "\tle = " << legrm[k] << "\tpr = " << prgrm[k] << endl; cout << endl << endl; 
s=s+ras01[k]; } //cout << "s ras01 = " << s << endl; 
double *ms=new double[1], *mm=new double[1];
if ((!ms) || (!mm)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
r=0.0; for (k=0; k<qg; k++) r=r+ht; k=0; mm[k]=r;
for (k=0; k<qg; k++) srra[k]=(legrm[k]+prgrm[k])/2e0;
k=0; mu[k]=ras01; k++; mu[k]=srra; k++; mu[k]=prgrm; k++; mu[k]=legrm; k++; mu[k]=ms; k++; mu[k]=mm;
return mu; }
double *NapMasRaspitom(int vvi, int n, double *r)
{ 
int k=0;
if (!vvi) r=rasPorpoRazmitom440(r, n); //ИТОМ-440
else if (vvi==1) r=rasPorpoRazmitom620(r, n); //ИТОМ-620
else if (vvi==2) r=rasPorpoRazmitom860(r, n); //ИТОМ-860
else if (vvi==3) r=rasPorpoRazmitom1000(r, n); //ИТОМ-1000
return r;
}

//Шамот - ШБ-1 №2 2-1
double **vybFunRasPorpoRazSha(double porsha, int vrsh, int vystsha)
{ double *po; int k=0, n=0, m=cvym;
double ko=1e-2, rpn=1e0, dp=1e0, p0=porsha, r=1e-6; 
double *legr=NULL, *prgr=NULL, *raspr=NULL;
double *srra=NULL, *rpr=NULL, **mu=new double*[m];
legr=new double[n]; prgr=new double[n]; raspr=new double[n]; 
srra=new double[n]; rpr=new double[n];
if ((!legr) || (!prgr) || (!raspr) || (!srra) || (!rpr) || (!mu)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { raspr[k]=0.0; prgr[k]=0.0; legr[k]=0.0; srra[k]=0.0; rpr[k]=0.0; } 
k=0; mu[k]=legr; k++; mu[k]=prgr; k++; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=rpr;
if (!vrsh) { n=dmsrps; mu=NapMasRasLePrSha(n, legr, prgr, raspr, srra, rpr, mu); p0=raspr[0]; }
else if (vrsh==1) { n=dmsrpv; mu=NapMasRasLePrSha2(n, legr, prgr, raspr, srra, rpr, mu); p0=raspr[0]; }
else if (vrsh==2) { n=dmrprs3; mu=NapMasRasLePrSha3(n, legr, prgr, raspr, srra, rpr, mu); p0=25.2*ko; }
else if (vrsh==3) { n=dmrprs4; mu=NapMasRasLePrSha4(n, legr, prgr, raspr, srra, rpr, mu); p0=26.5*ko; }
else if (vrsh==4) { n=dmrprs5; mu=NapMasRasLePrSha5(n, legr, prgr, raspr, srra, rpr, mu); p0=11.5*ko; }
else if (vrsh==5) { n=dmrprs6; mu=NapMasRasLePrSha6(n, legr, prgr, raspr, srra, rpr, mu); p0=raspr[0]; }
if (vrsh==1) ko=9.962*r; else if (vrsh==2) ko=14.0*r; else if (vrsh==3) ko=18.0*r; 
else if (vrsh==4) ko=15.0*r; else if (vrsh==5) ko=13.0*r; cout << "Sred razm por = " << ko << "\t";
k=0; legr=mu[k]; k++; prgr=mu[k]; k++; raspr=mu[k]; k++; srra=mu[k]; k++; rpr=mu[k]; 
if (rpr) delete[]rpr; if (srra) delete[]srra;
mu=vybRasPorRazSha(porsha, p0, n, k, prgr, legr, raspr, mu);
if (raspr) delete[]raspr; if (prgr) delete[]prgr; if (legr) delete[]legr;
return mu;
}
double *NapMasRaspSha3(double *raspr, int n)
{ int k=0; double ko=1e-2;
raspr[k]=0.0;     k++;	raspr[k]=10.6294; k++; raspr[k]=14.0694; k++; raspr[k]=17.8574; k++; raspr[k]=23.8485; k++;
raspr[k]=33.3956; k++;	raspr[k]=39.3094; k++; raspr[k]=44.8367; k++; raspr[k]=50.0161; k++; raspr[k]=53.8813; k++;
raspr[k]=59.6405; k++;	raspr[k]=61.4827; k++; raspr[k]=63.8741; k++; raspr[k]=68.1047; k++; raspr[k]=71.3586; k++;
raspr[k]=75.9268; k++;	raspr[k]=79.3882; k++; raspr[k]=83.2566; k++; raspr[k]=86.9122; k++; raspr[k]=89.9145; k++;
raspr[k]=92.5579; k++;	raspr[k]=95.1927; k++; raspr[k]=96.8355; k++; raspr[k]=97.998;  k++; raspr[k]=98.8324; k++;
raspr[k]=99.3257; k++;	raspr[k]=99.6412; k++; raspr[k]=99.9201; k++; raspr[k]=99.9442; k++; raspr[k]=99.9575; k++;
raspr[k]=1e2; k++;
for (k=0; k<n; k++) raspr[k]=ko*raspr[k];
int j=(n-(n%2))/2; double t;
for (k=0; k<j; k++) { t=raspr[k]; raspr[k]=raspr[n-1-k]; raspr[n-1-k]=t; }
return raspr;
}
double *NapMasPrGrSha3(double *rapo, int n)
{ int k=0, j=(n-(n%2))/2; double ko=1e-6, t=0.0;
rapo[k]=340.838;  k++;  rapo[k]=200.63933; k++; rapo[k]=45.34257; k++; rapo[k]=26.45502; k++; rapo[k]=20.39776; k++;
rapo[k]=16.04155; k++;  rapo[k]=12.64329;  k++; rapo[k]=10.19442; k++; rapo[k]=8.15041;  k++; rapo[k]=6.63919;  k++;
rapo[k]=5.27988;  k++;  rapo[k]=4.16287;   k++; rapo[k]=3.46842;  k++; rapo[k]=2.80983;  k++; rapo[k]=2.2606;   k++;
rapo[k]=1.81613;  k++;  rapo[k]=1.44981;   k++; rapo[k]=1.17497;  k++; rapo[k]=0.93713;  k++; rapo[k]=0.75399;  k++;
rapo[k]=0.61365;  k++;  rapo[k]=0.49093;   k++; rapo[k]=0.39032;  k++; rapo[k]=0.31608;  k++; rapo[k]=0.2547;   k++;
rapo[k]=0.20469;  k++;  rapo[k]=0.16692;   k++; rapo[k]=0.1359;   k++; rapo[k]=0.10807;  k++; rapo[k]=0.0862;   k++;
rapo[k]=0.06979;  k++;
for (k=0; k<n; k++) rapo[k]=ko*rapo[k]; 
for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; } 
return rapo; }
double **NapMasRasLePrSha3(int n, double *legr, double *prgr, double *raspr, double *srra, double *rpr, double **mu)
{ int k=5; double p0=1e-2, ko=1e-6; 
raspr=NapMasRaspSha3(raspr,n); 
for (k=1; k<n; k++) rpr[k]=fabs(raspr[k-1]-raspr[k]);
prgr=NapMasPrGrSha3(prgr,n); 
k=0; legr[k]=ko*p0; for (k=1; k<n; k++) legr[k]=prgr[k-1];
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2e0; /*средний размер каждого из диапазонов*/
k=0; mu[k]=legr; k++; mu[k]=prgr; k++; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=rpr;
return mu; }

//ШВ-1 № 1 1-1
double *NapMasRaspSha4(double *raspr, int n)
{ int k=0; double ko=1e-2;
raspr[k]=0.0;	  k++; raspr[k]=11.5833; k++; raspr[k]=16.7564; k++; raspr[k]=23.6914; k++; raspr[k]=33.213;  k++;
raspr[k]=42.772;  k++; raspr[k]=49.4821; k++; raspr[k]=55.7423; k++; raspr[k]=60.9905; k++; raspr[k]=65.1889; k++;
raspr[k]=70.437;  k++; raspr[k]=71.3838; k++; raspr[k]=72.3753; k++; raspr[k]=73.7401; k++; raspr[k]=76.0907; k++;
raspr[k]=79.529;  k++; raspr[k]=82.1748; k++; raspr[k]=85.175;  k++; raspr[k]=88.2587; k++; raspr[k]=90.8563; k++;
raspr[k]=93.1275; k++; raspr[k]=95.443;  k++; raspr[k]=96.8905; k++; raspr[k]=97.874;  k++; raspr[k]=98.6652; k++;
raspr[k]=99.0294; k++; raspr[k]=99.3704; k++; raspr[k]=99.6077; k++; raspr[k]=99.7959; k++; raspr[k]=99.9616; k++;
raspr[k]=99.9757; k++; raspr[k]=1e2; k++;
for (k=0; k<n; k++) raspr[k]=ko*raspr[k];
int j=(n-(n%2))/2; double t;
for (k=0; k<j; k++) { t=raspr[k]; raspr[k]=raspr[n-1-k]; raspr[n-1-k]=t; }
return raspr;
}
double *NapMasPrGrSha4(double *rapo, int n)
{ int k=0; double ko=1e-6;
rapo[k]=341.24803; k++; rapo[k]=200.87079; k++; rapo[k]=45.36211; k++; rapo[k]=26.43548; k++; rapo[k]=20.38867; k++;
rapo[k]=16.04315;  k++; rapo[k]=12.64066;  k++; rapo[k]=10.19261; k++; rapo[k]=8.14948;  k++; rapo[k]=6.63862;  k++;
rapo[k]=5.27717;   k++; rapo[k]=4.144;	   k++; rapo[k]=3.41414;  k++; rapo[k]=2.76709;  k++; rapo[k]=2.27202;  k++;
rapo[k]=1.82495;   k++; rapo[k]=1.44478;   k++; rapo[k]=1.18063;  k++; rapo[k]=0.93948;  k++; rapo[k]=0.75206;  k++;
rapo[k]=0.61279;   k++; rapo[k]=0.49187;   k++; rapo[k]=0.39041;  k++; rapo[k]=0.31641;  k++; rapo[k]=0.25516;  k++;
rapo[k]=0.20479;   k++; rapo[k]=0.16705;   k++; rapo[k]=0.13588;  k++; rapo[k]=0.10807;  k++; rapo[k]=0.08619;  k++;
rapo[k]=0.06981;   k++; rapo[k]=0.05648;  k++;
for (k=0; k<n; k++) rapo[k]=ko*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; } 
return rapo; }
double **NapMasRasLePrSha4(int n, double *legr, double *prgr, double *raspr, double *srra, double *rpr, double **mu)
{ int k=5; double p0=1e-6, ko=1e-2;
raspr=NapMasRaspSha4(raspr,n); for (k=1; k<n; k++) rpr[k]=fabs(raspr[k-1]-raspr[k]);
prgr=NapMasPrGrSha4(prgr,n); 
k=0; legr[k]=p0*ko; for (k=1; k<n; k++) legr[k]=prgr[k-1];
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2e0; /*средний размер каждого из диапазонов*/ 
k=0; mu[k]=legr; k++; mu[k]=prgr; k++; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=rpr;
return mu;
}

//ШПД
double *NapMasRaspSha5(double *raspr, int n)
{ int k=0;
raspr[k]=0.0;     k++; raspr[k]=10.1713; k++; raspr[k]=18.4934; k++; raspr[k]=26.538;  k++; raspr[k]=36.3394; k++;
raspr[k]=48.4526; k++; raspr[k]=56.4047; k++; raspr[k]=65.0966; k++; raspr[k]=71.6618; k++; raspr[k]=76.5625; k++;
raspr[k]=83.4975; k++; raspr[k]=84.9653; k++; raspr[k]=87.1196; k++; raspr[k]=90.9463; k++; raspr[k]=93.8334; k++;
raspr[k]=96.7929; k++; raspr[k]=97.6675; k++; raspr[k]=98.5212; k++; raspr[k]=99.0279; k++; raspr[k]=99.396;  k++;
raspr[k]=99.7339; k++; raspr[k]=99.8082; k++; raspr[k]=99.9293; k++; raspr[k]=99.9293; k++; raspr[k]=1e2;     k++;
for (k=0; k<n; k++) raspr[k]=(1e-2)*raspr[k];
int j=(n-(n%2))/2; double t;
for (k=0; k<j; k++) { t=raspr[k]; raspr[k]=raspr[n-1-k]; raspr[n-1-k]=t; }
return raspr;
}
double *NapMasPrGrSha5(double *rapo, int n)
{ int k=0;
rapo[k]=343.22645; k++; rapo[k]=201.84496; k++; rapo[k]=45.34218; k++; rapo[k]=26.43562; k++;  rapo[k]=20.39193; k++;
rapo[k]=16.04629;  k++; rapo[k]=12.64759;  k++; rapo[k]=10.19305; k++; rapo[k]=8.14791;  k++;  rapo[k]=6.6392;   k++;
rapo[k]=5.27874;   k++; rapo[k]=4.12243;   k++; rapo[k]=3.40431;  k++; rapo[k]=2.77495;  k++;  rapo[k]=2.25745;  k++;
rapo[k]=1.82966;   k++; rapo[k]=1.45675;   k++; rapo[k]=1.1781;   k++; rapo[k]=0.94135;  k++;  rapo[k]=0.75467;  k++;
rapo[k]=0.61214;   k++; rapo[k]=0.49046;   k++;  rapo[k]=0.3903;  k++; rapo[k]=0.31611;  k++;  rapo[k]=0.2549;   k++;
for (k=0; k<n; k++) rapo[k]=(1e-6)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; } 
return rapo; }
double **NapMasRasLePrSha5(int n, double *legr, double *prgr, double *raspr, double *srra, double *rpr, double **mu)
{ int k=0; double p0=1e-6, ko=1e-2;
raspr=NapMasRaspSha5(raspr,n); 
for (k=1; k<n; k++) rpr[k]=fabs(raspr[k-1]-raspr[k]);
prgr=NapMasPrGrSha5(prgr,n); 
k=0; legr[0]=ko*p0; for (k=1; k<n; k++) legr[k]=prgr[k-1];
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2e0; /*средний размер каждого из диапазонов*/ 
k=0; mu[k]=legr; k++; mu[k]=prgr; k++; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=rpr;
return mu;
}

//ШКУ-32-3-1
double *NapMasRaspSha6(double *raspr, int n)
{ int k=0;
raspr[k]=0.0;     k++; raspr[k]=12.1030; k++; raspr[k]=14.4598; k++; raspr[k]=17.4058; k++; raspr[k]=22.6349; k++;
raspr[k]=31.3009; k++; raspr[k]=38.0275; k++; raspr[k]=44.2632; k++; raspr[k]=50.1060; k++; raspr[k]=54.9669; k++;
raspr[k]=62.0126; k++; raspr[k]=64.3756; k++; raspr[k]=66.5437; k++; raspr[k]=69.6054; k++; raspr[k]=73.0323; k++;
raspr[k]=77.4105; k++; raspr[k]=80.8063; k++; raspr[k]=84.1626; k++; raspr[k]=87.4303; k++; raspr[k]=89.9592; k++;
raspr[k]=91.9952; k++; raspr[k]=93.9316; k++; raspr[k]=95.2501; k++; raspr[k]=96.3643; k++; raspr[k]=97.3739; k++;
raspr[k]=98.0557; k++; raspr[k]=98.6241; k++; raspr[k]=99.0612; k++; raspr[k]=99.3657; k++; raspr[k]=99.6501; k++;
raspr[k]=99.8228; k++; raspr[k]=99.9420; k++; raspr[k]=1e2; 
for (k=0; k<n; k++) raspr[k]=(1e-2)*raspr[k];
int j=(n-(n%2))/2; double t;
for (k=0; k<j; k++) { t=raspr[k]; raspr[k]=raspr[n-1-k]; raspr[n-1-k]=t; }
return raspr;
}
double *NapMasPrGrSha6(double *rapo, int n)
{ int k=0;
rapo[k]=341.18288; k++; rapo[k]=200.85458; k++; rapo[k]=45.37898; k++; rapo[k]=26.44827; k++; rapo[k]=20.39319; k++;
rapo[k]=16.03650;  k++; rapo[k]=12.62806;  k++; rapo[k]=10.17641; k++; rapo[k]=8.14755;  k++; rapo[k]=6.63955;  k++;
rapo[k]=5.27949;   k++; rapo[k]=4.12115;   k++; rapo[k]=3.42943;  k++; rapo[k]=2.82385;  k++; rapo[k]=2.27122;  k++;
rapo[k]=1.81876;   k++; rapo[k]=1.44898;   k++; rapo[k]=1.17408;  k++; rapo[k]=0.93791;  k++; rapo[k]=0.75018;  k++;
rapo[k]=0.61214;   k++; rapo[k]=0.49068;   k++; rapo[k]=0.39039;  k++; rapo[k]=0.31661;  k++; rapo[k]=0.25516;  k++;
rapo[k]=0.20487;   k++; rapo[k]=0.16710;   k++; rapo[k]=0.13595;  k++; rapo[k]=0.10807;  k++; rapo[k]=0.08620;  k++;
rapo[k]=0.06981;   k++; rapo[k]=0.05647;   k++; rapo[k]=0.04533; k++;
for (k=0; k<n; k++) rapo[k]=(1e-6)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; } 
return rapo; }
double **NapMasRasLePrSha6(int n, double *legr, double *prgr, double *raspr, double *srra, double *rpr, double **mu)
{ int k=0; double p0=0.0, ko=0.0;
raspr=NapMasRaspSha6(raspr,n); 
for (k=1; k<n; k++) rpr[k]=fabs(raspr[k-1]-raspr[k]);
prgr=NapMasPrGrSha6(prgr,n); 
p0=1e-6; ko=1e-2; legr[0]=p0*ko; for (k=1; k<n; k++) legr[k]=prgr[k-1];
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2e0; /*средний размер каждого из диапазонов*/ 
k=0; mu[k]=legr; k++; mu[k]=prgr; k++; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=rpr;
return mu; }

double *NapMasRaspSha2(double *raspr, int n)
{ int k=0;
raspr[k]=0.0; k++;   raspr[k]=0.0; k++;   raspr[k]=0.0; k++;  raspr[k]=0.23; k++;
raspr[k]=0.29; k++;  raspr[k]=0.4; k++;   raspr[k]=0.46; k++; raspr[k]=0.75; k++;
raspr[k]=2.48; k++;  raspr[k]=2.48; k++;  raspr[k]=7.86; k++; raspr[k]=10.8; k++;
raspr[k]=23.92; k++; raspr[k]=41.88; k++; raspr[k]=1.79; k++; raspr[k]=1.04; k++;
raspr[k]=3.64; k++;  raspr[k]=1.96; k++;
for (k=0; k<n; k++) raspr[k]=(1e-2)*raspr[k]; 
int j=(n-(n%2))/2; double t; 
for (k=0; k<j; k++) { t=raspr[k]; raspr[k]=raspr[n-1-k]; raspr[n-1-k]=t; }
return raspr; }
double *NapMasPrGrSha2(double *rapo, int n)
{ int k=0; 
rapo[k]=13e1; k++; rapo[k]=12e1; k++; rapo[k]=11e1; k++; rapo[k]=1e2; k++; rapo[k]=9e1; k++;
rapo[k]=8e1;  k++; rapo[k]=7e1;  k++; rapo[k]=6e1;  k++; rapo[k]=5e1; k++; rapo[k]=4e1; k++;
rapo[k]=3e1;  k++; rapo[k]=2e1;  k++; rapo[k]=15.0; k++; rapo[k]=1e1; k++; rapo[k]=5.0; k++;
rapo[k]=3.0;  k++; rapo[k]=1.0;  k++; rapo[k]=1e-1;
for (k=0; k<n; k++) rapo[k]=(1e-6)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; } 
return rapo; }
double **NapMasRasLePrSha2(int n, double *legr, double *prgr, double *raspr, double *srra, double *rpr, double **mu)
{ int k=5; double p0=0.0, ko=0.0;
raspr=NapMasRaspSha2(raspr, n); 
for (k=0; k<n; k++) rpr[k]=raspr[k];
prgr=NapMasPrGrSha2(prgr,n); 
p0=1e-2; ko=1e-6; k=0; legr[k]=p0*ko; for (k=1; k<n; k++) legr[k]=prgr[k-1];
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2e0; /*средний размер каждого из диапазонов*/ 
k=0; mu[k]=legr; k++; mu[k]=prgr; k++; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=rpr;
return mu; }

double *NapMasRaspSha(double *raspr)
{ int k=0;
raspr[k]=21.8;  k++; raspr[k]=21.75;  k++; raspr[k]=21.625; k++; raspr[k]=21.25; k++; 
raspr[k]=20.75; k++; raspr[k]=20.25;  k++; raspr[k]=19.75;  k++; raspr[k]=19.0;  k++; 
raspr[k]=18.25; k++; raspr[k]=17.5;   k++; raspr[k]=12.0;   k++; raspr[k]=9.875; k++; 
raspr[k]=8.750; k++; raspr[k]=8.25;   k++; raspr[k]=7.25;   k++; raspr[k]=6.250; k++; 
raspr[k]=5.0;   k++; raspr[k]=4.25;   k++; raspr[k]=3.25;   k++; raspr[k]=0.875; k++; 
raspr[k]=0.325; k++; raspr[k]=0.0; k++; raspr[k]=0.0; 
return raspr; }
double *NapMasPrGrSha(double *prgr, int n)
{ int k=0; double ko=1e-6;
prgr[k]=5e-2; k++; prgr[k]=1e-1; k++; prgr[k]=4e-1; k++; prgr[k]=5e-1; k++; 
prgr[k]=6e-1; k++; prgr[k]=7e-1; k++; prgr[k]=8e-1; k++; prgr[k]=9e-1; k++; 
prgr[k]=1.00; k++; prgr[k]=2.00; k++; prgr[k]=3.00; k++; prgr[k]=4.00; k++; 
prgr[k]=5.00; k++; prgr[k]=6.00; k++; prgr[k]=7.00; k++; prgr[k]=8.00; k++; 
prgr[k]=9.00; k++; prgr[k]=1e1;  k++; prgr[k]=2e1;  k++; prgr[k]=3e1;  k++;
prgr[k]=4e1;  k++; prgr[k]=5e1;  k++; prgr[k]=13e1;
for (k=0; k<n; k++) prgr[k]=prgr[k]*ko; /*в метрах*/ 
return prgr; }
double **NapMasRasLePrSha(int n, double *legr, double *prgr, double *raspr, double *srra, double *rpr, double **mu)
{ int k=0; double p0;
raspr=NapMasRaspSha(raspr); k=0; p0=raspr[k]; p0=1e2/p0;
for (k=1; k<n; k++) rpr[k]=fabs(raspr[k-1]-raspr[k])*p0;
prgr=NapMasPrGrSha(prgr, n); 
k=0; legr[k]=0.0; for (k=1; k<n; k++) legr[k]=prgr[k-1];
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2.0; /*средний размер каждого из диапазонов*/ 
k=0; mu[k]=legr; k++; mu[k]=prgr; k++; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=rpr;
return mu; }

double **vybRasPorRazSha(double poris, double raspr0, int n, int vlpr, double *prgr, double *legr, double *rpr, double **mu)
{
double ko=1e-2, rpn=1e0, dp=1e0, ht=dp, p0=0.0; if (raspr0<ht) ko=ht; p0=raspr0*ko;
if ((poris<stpoSha20) && (poris>=stpoSha16)) mu=poisrasprpor1(poris, p0, n, rpn, dp, vlpr, prgr, legr, rpr, mu);
else if ((poris>stpoSha30) && (poris<stpoSha33)) mu=poisrasprpor2(poris, p0, n, rpn, dp, vlpr, prgr, legr, rpr, mu);
else if ((poris>=p0) && (poris<stpoSha30)) mu=poisrasprpor2(poris, p0, n, rpn, dp, vlpr, prgr, legr, rpr, mu);
else if ((poris<p0) && (poris>stpoSha20)) mu=poisrasprpor1(poris, p0, n, rpn, dp, vlpr, prgr, legr, rpr, mu);
else if ((poris<stpoSha16) && (poris>=stpoSha10)) mu=poisrasprpor1(poris, p0, n, rpn, dp, vlpr, prgr, legr, rpr, mu);
return mu;
}

//Вермикулит
double *PravGranPoromVer(double *rapo, int n)
{ int k=0; double ko=1e-6;
rapo[k]=13e1; k++; rapo[k]=12e1; k++; rapo[k]=11e1; k++; rapo[k]=1e2; k++; rapo[k]=9e1; k++;
rapo[k]=8e1;  k++; rapo[k]=7e1;  k++; rapo[k]=6e1;  k++; rapo[k]=5e1; k++; rapo[k]=4e1; k++;
rapo[k]=3e1;  k++; rapo[k]=2e1;  k++; rapo[k]=15.0; k++; rapo[k]=1e1; k++; rapo[k]=5.0; k++;
rapo[k]=3.0;  k++; rapo[k]=1.0;  k++; rapo[k]=1e-1; k++;
for (k=0; k<n; k++) rapo[k]=ko*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; } 
return rapo; }
double *LevGranPoromVer(double *pravgr, double *levgr, int n)
{ int k; double ko=1e-6, mr=1e-2; levgr[0]=ko*mr; for (k=1; k<n; k++) levgr[k]=pravgr[k-1]; return levgr; }
double *rasPorpoRazmVerIs84(double *rapo, int n)
{ int k=0; 
rapo[k]=0.00; k++; rapo[k]=0.00; k++; rapo[k]=2.14; k++; rapo[k]=1.000; k++; rapo[k]=1.06; k++;
rapo[k]=1.33; k++; rapo[k]=1.30; k++; rapo[k]=1.57; k++; rapo[k]=1.900; k++; rapo[k]=2.35; k++;
rapo[k]=3.83; k++; rapo[k]=3.77; k++; rapo[k]=8.94; k++; rapo[k]=27.38; k++; rapo[k]=7.70; k++;
rapo[k]=29.0; k++; rapo[k]=6.73; k++; rapo[k]=0.0;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *rasPorpoRazmVerIn84(double *rapo, int n)
{ int k=0; 
rapo[k]=0.0;   k++; rapo[k]=0.0;  k++; rapo[k]=0.34;  k++; rapo[k]=0.83;  k++; rapo[k]=0.86; k++;
rapo[k]=1.24;  k++; rapo[k]=1.31; k++; rapo[k]=1.66;  k++; rapo[k]=2.28;  k++; rapo[k]=3.55; k++;
rapo[k]=8.14;  k++; rapo[k]=9.38; k++; rapo[k]=16.11; k++; rapo[k]=25.01; k++; rapo[k]=5.62; k++;
rapo[k]=15.14; k++; rapo[k]=7.0;  k++; rapo[k]=1.52;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *rasPorpoRazmVerI207(double *rapo, int n)
{ int k=0;
rapo[k]=0.00; k++; rapo[k]=0.00; k++; rapo[k]=2e-2; k++; rapo[k]=1.860; k++;
rapo[k]=1.45; k++; rapo[k]=1.54; k++; rapo[k]=1.72; k++; rapo[k]=1.890; k++;
rapo[k]=2.28; k++; rapo[k]=2.57; k++; rapo[k]=3.79; k++; rapo[k]=2.990; k++;
rapo[k]=5.77; k++; rapo[k]=19.8; k++; rapo[k]=9.00; k++; rapo[k]=45.31; k++; 
rapo[k]=0.21; k++; rapo[k]=0.00; k++;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *rasPorpoRazmVerI16035(double *rapo, int n)
{ int k=0;
rapo[k]=0.00; k++; rapo[k]=0.00; k++; rapo[k]=2e-2; k++; rapo[k]=1.860; k++;
rapo[k]=1.45; k++; rapo[k]=1.54; k++; rapo[k]=1.72; k++; rapo[k]=1.890; k++;
rapo[k]=2.28; k++; rapo[k]=2.57; k++; rapo[k]=3.79; k++; rapo[k]=2.990; k++;
rapo[k]=5.77; k++; rapo[k]=19.8; k++; rapo[k]=9.00; k++; rapo[k]=45.31; k++; 
rapo[k]=0.21; k++; rapo[k]=0.00; k++;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *rasPorpoRazmVerO84(double *rapo, int n)
{ int k=0; 
rapo[k]=0.0;   k++; rapo[k]=0.0;  k++; rapo[k]=0.33;  k++; rapo[k]=0.93;  k++; rapo[k]=0.96; k++;
rapo[k]=1.19;  k++; rapo[k]=1.3;  k++; rapo[k]=1.7;   k++; rapo[k]=1.93;  k++; rapo[k]=2.67; k++; 
rapo[k]=4.97;  k++; rapo[k]=5.48; k++; rapo[k]=13.52; k++; rapo[k]=28.46; k++; rapo[k]=9.41; k++; 
rapo[k]=17.93; k++; rapo[k]=7.26; k++; rapo[k]=1.96;  k++;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo; }
double *rasPorpoRazmVerO16035(double *rapo, int n)
{ 
int k=0;
rapo[k]=0.00; k++; rapo[k]=0.00; k++; rapo[k]=2e-2; k++; rapo[k]=1.860; k++;
rapo[k]=1.45; k++; rapo[k]=1.54; k++; rapo[k]=1.72; k++; rapo[k]=1.890; k++;
rapo[k]=2.28; k++; rapo[k]=2.57; k++; rapo[k]=3.79; k++; rapo[k]=2.990; k++;
rapo[k]=5.77; k++; rapo[k]=19.8; k++; rapo[k]=9.00; k++; rapo[k]=45.31; k++; 
rapo[k]=0.21; k++; rapo[k]=0.00; k++;
for (k=0; k<n; k++) rapo[k]=(1e-2)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; }
return rapo;
}
double *NapMasRaspVer(int vfv, int n, int vysove, int vpkf, double *raspr)
{ int k=0;
if ((!vysove) || (vysove==2)) { 
if (!vfv) raspr=rasPorpoRazmVerI207(raspr, dmsrpv); //для фракции 2-0,7 мм, исходный
else if (vfv==1) { if (!vpkf) raspr=rasPorpoRazmVerIs84(raspr, dmsrpv); //для фракции 8-4 мм, исходный
else if (vpkf==1) raspr=rasPorpoRazmVerIn84(raspr, dmsrpv); } //для фракции 8-4 мм, исходный, другая пористость
else if (vfv==2) raspr=rasPorpoRazmVerI16035(raspr, dmsrpv); } //для фракции 1,6-0,35 мм, исходный
if (vysove==1) {
if ((vfv==2) || (!vfv)) raspr=rasPorpoRazmVerO16035(raspr, dmsrpv); //для фракции 1,6-0,35 мм и 2-0,7 мм, после обжига
else if (vfv==1) raspr=rasPorpoRazmVerO84(raspr, dmsrpv); //для фракции 8-4 мм, после обжига
else { cout << "Net takoy fraktsii!"; k=getchar(); exit(1); } } //только для фракции 8-4 мм
return raspr; }
double **NapMasRasLePrVer(int vfv, int vysove, int vpkf, int nom, double **mu)
{ int k, n=dmsrpv, q;
double *legr=NULL, *prgr=NULL, *raspr=NULL;
double *srra=NULL, ht=1e0, t, r;
legr=new double[n]; prgr=new double[n]; raspr=new double[n]; srra=new double[n]; 
if ((!legr) || (!prgr) || (!raspr) || (!srra) || (!mu)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { raspr[k]=0.0; prgr[k]=0.0; legr[k]=0.0; srra[k]=0.0; }
raspr=NapMasRaspVer(vfv, n, vysove, vpkf, raspr);
prgr=PravGranPoromVer(prgr, n);
legr=LevGranPoromVer(prgr, legr, n);
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2e0;
q=1; t=0.0; double m=0.0, *mm=new double[q], *ms=new double[q];
m=0.0; for (k=0; k<n; k++) m=m+ht; k=0; mm[k]=m;
t=0.0; r=0.0; for (k=0; k<n; k++) { t=t+srra[k]*raspr[k]; r=r+raspr[k]; } 
k=0; if (fabs(r)>0.0) t=t/r; else t=0.0; ms[k]=t;
k=0; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=prgr; k++; mu[k]=legr; k++; mu[k]=ms; k++; mu[k]=mm;
return mu; } //средний размер каждого из диапазонов
double *vybRasPorRazmVer(double *rpr, int vysove, int vfv, int vpkf, int n)
{
if ((!vysove) || (vysove==1)) {
if (!vfv) rpr=rasPorpoRazmVerI207(rpr,n);
else if (vfv==1) { if (!vpkf) rpr=rasPorpoRazmVerIs84(rpr,n); 
else if (vpkf==1) rpr=rasPorpoRazmVerIn84(rpr,n); }
else if (vfv==2) rpr=rasPorpoRazmVerI16035(rpr,n); }
if (vysove==2) { if ((!vfv) || (vfv==2)) rpr=rasPorpoRazmVerO16035(rpr,n);
else if (vfv==1) rpr=rasPorpoRazmVerO84(rpr,n); }
return rpr;
}
double **rasPorpoRazVer(double poris, int vfv, int vyvyzn, int vysove, int isrp, int vpkf) //0 - старые, 1 - новые значения
{ int k=0, n=0, j=0, qg=0, nm=cvym, q=1;
double *legr=NULL, *prgr=NULL, *raspr=NULL, *srra=NULL, ht=1e0, s=0.0, e=1e-1;
double rpn=1e0, dp=1e0, ras0=0.0, marp=0.0, srp=0.0, **mu=NULL, *po=NULL, r=0.0;
double *ms, *mm, koef=1e-6;
mu=new double*[nm]; if (!mu) { cout << "No memory" << endl; k=getchar(); exit(1); } 
mu=NapMasRasLePrVer(vfv, vysove, vpkf, k, mu);
k=0; raspr=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++;
po=mu[k]; j=0; srp=po[j]; k++; if (po) delete[]po; 
po=mu[k]; marp=po[j]; if (po) delete[]po; 
k=0; s=e; while (s<marp) { s=s+ht; k++; } n=k; //cout << "srp = " << srp << "\tmarp = " << marp << "\tn = " << n << endl; //for (k=0; k<n; k++) cout << "rpr ( " << k << " ) = " << raspr[k] << "\tsrra ( " << k << " ) = " << srra[k] << endl;
k=0; delete[]mu; mu=new double*[nm]; if (!mu) { cout << "No memory" << endl; k=getchar(); exit(1); } 
mu=rasrasporporaz0(prgr, legr, raspr, rpn, k, n, dp, mu);
s=0.0; double l=s, p=1e-6, h=p, ko=ht/h; 
if (legr) delete[]legr; if (prgr) delete[]prgr; 
if (raspr) delete[]raspr; if (srra) delete[]srra; 
k=0; raspr=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++;
ms=mu[k]; j=0; srp=ms[j]; k++; mm=mu[k]; marp=mm[j]; 
k=0; s=e; while (s<marp) { s=s+ht; k++; } qg=k; s=srp; l=marp; 
if (isrp==1) { //расчет характеристик пористой структуры при условии шарообразных пор
double sz=0.0, sch=sz, dv=sz, vpre=sz, vtek=sz, rtek=sz, rpre=sz; l=0.0; p=1e-6;
for (k=0; k<qg; k++) { rtek=(p+l)/2e0; srra[k]=rtek; 
vtek=(pi/6e0)*pow(rtek,3e0); vpre=(pi/6e0)*pow(rpre,3e0); dv=fabs(vtek-vpre);
sz=sz+raspr[k]; sch=sch+raspr[k]*dv; rpre=rtek;
legr[k]=l; p=p+h; l=l+h; } s=sch/sz; } 
else if (isrp>1) { cout << "Nepravilno vybran nomer" << endl; k=getchar(); exit(1); }
for (k=0; k<qg; k++) { srra[k]=srra[k]*koef; legr[k]=legr[k]*koef; prgr[k]=prgr[k]*koef; }
k=0; s=s*koef; ms[k]=s; l=l*koef; mm[k]=l; //cout << "Sr raz por = " << s << "\tMax Raz Por = " << l << endl;
k=0; mu[k]=raspr; k++; mu[k]=srra; k++; mu[k]=legr; k++; mu[k]=prgr; k++; mu[k]=ms; k++; mu[k]=mm;
return mu; }

double **rasrasporporaz1(double *raspr0, double rpn, int vyb, int nn, double *ras11, double *prgr01, double *legr01, double p1, double p0, double de, int vlpr, double **mu)
{ 	int k, q, j, n=nn, qg=0;
	double *legr11=new double[n], *prgr11=new double[n], kot, dob1, dob2, dob3, xp, kop, xt, s, e=1e-1, ht=1e0;
	if ((!legr11) || (!prgr11)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) prgr11[k]=prgr01[k]*p1/p0; legr11[0]=0.0; for (k=1; k<n; k++) legr11[k]=prgr11[k-1]; //for (k=0; k<n; k++) if (k<50) cout << "k = " << k << "\tr = " << ras11[k] << endl;
	struct raspporporazmm *rpr21=NULL, *rprp21=NULL, *roo=new struct raspporporazmm;
	rpr21=roo; q=0; xp=0.0; xt=0.0;
	for (k=0; k<n; k++) {
    xt=prgr01[k];
    for (j=0; j<n; j++) {
        if ((legr11[j]<xt) && (prgr11[j]>xt)) {
        kot=ras11[j]/fabs(prgr11[j]-legr11[j]);
        dob1=kot*(xt-legr11[j]);
		rpr21->rprm=dob1;
		        if (xp<legr11[j]) {
					if (j>0) {
                kop=ras11[j-1]/fabs(prgr11[j-1]-legr11[j-1]);
                dob2=(legr11[j]-xp)*kop; 
				rpr21->rprm=rpr21->rprm+dob2; } } //if ((xp>legr11[j]) && (xp<prgr11[j])) { dob3=kot*(xt-xp); rpr21->rprm=dob3; } 
		} } rprp21=rpr21;
		rpr21=new raspporporazmm;
		rpr21->rprm=0.0;
		rprp21->slel=rpr21;
		rpr21->slel=NULL;
		q++; 
	xp=xt; }
double *ras21=new double[q], *prgrm=new double[q], *legrm=new double[q];
k=1; double prgr001=rpn, t, *srra21=new double[q], *tm=new double[k], *ts=new double[k], r;
if ((!ras21) || (!prgrm) || (!legrm) || (!srra21) || (!tm) || (!ts)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr21=roo; rprp21=roo; k=0; legrm[0]=0.0; 
while ((rpr21) && (k<q)) { 
	prgrm[k]=prgr001;
	if (k) legrm[k]=prgrm[k-1];
	ras21[k]=rprp21->rprm; 
	rpr21=rprp21->slel; 
	delete rprp21;
	rprp21=rpr21;
	prgr001=prgr001+de; k++; }
q=k; qg=q; s=0.0; for (k=0; k<q; k++) { s=s+ras21[k]; srra21[k]=(legrm[k]+prgrm[k])/2e0; }
if (s<(ht+e)) t=ht; else t=1e2; 
for (k=0; k<qg; k++) ras21[k]=ras21[k]*t/s;
t=0.0; for (k=0; k<q; k++) t=t+ht; //for (k=0; k<q; k++) if (k<nk) printf("r = %0.4lf\tpr = %0.0lf\n",ras21[k],prgrm[k]);
if (prgr11) delete []prgr11; if (legr11) delete []legr11; 
s=0.0; r=0.0; for (k=0; k<qg; k++) { s=s+srra21[k]*ras21[k]; r=r+ras21[k]; } 
k=0; if (fabs(r)>0.0) s=s/r; else s=0.0;
k=0; tm[k]=t; ts[k]=s;
k=0; mu[k]=ras21; k++; mu[k]=srra21; k++; mu[k]=prgrm; k++; 
mu[k]=legrm; k++; mu[k]=ts; k++; mu[k]=tm;
return mu; }
double **poisrasprpor1(double p1, double p0, int nn, double rpn, double dp, int vlpr, double *prgr, double *legr, double *rpr, double **mu)
{ int n=nn, k=0, nk=0, j=0, qg=0;
mu=rasrasporporaz0(prgr, legr, rpr, rpn, k, n, dp, mu);
double ko=1e-6, kk=1e-2, s=0.0, ht=1e0;
double *ras01=NULL, *prgr01=NULL, *legr01=NULL, *qgp=NULL, qgf=0.0, *srra=NULL, *ms=NULL, srp=0.0;
k=0; j=k; ras01=mu[k]; k++; srra=mu[k]; k++; prgr01=mu[k]; k++; legr01=mu[k]; k++;
ms=mu[k]; srp=ms[j]; k++; if (ms) delete[]ms; 
qgp=mu[k]; qgf=qgp[j]; if (qgp) delete[]qgp; if (mu) delete[]mu;
s=-kk; k=0; while (s<qgf) { k++; s=s+ht; } qg=k;
mu=rasrasporporaz1(ras01, rpn, 0, n, ras01, prgr01, legr01, p1, p0, dp, vlpr, mu);
if (srra) delete[]srra; if (ras01) delete[]ras01; 
if (prgr01) delete[]prgr01; if (legr01) delete[]legr01;
return mu; }
double **rasppor2(double p2, double p0, double *legr01, double *prgr01, double *ras21, int nn, double rpn, double de, int vlpr, double **mu)
{ int k, n=nn, qg; 
double *legr21=new double[n], *prgr21=new double[n], e=1e-1, ht=1e0;
if ((!prgr21) || (!legr21)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
struct raspporporazmm *rpr22, *rprp22, *roo=new struct raspporporazmm;
for (k=0; k<n; k++) { legr21[k]=legr01[k]*p2/p0; prgr21[k]=prgr01[k]*p2/p0; }
int q=0, j; double xp=0.0, xt=0.0, dob=0.0, kot=0.0, kop=0.0, s;
rpr22=roo; 
for (k=0; k<n; k++) {
    xt=prgr01[k];
    for (j=0; j<(n-1); j++) {
    if ((legr21[j]<xt) && (prgr21[j]>xt)) {
        kot=ras21[j]/fabs(prgr21[j]-legr21[j]);
		rpr22->rprm=kot*(xt-legr21[j]);
		if (xp>0.0) { 
			if (xp>legr21[j]) { 
				if (xp<prgr21[j]) {
                dob=kot*(xt-xp); 
				rpr22->rprm=dob; } } 
            if (xp<=legr21[j]) { if (j>0) {
                kop=ras21[j-1]/fabs(prgr21[j-1]-legr21[j-1]);
                rpr22->rprm=rpr22->rprm+(legr21[j]-xp)*kop; } } }
		rprp22=rpr22; 
		rpr22=new raspporporazmm;
		rpr22->rprm=0.0;
		rprp22->slel=rpr22;
		rpr22->slel=NULL;
		q++; } }
    xp=xt; kop=kot; }
double *ras22=new double[q], *prgrm=new double[q], *legrm=new double[q], *srra22=new double[q], prgr002=rpn, t; 
if ((!ras22) || (!prgrm) || (!legrm) || (!srra22)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr22=roo; rprp22=roo; k=0; legrm[0]=0.0; n=q; qg=n;
while ((rpr22) && (k<q)) { 
	prgrm[k]=prgr002;
	if (k>0) legrm[k]=prgrm[k-1];
	ras22[k]=rprp22->rprm; 
	rpr22=rprp22->slel; 
	delete rprp22;
	rprp22=rpr22;
	prgr002=prgr002+de; k++; } 
s=0.0; for (k=0; k<n; k++) { s=s+ras22[k]; srra22[k]=(legrm[k]+prgrm[k])/2e0; } 
if (s<(ht+e)) t=ht; else t=1e2; for (k=0; k<n; k++) ras22[k]=ras22[k]*t/s; 
t=0.0; for (k=0; k<n; k++) t=t+ht; //for (k=0; k<q; k++) printf("r22 = %0.4lf\tpr = %0.1lf\tr21 = %0.4lf\n",ras22[k],srra22[k],ras21[k]);
if (legr21) delete []legr21; if (prgr21) delete []prgr21; 
q=1; double *tm=new double[q], *ts=new double[q], st=0.0, r=0.0;
if ((!tm) || (!ts)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { r=r+ras22[k]; st=st+srra22[k]*ras22[k]; } 
if (fabs(r)>0.0) st=st/r; else st=0.0;
k=0; tm[k]=t; ts[k]=st;
k=0; mu[k]=ras22; k++; mu[k]=srra22; k++; mu[k]=prgrm; k++; 
mu[k]=legrm; k++; mu[k]=ts; k++; mu[k]=tm;
return mu; }
double **poisrasprpor2(double p2, double p0, int nn, double rpn, double h, int vlpr, double *prgr, double *legr, double *rpr, double **mu)
{ int n=nn, k=0, nk=0, j=0, qg=0;
mu=rasrasporporaz0(prgr, legr, rpr, rpn, k, n, h, mu);
double kk=1e-2, s=0.0, ht=1e0;
double *ras01=NULL, *srra01=NULL, *prgr01=NULL, *legr01=NULL, *qgp=NULL, qgf=0.0, *ts=NULL, srp=0.0;
k=0; j=k; ras01=mu[k]; k++; srra01=mu[k]; k++; prgr01=mu[k]; k++; legr01=mu[k]; k++; 
ts=mu[k]; srp=ts[j]; k++; if (ts) delete[]ts;
qgp=mu[k]; qgf=qgp[j]; if (qgp) delete[]qgp; if (mu) delete[]mu;
s=kk; k=0; while (s<qgf) { k++; s=s+ht; } qg=k;
mu=rasppor2(p2, p0, legr01, prgr01, ras01, qg, rpn, h, k, mu);
if (ras01) delete[]ras01; if (prgr01) delete[]prgr01; 
if (legr01) delete[]legr01; if (srra01) delete[]srra01;
return mu; }
double **rasrasporporaz0(double *prgr0, double *legr0, double *raspr0, double rpn, int vyb, int n, double de, double **mu)
{ double s=0.0, e=1e-1, r=0.0, koef=1e6, *prgr00=new double[n];
double *legr00=new double[n], m=0.0, t=0.0, ht=1e0; 
int k=0, q=0, p=0, j=0, jk=0, qg=0;
for (k=0; k<n; k++) prgr00[k]=prgr0[k]*koef;
legr00[0]=0.0; for (k=1; k<n; k++) legr00[k]=prgr00[k-1]; 
m=0.0; for (k=0; k<n; k++) if (m<prgr00[k]) m=prgr00[k]; s=e; k=0; while (s<m) { k++; s=s+ht; } qg=k;
s=0.0; for (k=0; k<n; k++) if (prgr00[k]<(ht+e)) { p=k; s=s+raspr0[k]; } //размер пор до 1 мкм
double *ras01=new double[qg], *prgrm=new double[qg], *legrm=new double[qg];
double *mm=new double[q], *ms=new double[q], *srra01=new double[qg], prgr01=rpn;
if ((!ras01) || (!prgrm) || (!legrm) || (!srra01) || (!ms) || (!mm))
{ cout << "No memory!" << endl; k=getchar(); exit(1); }
k=0; prgr01=de; legrm[k]=0.0; prgrm[k]=prgr01; 
for (k=1; k<qg; k++) { ras01[k]=0.0; prgrm[k]=prgr01; legrm[k]=prgrm[k-1]; prgr01=prgr01+de; }
q=0; ras01[q]=s; q++;
for (k=p; k<n; k++) {
    r=prgr00[k]-legr00[k];
    if ((prgr00[k]>(ht+e)) && (r<(ht+e)) && (r>e)) { ras01[q]=raspr0[k]; q++; }
    if (r>(ht+e)) {
        m=raspr0[k]/r;
		t=r; j=0; while (t>e) { j++; t=t-ht; } jk=j; 
        for (j=0; j<jk; j++) { ras01[q]=m; q++; } } }
if (prgr00) delete []prgr00; if (legr00) delete []legr00; //cout << "q = " << q << "\tqg = " << qg << endl;
t=0.0; for (k=0; k<qg; k++) t=t+ht; k=0; mm[k]=t;
for (k=0; k<qg; k++) srra01[k]=(prgrm[k]+legrm[k])/2e0;
s=0.0; r=s; for (k=0; k<qg; k++) { s=s+ras01[k]*srra01[k]; r=r+ras01[k]; }
if (fabs(r)>0.0) s=s/r; else s=0.0; k=0; ms[k]=s;
k=0; mu[k]=ras01; k++; mu[k]=srra01; k++; mu[k]=prgrm; k++; mu[k]=legrm; k++; 
mu[k]=ms; k++; mu[k]=mm; 
return mu; }

//КВИ
double **rasPorpoRazkvi(int vypl)
{ 
int k=0, f=cvym; //число выходных массивов
double **mu=new double*[f]; if (!mu) { cout << "No memory!" << endl; k=getchar(); exit(1); } //0 - массив распределения пор, 1 - массив средних размеров пор, 2, 3 - левая и правая границы, 4 - средний размер пор, 5 - максимальный
mu=poisrasprporkvi(vypl, mu);
return mu; } //КВ
double *NapMasRaspkvi(int vvk, int n)
{
int k=0; double *raspr;
if (vvk==4) raspr=rasprPorpoRazmAbskvi400();
if (vvk==5) raspr=rasprPorpoRazmAbskvi500();
if (vvk==6) raspr=rasprPorpoRazmAbskvi600();
if (vvk==7) raspr=rasprPorpoRazmAbskvi700();
if (vvk==8) raspr=rasprPorpoRazmAbskvi800();
if (vvk==9) raspr=rasprPorpoRazmAbskvi900();
if (vvk==10) raspr=rasprPorpoRazmAbskvi1000();
return raspr;
}
double **poisrasprporkvi(int n, double **mu)
{ double *raspr0, *prgr0; int n1;
	if (n==4) { raspr0=rasprPorpoRazmAbskvi400(); prgr0=LegrRasPorpoRazmkvi400(); n1=kv4; } //400
	else if (n==5) { raspr0=rasprPorpoRazmAbskvi500(); prgr0=LegrRasPorpoRazmkvi500(); n1=kv5; } //500
	else if (n==6) { raspr0=rasprPorpoRazmAbskvi600(); prgr0=LegrRasPorpoRazmkvi600(); n1=kv6; } //600        
	else if (n==7) { raspr0=rasprPorpoRazmAbskvi700(); prgr0=LegrRasPorpoRazmkvi700(); n1=kv7; } //700
	else if (n==8) { raspr0=rasprPorpoRazmAbskvi800(); prgr0=LegrRasPorpoRazmkvi800(); n1=kv8; } //800
	else if (n==9) { raspr0=rasprPorpoRazmAbskvi900(); prgr0=LegrRasPorpoRazmkvi900(); n1=kv9; } //900
	else if (n==10) { raspr0=rasprPorpoRazmAbskvi1000(); prgr0=LegrRasPorpoRazmkvi1000(); n1=kv10; } //1000
	mu=obrabMaskvi(raspr0, prgr0, n1, mu);
return mu; }
double **obrabMaskvi(double *raspr0, double *prgr0, int n1, double **mu)
{ int v=0, k=0, j=0, n=0, p=0, ri=1, u=0, f=0, q=0, x=0, v0, qg;
double prmi=1e0, h=1e0, nh=0.0, r, pr, s, ko, w, w0, k0=1e2, ht=1e-6, mrp0, mpg0, mpg01, e=1e-1;
nh=0.0; while (nh<prk) { nh=nh+h; k++; } n=k; qg=n; //cout << "qg = " << qg << endl;
double *prgr01=new double[n], *ras01=new double[n], *ras02=new double[n];
double *legr0=new double[n1], *legr01=new double[n], *srra=new double[n], srp=0.0; 
if ((!prgr01) || (!legr0) || (!ras01) || (!ras02)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
k=0; mpg0=prgr0[k]; mrp0=raspr0[k]; for (k=0; k<n1; k++) { if (mpg0<prgr0[k]) mpg0=prgr0[k]; if (mrp0>raspr0[k]) mrp0=raspr0[k]; }
for (k=0; k<n; k++) { prgr01[k]=0.0; ras01[k]=0.0; ras02[k]=0.0; }
k=0; legr0[k]=0.0; for (k=1; k<n1; k++) legr0[k]=prgr0[k-1];
nh=0.0; for (k=0; k<qg; k++) { nh=nh+h; prgr01[k]=nh; } //for (k=0; k<n1; k++) cout << "rasrp ( " << k << " ) = " << raspr0[k] << "\tprgr ( " << k << " ) = " << prgr0[k] << endl;
k=0; mpg01=prgr01[k]; for (k=0; k<qg; k++) if (mpg01<prgr01[k]) mpg01=prgr01[k]; //cout << "max prgr0 = " << mpg0 << "\tmax prgr01 = " << mpg01 << "\tmax rp0 = " << mrp0 << endl;
for (k=0; k<qg; k++) {
    if ((prgr0[k]>=prmi) && (ri>0)) {
        p=k; ri=-1; break; } }
pr=raspr0[p]; w=raspr0[p-1]; pr=pr-w;
s=prgr0[p]; ko=prgr0[p-1]; s=s-ko;
s=pr/s; r=w+s*(prmi-ko);
v=0; w0=raspr0[v]; ras01[v]=w0; v++; 
ras01[v]=r; pr=r; v++; v0=v;
u=0; r=(w0-r)/w0; ras02[u]=r*k0; u++; //cout << "w0 = " << w0 << "\tras01 = " << pr << "\tras02 = " << r*k0 << "\tr = " << r << endl; 
for (k=v0-1; k<qg; k++) { //идем по prgr01, legr01
    f=1; x=-1;
        for (q=p; q<n1; q++) { //идем по prgr0
        if ((legr0[q]<=prgr01[k]) && (f>0) && (prgr0[q]>=prgr01[k])) {
			x=q; f=-1; break; } }
if ((x>0) && (f<0)) {
s=raspr0[x]; r=raspr0[x-1]; s=s-r; 
r=prgr0[x]; w=legr0[x]; r=r-w; ko=s/r;
r=prgr01[k]; w=legr0[x]; r=r-w;
w=raspr0[x-1]; s=w+ko*r; }
else {
r=mpg01-mpg0; w=-mrp0; r=w/r;
w=prgr01[k]-mpg0;
s=mrp0+r*w; }
w=pr; w=(w-s)*k0/w0;
if (v<qg) { ras01[v]=s; v++; } pr=s; 
if (u<qg) { ras02[u]=w; u++; } }  //for (k=0; k<qg; k++) if (k>290) cout << "ras01 ( " << k << " ) = " << ras01[k] << "\tras02 ( " << k << " ) = " << ras02[k] << endl;
s=0.0; for (k=0; k<qg; k++) s=s+ras02[k]; //cout << "v = " << v << "\tu = " << u << "\ts = " << s << endl;
h=1e0; e=1e-1; if (s<(h+e)) w=h; else w=1e2;
for (k=0; k<qg; k++) { r=ras02[k]; r=r*w/s; ras02[k]=r; prgr01[k]=prgr01[k]*ht; }
k=0; legr01[k]=0.0; for (k=1; k<qg; k++) legr01[k]=prgr01[k-1];
for (k=0; k<qg; k++) srra[k]=(legr01[k]+prgr01[k])/2e0;
srp=0.0; w=srp; for (k=0; k<qg; k++) { srp=srp+srra[k]*ras02[k]; w=w+ras02[k]; }
if (fabs(w)>0.0) srp=srp/w; else srp=0.0;
u=cvym; k=1; double *mm=new double[k], *ms=new double[k];
if (legr0) delete[]legr0; if (prgr0) delete[]prgr0; if (raspr0) delete[]raspr0;
s=0.0; for (k=0; k<qg; k++) if (s<prgr01[k]) s=prgr01[k];
k=0; mm[k]=s; ms[k]=srp;
k=0; mu[k]=ras02; k++; mu[k]=srra; k++; 
mu[k]=prgr01; k++; mu[k]=legr01; k++; mu[k]=ms; k++; mu[k]=mm;
if (ras01) delete[]ras01;
return mu; }
double *rasprPorpoRazmAbskvi400()
{ int k=0; double *rpr4=new double[kv4];
if (!rpr4) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr4[k]=12.936; k++; rpr4[k]=12.851; k++; rpr4[k]=12.667; k++; rpr4[k]=12.500; k++; rpr4[k]=12.000; k++;
rpr4[k]=10.320; k++; rpr4[k]=10.000; k++; rpr4[k]=9.125;  k++; rpr4[k]=8.000;  k++; rpr4[k]=6.980;  k++;
rpr4[k]=6.000;  k++; rpr4[k]=4.917;  k++; rpr4[k]=4.000;  k++; rpr4[k]=3.125;  k++; rpr4[k]=2.000;  k++;
rpr4[k]=1.423;  k++; rpr4[k]=1.115;  k++; rpr4[k]=0.846;  k++; rpr4[k]=0.615;  k++; rpr4[k]=0.471;  k++;
rpr4[k]=0.392;  k++; rpr4[k]=0.314;  k++; rpr4[k]=0.275;
return rpr4; }
double *LegrRasPorpoRazmkvi400()
{ int k=0; double *rp4=new double[kv4];
if (!rp4) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rp4[k]=0.053;  k++; rp4[k]=0.631;   k++; rp4[k]=0.774;   k++; rp4[k]=1.000;   k++; 
rp4[k]=1.442;  k++; rp4[k]=5.374;   k++; rp4[k]=6.610;   k++; rp4[k]=10.000;  k++; 
rp4[k]=13.092; k++; rp4[k]=16.408;  k++; rp4[k]=20.347;  k++; rp4[k]=26.640;  k++; 
rp4[k]=32.407; k++; rp4[k]=43.481;  k++; rp4[k]=59.785;  k++; rp4[k]=71.225;  k++; 
rp4[k]=84.243; k++; rp4[k]=100.00;  k++; rp4[k]=128.792; k++; rp4[k]=149.908; k++;
rp4[k]=193.07; k++; rp4[k]=213.634; k++; rp4[k]=261.567;
return rp4; }
double *rasprPorpoRazmAbskvi500()
{ int k=0; double *opr5=new double[kv5];
if (!opr5) { cout << "No memory!" << endl; k=getchar(); exit(1); }
opr5[k]=10.414; k++; opr5[k]=10.339; k++; opr5[k]=10.207; k++; opr5[k]=10.000; k++; opr5[k]=9.822; k++; 
opr5[k]=9.663;  k++; opr5[k]=9.465;  k++; opr5[k]=9.212;  k++; opr5[k]=9.059;  k++; opr5[k]=8.870; k++; 
opr5[k]=8.664;  k++; opr5[k]=8.458;  k++; opr5[k]=8.335;  k++; opr5[k]=8.190;  k++; opr5[k]=8.000; k++; 
opr5[k]=7.695;  k++; opr5[k]=7.477;  k++; opr5[k]=7.172;  k++; opr5[k]=6.859;  k++; opr5[k]=6.547; k++; 
opr5[k]=6.105;  k++; opr5[k]=6.000;  k++; opr5[k]=4.929;  k++; opr5[k]=4.469;  k++; opr5[k]=4.000; k++; 
opr5[k]=3.716;  k++; opr5[k]=3.406;  k++; opr5[k]=2.916;  k++; opr5[k]=2.699;  k++; opr5[k]=2.408; k++; 
opr5[k]=2.168;  k++; opr5[k]=2.000;  k++; opr5[k]=1.814;  k++; opr5[k]=1.616;  k++; opr5[k]=1.414; k++; 
opr5[k]=1.201;  k++; opr5[k]=0.935;  k++; opr5[k]=0.755;  k++; opr5[k]=0.667;  k++; opr5[k]=0.561; k++; 
opr5[k]=0.464;  k++; opr5[k]=0.360;  k++; opr5[k]=0.298;  k++; opr5[k]=0.286;  k++; opr5[k]=0.280; k++; 
opr5[k]=0.273;
return opr5; }
double *LegrRasPorpoRazmkvi500()
{ int k=0; double *rpr5=new double[kv5];
if (!rpr5) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr5[k]=0.05500; k++; rpr5[k]=0.44500; k++; rpr5[k]=0.55000; k++; rpr5[k]=0.64000; k++; rpr5[k]=0.71000; k++;
rpr5[k]=0.80200; k++; rpr5[k]=1.00000; k++; rpr5[k]=1.48200; k++; rpr5[k]=1.94400; k++; rpr5[k]=2.34800; k++;
rpr5[k]=2.86100; k++; rpr5[k]=3.52100; k++; rpr5[k]=3.86300; k++; rpr5[k]=4.66500; k++; rpr5[k]=4.94500; k++;
rpr5[k]=6.01500; k++; rpr5[k]=6.70000; k++; rpr5[k]=7.78600; k++; rpr5[k]=8.87300; k++; rpr5[k]=10.0000; k++;
rpr5[k]=11.0800; k++; rpr5[k]=11.2560; k++; rpr5[k]=11.9810; k++; rpr5[k]=13.1270; k++; rpr5[k]=14.4290; k++;
rpr5[k]=15.9240; k++; rpr5[k]=18.0650; k++; rpr5[k]=22.5290; k++; rpr5[k]=25.6180; k++; rpr5[k]=30.1050; k++;
rpr5[k]=33.8780; k++; rpr5[k]=37.3810; k++; rpr5[k]=41.7360; k++; rpr5[k]=46.4770; k++; rpr5[k]=53.0620; k++;
rpr5[k]=61.3300; k++; rpr5[k]=71.8070; k++; rpr5[k]=84.0730; k++; rpr5[k]=91.3160; k++; rpr5[k]=100.000; k++;
rpr5[k]=117.438; k++; rpr5[k]=137.916; k++; rpr5[k]=162.240; k++; rpr5[k]=190.638; k++; rpr5[k]=223.378; k++;
rpr5[k]=252.815;
return rpr5; }
double *rasprPorpoRazmAbskvi600()
{ int k=0; double *rpr6=new double[kv6];
if (!rpr6) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr6[k]=8.585; k++; rpr6[k]=8.558; k++; rpr6[k]=8.400; k++; rpr6[k]=8.246; k++; rpr6[k]=8.173; k++;
rpr6[k]=8.000; k++; rpr6[k]=7.849; k++; rpr6[k]=7.652; k++; rpr6[k]=7.410; k++; rpr6[k]=7.115; k++;
rpr6[k]=6.938; k++; rpr6[k]=6.721; k++; rpr6[k]=6.477; k++; rpr6[k]=6.222; k++; rpr6[k]=6.000; k++; 
rpr6[k]=5.886; k++; rpr6[k]=5.549; k++; rpr6[k]=5.338; k++; rpr6[k]=5.014; k++; rpr6[k]=4.593; k++; 
rpr6[k]=4.283; k++; rpr6[k]=4.000; k++; rpr6[k]=3.715; k++; rpr6[k]=3.178; k++; rpr6[k]=2.822; k++; 
rpr6[k]=2.362; k++; rpr6[k]=2.233; k++; rpr6[k]=2.000; k++; rpr6[k]=1.747; k++; rpr6[k]=1.531; k++; 
rpr6[k]=1.302; k++; rpr6[k]=1.093; k++; rpr6[k]=0.883; k++; rpr6[k]=0.739; k++; rpr6[k]=0.617; k++; 
rpr6[k]=0.528; k++; rpr6[k]=0.456; k++; rpr6[k]=0.313; k++; rpr6[k]=0.306; k++; rpr6[k]=0.275; k++; 
rpr6[k]=0.269; return rpr6; }
double *LegrRasPorpoRazmkvi600()
{ int k=0; double *rp6=new double[kv6];
if (!rp6) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rp6[k]=0.061;  k++; rp6[k]=0.606;  k++; rp6[k]=0.784;  k++; rp6[k]=0.921;   k++; rp6[k]=1.000;   k++; 
rp6[k]=1.238;  k++; rp6[k]=1.648;  k++; rp6[k]=2.296;  k++; rp6[k]=3.209;   k++; rp6[k]=4.477;   k++; 
rp6[k]=5.288;  k++; rp6[k]=6.247;  k++; rp6[k]=7.352;  k++; rp6[k]=8.682;   k++; rp6[k]=9.633;   k++; 
rp6[k]=10.321; k++; rp6[k]=11.711; k++; rp6[k]=12.674; k++; rp6[k]=14.260;  k++; rp6[k]=16.050;  k++; 
rp6[k]=17.383; k++; rp6[k]=18.013; k++; rp6[k]=18.738; k++; rp6[k]=19.882;  k++; rp6[k]=21.687;  k++; 
rp6[k]=24.608; k++; rp6[k]=25.803; k++; rp6[k]=27.377; k++; rp6[k]=31.436;  k++; rp6[k]=35.390;  k++; 
rp6[k]=41.447; k++; rp6[k]=48.541; k++; rp6[k]=56.848; k++; rp6[k]=66.578;  k++; rp6[k]=77.972;  k++; 
rp6[k]=84.381; k++; rp6[k]=100.00; k++; rp6[k]=123.24; k++; rp6[k]=149.563; k++; rp6[k]=190.423; k++; 
rp6[k]=252.405; return rp6; }
double *rasprPorpoRazmAbskvi700()
{ int k=0; double *rpr7=new double[kv7];
if (!rpr7) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr7[k]=5.678; k++; rpr7[k]=5.626; k++; rpr7[k]=5.579; k++; rpr7[k]=5.385; k++; rpr7[k]=5.142; k++;  
rpr7[k]=4.865; k++; rpr7[k]=4.709; k++; rpr7[k]=4.615; k++; rpr7[k]=4.561; k++; rpr7[k]=4.493; k++; 
rpr7[k]=4.446; k++; rpr7[k]=4.361; k++; rpr7[k]=4.242; k++; rpr7[k]=4.074; k++; rpr7[k]=4.000; k++; 
rpr7[k]=3.883; k++; rpr7[k]=3.440; k++; rpr7[k]=3.149; k++; rpr7[k]=2.864; k++; rpr7[k]=2.550; k++; 
rpr7[k]=2.259; k++; rpr7[k]=2.000; k++; rpr7[k]=1.904; k++; rpr7[k]=1.741; k++; rpr7[k]=1.617; k++;  
rpr7[k]=1.486; k++; rpr7[k]=1.381; k++; rpr7[k]=1.292; k++; rpr7[k]=1.205; k++; rpr7[k]=1.143; k++; 
rpr7[k]=1.043; k++; rpr7[k]=0.957; k++; rpr7[k]=0.876; k++; rpr7[k]=0.807; k++; rpr7[k]=0.725; k++;  
rpr7[k]=0.663; k++; rpr7[k]=0.594; k++; rpr7[k]=0.544; k++; rpr7[k]=0.469; k++; rpr7[k]=0.406; k++; 
rpr7[k]=0.384; k++; rpr7[k]=0.347; k++; rpr7[k]=0.285; k++; rpr7[k]=0.173; k++; rpr7[k]=0.149; k++; 
return rpr7; }
double *LegrRasPorpoRazmkvi700()
{ int k=0; double *rp7=new double[kv7];
if (!rp7) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rp7[k]=0.060;   k++; rp7[k]=1.000;   k++; rp7[k]=1.105;   k++; rp7[k]=1.452;   k++; rp7[k]=1.865;  k++;
rp7[k]=2.394;   k++; rp7[k]=2.947;   k++; rp7[k]=3.479;   k++; rp7[k]=4.109;   k++; rp7[k]=4.852;  k++;
rp7[k]=5.730;   k++; rp7[k]=6.766;   k++; rp7[k]=7.990;   k++; rp7[k]=9.435;   k++; rp7[k]=10.00;  k++;
rp7[k]=11.081;  k++; rp7[k]=13.716;  k++; rp7[k]=15.441;  k++; rp7[k]=17.367;  k++; rp7[k]=19.548; k++;
rp7[k]=22.002;  k++; rp7[k]=24.340;  k++; rp7[k]=25.761;  k++; rp7[k]=27.923;  k++; rp7[k]=30.219; k++;
rp7[k]=32.702;  k++; rp7[k]=35.390;  k++; rp7[k]=38.299;  k++; rp7[k]=41.447;  k++; rp7[k]=44.854; k++;
rp7[k]=48.541;  k++; rp7[k]=52.531;  k++; rp7[k]=56.848;  k++; rp7[k]=61.521;  k++; rp7[k]=66.578; k++;
rp7[k]=72.050;  k++; rp7[k]=77.972;  k++; rp7[k]=84.381;  k++; rp7[k]=91.316;  k++; rp7[k]=100.00; k++;
rp7[k]=107.542; k++; rp7[k]=117.537; k++; rp7[k]=138.150; k++; rp7[k]=176.041; k++; rp7[k]=233.572;
return rp7; }
double *rasprPorpoRazmAbskvi800()
{ int k=0; double *rpr8=new double[kv8];
if (!rpr8) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr8[k]=5.106; k++; rpr8[k]=5.078; k++; rpr8[k]=5.044; k++; rpr8[k]=4.935; k++;
rpr8[k]=4.792; k++; rpr8[k]=4.601; k++; rpr8[k]=4.348; k++; rpr8[k]=4.000; k++;
rpr8[k]=3.929; k++; rpr8[k]=3.811; k++; rpr8[k]=3.718; k++; rpr8[k]=3.686; k++;
rpr8[k]=3.663; k++; rpr8[k]=3.631; k++; rpr8[k]=3.608; k++; rpr8[k]=3.576; k++;
rpr8[k]=3.545; k++; rpr8[k]=3.467; k++; rpr8[k]=3.357; k++; rpr8[k]=3.224; k++;
rpr8[k]=3.055; k++; rpr8[k]=2.977; k++; rpr8[k]=2.852; k++; rpr8[k]=2.779; k++;
rpr8[k]=2.669; k++; rpr8[k]=2.565; k++; rpr8[k]=2.468; k++; rpr8[k]=2.383; k++;
rpr8[k]=2.279; k++; rpr8[k]=2.156; k++; rpr8[k]=2.084; k++; rpr8[k]=2.000; k++;
rpr8[k]=1.846; k++; rpr8[k]=1.611; k++; rpr8[k]=1.481; k++; rpr8[k]=1.309; k++;
rpr8[k]=1.173; k++; rpr8[k]=1.074; k++; rpr8[k]=0.932; k++; rpr8[k]=0.765; k++;
rpr8[k]=0.630; k++; rpr8[k]=0.568; k++; rpr8[k]=0.457; k++; rpr8[k]=0.358; k++;
rpr8[k]=0.309; k++; rpr8[k]=0.272; k++; rpr8[k]=0.225; k++; rpr8[k]=0.210; k++;
rpr8[k]=0.180; k++; rpr8[k]=0.172; k++; rpr8[k]=0.165; k++; rpr8[k]=0.157; k++;
rpr8[k]=0.150; k++; rpr8[k]=0.142;
return rpr8; }
double *LegrRasPorpoRazmkvi800()
{ int k=0; double *rp8=new double[kv8];
if (!rp8) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rp8[k]=0.061;  k++; rp8[k]=0.348;  k++; rp8[k]=0.417;  k++; rp8[k]=0.471;  k++;
rp8[k]=0.554;  k++; rp8[k]=0.648;  k++; rp8[k]=0.762;  k++; rp8[k]=0.952;  k++;
rp8[k]=1.000;  k++; rp8[k]=1.161;  k++; rp8[k]=1.348;  k++; rp8[k]=1.489;  k++;
rp8[k]=2.216;  k++; rp8[k]=2.704;  k++; rp8[k]=2.986;  k++; rp8[k]=3.299;  k++;
rp8[k]=3.644;  k++; rp8[k]=4.446;  k++; rp8[k]=5.424;  k++; rp8[k]=6.618;  k++;
rp8[k]=8.075;  k++; rp8[k]=8.919;  k++; rp8[k]=10.00;  k++; rp8[k]=10.82;  k++;
rp8[k]=11.708; k++; rp8[k]=12.669; k++; rp8[k]=13.708; k++; rp8[k]=14.833; k++;
rp8[k]=16.050; k++; rp8[k]=17.367; k++; rp8[k]=18.792; k++; rp8[k]=19.242; k++;
rp8[k]=22.002; k++; rp8[k]=25.761; k++; rp8[k]=27.874; k++; rp8[k]=30.161; k++;
rp8[k]=32.636; k++; rp8[k]=35.314; k++; rp8[k]=41.346; k++; rp8[k]=48.410; k++;
rp8[k]=56.679; k++; rp8[k]=61.330; k++; rp8[k]=71.807; k++; rp8[k]=84.073; k++;
rp8[k]=90.971; k++; rp8[k]=100.00; k++; rp8[k]=110,181;k++; rp8[k]=121.398;k++;
rp8[k]=133.757;k++; rp8[k]=147.374;k++; rp8[k]=162.378;k++; rp8[k]=178.909;k++;
rp8[k]=197.123;k++; rp8[k]=217.191;
return rp8; }
double *rasprPorpoRazmAbskvi900()
{ int k=0; double *rpr9=new double[kv9];
if (!rpr9) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr9[k]=4.374; k++; rpr9[k]=4.268; k++; rpr9[k]=4.236; k++; rpr9[k]=4.140; k++;
rpr9[k]=4.000; k++; rpr9[k]=3.868; k++; rpr9[k]=3.660; k++; rpr9[k]=3.538; k++;
rpr9[k]=3.349; k++; rpr9[k]=3.236; k++; rpr9[k]=3.094; k++; rpr9[k]=2.934; k++;
rpr9[k]=2.841; k++; rpr9[k]=2.710; k++; rpr9[k]=2.516; k++; rpr9[k]=2.302; k++;
rpr9[k]=2.151; k++; rpr9[k]=2.000; k++; rpr9[k]=1.778; k++; rpr9[k]=1.591; k++;
rpr9[k]=1.404; k++; rpr9[k]=1.218; k++; rpr9[k]=1.076; k++; rpr9[k]=0.932; k++; 
rpr9[k]=0.765; k++; rpr9[k]=0.630; k++; rpr9[k]=0.568; k++; rpr9[k]=0.457; k++; 
rpr9[k]=0.358; k++; rpr9[k]=0.309; k++; rpr9[k]=0.272; k++; rpr9[k]=0.225; k++; 
rpr9[k]=0.210; k++; rpr9[k]=0.180; k++; rpr9[k]=0.172; k++; rpr9[k]=0.165; k++; 
rpr9[k]=0.157; k++; rpr9[k]=0.150; k++; rpr9[k]=0.142;
return rpr9; }
double *LegrRasPorpoRazmkvi900()
{ int k=0; double *rp9=new double[kv9];
if (!rp9) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rp9[k]=0.063;  k++; rp9[k]=0.823;  k++; rp9[k]=1.000;  k++; rp9[k]=1.127;  k++;
rp9[k]=1.254;  k++; rp9[k]=1.430;  k++; rp9[k]=1.711;  k++; rp9[k]=2.046;  k++;
rp9[k]=2.926;  k++; rp9[k]=3.715;  k++; rp9[k]=4.716;  k++; rp9[k]=5.987;  k++;
rp9[k]=6.746;  k++; rp9[k]=8.067;  k++; rp9[k]=10.000; k++; rp9[k]=12.546; k++;
rp9[k]=14.053; k++; rp9[k]=15.759; k++; rp9[k]=18.690; k++; rp9[k]=22.166; k++;
rp9[k]=26.288; k++; rp9[k]=31.176; k++; rp9[k]=34.931; k++; rp9[k]=41.346; k++; 
rp9[k]=48.410; k++; rp9[k]=56.679; k++; rp9[k]=61.330; k++; rp9[k]=71.807; k++; 
rp9[k]=84.073; k++; rp9[k]=90.971; k++; rp9[k]=100.00; k++; rp9[k]=110,181;k++; 
rp9[k]=121.398;k++; rp9[k]=133.757;k++; rp9[k]=147.374;k++; rp9[k]=162.378;k++; 
rp9[k]=178.909;k++; rp9[k]=197.123;k++; rp9[k]=217.191;
return rp9; }
double *rasprPorpoRazmAbskvi1000()
{ int k=0; double *rpr10=new double[kv10];
if (!rpr10) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr10[k]=3.607; k++; rpr10[k]=3.549; k++; rpr10[k]=3.446; k++;
rpr10[k]=3.066; k++; rpr10[k]=2.616; k++; rpr10[k]=2.151; k++;
rpr10[k]=2.047; k++; rpr10[k]=1.663; k++; rpr10[k]=1.328; k++;
rpr10[k]=1.092; k++; rpr10[k]=0.893; k++; rpr10[k]=0.713; k++;
rpr10[k]=0.515; k++; rpr10[k]=0.353; k++; rpr10[k]=0.250; k++;
rpr10[k]=0.140; k++; rpr10[k]=0.059;
return rpr10; }
double *LegrRasPorpoRazmkvi1000()
{ int k=0; double *rp10=new double[kv10];
if (!rp10) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rp10[k]=0.063;  k++; rp10[k]=0.820;  k++; rp10[k]=1.000;  k++;
rp10[k]=2.046;  k++; rp10[k]=4.186;  k++; rp10[k]=8.469;  k++;
rp10[k]=10.00;  k++; rp10[k]=14.609; k++; rp10[k]=19.412; k++;
rp10[k]=23.462; k++; rp10[k]=28.358; k++; rp10[k]=34.275; k++;
rp10[k]=41.427; k++; rp10[k]=50.071; k++; rp10[k]=60.519; k++;
rp10[k]=73.147; k++; rp10[k]=88.410;
return rp10; }