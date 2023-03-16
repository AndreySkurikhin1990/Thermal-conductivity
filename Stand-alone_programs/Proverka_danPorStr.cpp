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
const double pi = acos(-1e0), tem0=273.15, prk=3e2, pri=16e1;
const int dsnftk=50, dmvtk=28, cvym=6, dmsrpv=18;
struct raspporporazmm { double rprm; struct raspporporazmm *slel; struct raspporporazmm *prel; };
double **poisrasprpor1(double, double, int, double, double, int, double *, double *, double *, double **);
double **poisrasprpor2(double, double, int, double, double, int, double *, double *, double *, double **);
double **rasppor2(double, double, double *, double *, double *, int, double, double, int, double **);
double *rasPorpoRazmVerI207(double *, int);
double *rasPorpoRazmVerIs84(double *, int);
double *rasPorpoRazmVerIn84(double *, int);
double *rasPorpoRazmVerI16035(double *, int);
double *rasPorpoRazmVerO84(double *, int);
double *rasPorpoRazmVerO16035(double *, int);
double *NapMasRaspVer(int, int, int, int, double *);
double **NapMasRasLePrVer(int, int, int, int, double **);
double **rasPorpoRazVer(double, int, int, int, int, int);
double **rasrasporporaz1(double *, double, int, int, double *, double *, double *, double, double, double, int, double **);
double **rasrasporporaz0(double *, double *, double *, double, int, int, double, double **);
double *rasPorpoRazmitom440(double *, int);
double *rasPorpoRazmitom620(double *, int);
double *rasPorpoRazmitom860(double *, int);
double *rasPorpoRazmitom1000(double *, int);
double *PravGranPoromVer(double *, int);
double *LevGranPoromVer(double *, double *, int);
double **vybRasPorRazVer(double *, int, int, int, int);
double novNapMas(int, int, int, int, int, int, int, int, int);
//----------------
double novNapMas(int vyve, int vymave, int vmimfv, int vyfrve, int vyukve, int vysove, int vpmf, int vpkf, int cemv)
{	double por=0.0, ko=1e-2;
	if (!vyve) //если выбран шамот
	{
	double porsha0=21.8*ko, porsha1=11.014*ko, porsha2=25.2*ko, porsha3=26.5*ko, porsha4=11.5*ko, porsha5=16.5*ko; //0 - Роучка, 1 - ШПД-41, 2 - ШБ-1 2-1, 3 - ШВ-1 1-1, 4 - ШПД, 5 - ШКУ-32 3-1
	if (!vymave) por=porsha0; else if (vymave==1) por=porsha1; else if (vymave==2) por=porsha2;
	else if (vymave==3) por=porsha3; else if (vymave==4) por=porsha4; else if (vymave==5) por=porsha5;
	}
	else if (vyve==1) 
	{ //если выбран вермикулит
	double por207=66.35*ko, poris84=55.75*ko, porin84=81.53*ko;
	double poro84=86.61*ko, poro16035=84.36*ko, por16035=83.97*ko;
	if ((!vysove) || (vysove==1)) 
	{  //исходный или после повторных измерений
	if (!vyfrve) 
	{ 
	if (!vpmf) por = por207; 
	else if (vpmf==1) por=por16035; 
	} //выбор пористости мелкой фракции
	else if (vyfrve==1) 
	{ 
	if (!vpkf) por = poris84; 
	else if (vpkf==1) por=porin84; 
	} //выбор пористости крупной фракции
	else if (vyfrve==2) por = poro16035; 
	}
	if (vysove==2) 
	{ //после обжига
	if ((!vyfrve) || (vyfrve==2)) por=poro16035;
	else if (vyfrve==1) por=poro84; 
	} 
	}
	else if (vyve==2) 
	{ //если выбран ИТОМ
	double pori440=(8e1+82e0)*ko/2e0, pori620=(75e0+78e0)*ko/2e0;
	double pori860=(65e0+68e0)*ko/2, pori1000=(62e0+65e0)*ko/2e0;
	if (!vymave) por=pori440; 
	else if (vymave==1) por=pori620; 
	else if (vymave==2) por=pori860; 
	else if (vymave==3) por=pori1000;
	}
	else if (vyve==3) 
	{ //если выбран КВИ
	double pork400=51.74*ko, pork500=52.07*ko, pork600=51.51*ko, pork700=39.75*ko;
	double pork800=40.85*ko, pork900=39.37*ko, pork1000=36.07*ko;
	if (vymave==4) por=pork400; 
	else if (vymave==5) por=pork500; 
	else if (vymave==6) por=pork600; 
	else if (vymave==7) por=pork700; 
	else if (vymave==8) por=pork800; 
	else if (vymave==9) por=pork900;
	else if (vymave==10) por=pork1000;
	}
	return por;
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
{ int k=0, n=0, j=0, qg=0, nm=cvym;
double *legr=NULL, *prgr=NULL, *raspr=NULL, *srra=NULL, ht=1e0, s=0.0, e=1e-1;
double rpn=1e0, dp=1e0, ras0=0.0, marp=0.0, srp=0.0, **mu=NULL, *po=NULL, r=0.0;
mu=new double*[nm]; if (!mu) { cout << "No memory" << endl; k=getchar(); exit(1); } 
mu=NapMasRasLePrVer(vfv, vysove, vpkf, k, mu);
k=0; raspr=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++;
po=mu[k]; j=0; srp=po[j]; k++; if (po) delete[]po; 
po=mu[k]; marp=po[j]; if (po) delete[]po; 
k=0; s=e; while (s<marp) { s=s+ht; k++; } n=k; //cout << "srp = " << srp << "\tmarp = " << marp << "\tn = " << n << endl; //for (k=0; k<n; k++) cout << "rpr ( " << k << " ) = " << raspr[k] << "\tsrra ( " << k << " ) = " << srra[k] << endl;
k=0; mu=rasrasporporaz0(prgr, legr, raspr, rpn, k, n, dp, mu);
s=0.0; double l=s, p=1e-6, h=p, ko=ht/h;
if (legr) delete []legr; if (prgr) delete []prgr; 
if (raspr) delete[]raspr; if (srra) delete[]srra;
k=0; raspr=mu[k]; k++; srra=mu[k]; k++; prgr=mu[k]; k++; legr=mu[k]; k++; 
po=mu[k]; j=0; srp=po[j]; k++; if (po) delete[]po;
po=mu[k]; marp=po[j]; if (po) delete[]po; 
k=0; s=e; while (s<marp) { s=s+ht; k++; } qg=k; cout << "qg = " << qg << "\tsrp = " << srp << "\tisrp = " << isrp << endl;
k=1; double *ms=new double[k], *mm=new double[k];
if (!isrp) { s=0.0; r=0.0;
for (k=0; k<qg; k++) { srra[k]=(p+l)/2e0; //расчет характеристик пористой структуры при условии пор в форме прямоугольного параллелепипеда
s=s+srra[k]*raspr[k]; r=r+raspr[k]; legr[k]=l; p=p+h; l=l+h; } 
k=0; if (fabs(r)>0.0) s=s/r; else r=0.0; } 
else if (isrp==1) { //расчет характеристик пористой структуры при условии шарообразных пор
double sz=0.0, sch=sz, dv=sz, vpre=sz, vtek=sz, rtek=sz, rpre=sz; l=0.0; p=1e-6;
for (k=0; k<qg; k++) { rtek=(p+l)/2e0; srra[k]=rtek; 
vtek=(pi/6e0)*pow(rtek,3e0); vpre=(pi/6e0)*pow(rpre,3e0); dv=fabs(vtek-vpre);
sz=sz+raspr[k]; sch=sch+raspr[k]*dv; rpre=rtek;
legr[k]=l; p=p+h; l=l+h; } s=sch/sz; } 
else { cout << "Nepravilno vybran nomer" << endl; k=getchar(); exit(1); }
k=0; ms[k]=s; l=l*ko; mm[k]=l;
cout << "Sr raz por = " << s << "\tMax Raz Por = " << l << endl;
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
{ struct raspporporazmm *rpr01, *rprp01, *roo=new struct raspporporazmm, *totr;
double s=0.0, e=1e-1, r=0.0, koef=1e6, *prgr00=new double[n];
double *legr00=new double[n], m=0.0, t=0.0, ht=1e0; 
int k=0, q=0, p=0, j=0, jk=0, qg=0;
for (k=0; k<n; k++) prgr00[k]=prgr0[k]*koef;
legr00[0]=0.0; for (k=1; k<n; k++) legr00[k]=prgr00[k-1]; //for (k=0; k<n; k++) cout << "prgr ( " << k << " ) = " << prgr00[k] << "\tras = " << raspr0[k] << "\t"; 
s=0.0; for (k=0; k<n; k++) if (prgr00[k]<(ht+e)) { p=k; s=s+raspr0[k]; } //cout << "\np = " << p << "\ts = " << s << "\t"; //размер пор до 1 мкм
roo->rprm=s; roo->slel=NULL; rpr01=roo; q=0; rprp01=roo; roo->prel=NULL;
for (k=p; k<n; k++) { 
    r=prgr00[k]-legr00[k];
    if ((prgr00[k]>(ht+e)) && (r<(ht+e)) && (r>e)) {
		rpr01=new struct raspporporazmm; rpr01->prel=rprp01; rpr01->rprm=raspr0[k]; 
		rprp01->slel=rpr01; rpr01->slel=NULL; rprp01=rpr01; cout << "q_1 = " << q << endl; q++; }
    if (r>(ht+e)) {
        m=raspr0[k]/r;
		t=r; j=0; while (t>e) { j++; t=t-ht; } jk=j; 
        for (j=0; j<jk; j++) {
		rpr01=new struct raspporporazmm;
		rpr01->prel=rprp01;
		rpr01->rprm=m;
		rprp01->slel=rpr01; 
		rpr01->slel=NULL;
		rprp01=rpr01;
		q++; } } }
qg=q; totr=rpr01; double *ras01=new double[qg], *prgrm=new double[qg], *legrm=new double[qg];
double prgr01=rpn, *srra01=new double[qg]; cout << "qg = " << qg << "\trpn = " << rpn << "\tde = " << de << endl;
q=1; double *mm=new double[q], *ms=new double[q]; s=0.0; r=0.0;
if ((!ras01) || (!prgrm) || (!legrm) || (!srra01) || (!ms) || (!mm))
{ cout << "No memory!" << endl; k=getchar(); exit(1); }
k=0; legrm[k]=0.0; 
while ((rpr01) && (k<qg)) 
{ prgrm[k]=prgr01; cout << "\tras ( " << k << " ) = " << rpr01->rprm;
if (k>0) legrm[k]=prgrm[k-1];
ras01[k]=rpr01->rprm;
rprp01=rpr01;
rpr01=rpr01->prel;
if (rprp01) delete rprp01;
k++;
prgr01=prgr01+de; }
if (prgr00) delete []prgr00; if (legr00) delete []legr00;
t=0.0; for (k=0; k<qg; k++) t=t+ht; k=0; mm[k]=t;
for (k=0; k<qg; k++) srra01[k]=(prgrm[k]+legrm[k])/2e0;
s=0.0; r=s; for (k=0; k<qg; k++) { s=s+ras01[k]*srra01[k]; r=r+ras01[k]; }
if (fabs(r)>0.0) s=s/r; else s=0.0; k=0; ms[k]=s;
k=0; mu[k]=ras01; k++; mu[k]=srra01; k++; mu[k]=prgrm; k++; mu[k]=legrm; k++; 
mu[k]=ms; k++; mu[k]=mm; for (k=0; k<qg; k++) cout << "rp ( " << k << " ) = " << ras01[k] << "\t";
return mu; }