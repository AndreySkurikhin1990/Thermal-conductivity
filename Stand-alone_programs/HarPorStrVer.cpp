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
const int dsnftk=50, dmvtk=28, dmsrpv=18;
struct raspporporazmm { double rprm; struct raspporporazmm *slel; };
double *prgr=NULL, *raspr=NULL, *legr=NULL, *srra=NULL, *rpr=NULL; 
int qg;
double *poisrasprpor1(double, double, int, double, double, int);
double *poisrasprpor2(double, double, int, double, double, int);
double *poisrasprpor0(int, double, double, int);
double *rasppor1(double, double, double *, double *, double *, int, double, double, int);
double *rasppor2(double, double, double *, double *, double *, int, double, double, int);
double *rasPorpoRazmVerI207(double *, int);
double *rasPorpoRazmVerIs84(double *, int);
double *rasPorpoRazmVerIn84(double *, int);
double *rasPorpoRazmVerI16035(double *, int);
double *rasPorpoRazmVerO84(double *, int);
double *rasPorpoRazmVerO16035(double *, int);
void NapMasRaspVer(int, int, int, int);
double *rasrasporporaz1(double *, double, int, int, double *, double *, double *, double, double, double, int);
double *rasrasporporaz0(double *, double *, double *, double, int, int, double);
double *PravGranPoromVer(double *, int);
double *LevGranPoromVer(double *, double *, int);
double srRaPorVerm(int, int, int, int);
//-------------Вермикулит----------------------
double srRaPorVerm(int n, int vfv, int vysove, int vpkf)
{ int k; double s; cout << "fr = " << vfv << "\tvsv = " << vysove << "\tvpkf = " << vpkf << endl;
prgr=PravGranPoromVer(prgr, n); legr=LevGranPoromVer(prgr, legr, n);
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2e0; //средний размер пор
if ((!vysove) || (vysove==2)) { 
if (!vfv) raspr=rasPorpoRazmVerI207(raspr, dmsrpv); //для фракции 2-0,7 мм, исходный
else if (vfv==1) { if (!vpkf) raspr=rasPorpoRazmVerIs84(raspr, dmsrpv); 
else if (vpkf==1) raspr=rasPorpoRazmVerIn84(raspr, dmsrpv); } //для фракции 8-4 мм, исходный
else if (vfv==2) raspr=rasPorpoRazmVerI16035(raspr, dmsrpv); } //для фракции 1,6-0,35 мм, исходный
if (vysove==1) {
if ((vfv==2) || (!vfv)) raspr=rasPorpoRazmVerO16035(raspr, dmsrpv); //для фракции 1,6-0,35 мм и 2-0,7 мм, после обжига
else if (vfv==1) raspr=rasPorpoRazmVerO84(raspr, dmsrpv); //для фракции 8-4 мм, после обжига
else { cout << "Net takoy fraktsii!"; k=getchar(); exit(1); } } //только для фракции 8-4 мм
s=0.0; for (k=0; k<n; k++) s=s+srra[k]*raspr[k];
/*cout << "Sredniy razmer por vermikulita: " << s << endl;*/ return s; }
double *PravGranPoromVer(double *rapo, int n)
{ int k=0; 
rapo[k]=13e1; k++; rapo[k]=12e1; k++; rapo[k]=11e1; k++; rapo[k]=1e2; k++; rapo[k]=9e1; k++;
rapo[k]=8e1;  k++; rapo[k]=7e1;  k++; rapo[k]=6e1;  k++; rapo[k]=5e1; k++; rapo[k]=4e1; k++;
rapo[k]=3e1;  k++; rapo[k]=2e1;  k++; rapo[k]=15.0; k++; rapo[k]=1e1; k++; rapo[k]=5.0; k++;
rapo[k]=3.0;  k++; rapo[k]=1.0;  k++; rapo[k]=1e-1; k++;
for (k=0; k<n; k++) rapo[k]=(1e-6)*rapo[k]; 
int j=(n-(n%2))/2; double t; for (k=0; k<j; k++) { t=rapo[k]; rapo[k]=rapo[n-1-k]; rapo[n-1-k]=t; } 
return rapo; }
double *LevGranPoromVer(double *rapo, double *levgr, int n)
{ int k; levgr[0]=(1e-6)*(1e-2); for (k=1; k<n; k++) levgr[k]=rapo[k-1]; return levgr; }
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
void NapMasRaspVer(int vfv, int n, int vysove, int vpkf)
{ int k=0;
if (!vysove) { if (!vfv) raspr=rasPorpoRazmVerI207(raspr, n);
else if (vfv==1) { if (!vpkf) raspr=rasPorpoRazmVerIs84(raspr, n); 
else if (vpkf==1) raspr=rasPorpoRazmVerIn84(raspr, n); }
else if (vfv==2) raspr=rasPorpoRazmVerI16035(raspr, n); }
if (vysove==1) { if (vfv==1) raspr=rasPorpoRazmVerO84(raspr, n); 
if (vfv==2) raspr=rasPorpoRazmVerO16035(raspr, n); } }
double *NapMasRasLePrVer(int vfv, int vysove, int vpkf, int nom)
{ int k, n=dmsrpv, q; 
legr=new double[n]; prgr=new double[n]; raspr=new double[n]; srra=new double[n]; rpr=new double[n]; 
if ((!legr) || (!prgr) || (!raspr) || (!srra) || (!rpr)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) { raspr[k]=0.0; prgr[k]=0.0; legr[k]=0.0; srra[k]=0.0; rpr[k]=0.0; }
NapMasRaspVer(vfv, n, vysove, vpkf); for (k=0; k<n; k++) if (raspr[k]>0.0) rpr[k]=raspr[k];
prgr=PravGranPoromVer(prgr, n); legr=LevGranPoromVer(prgr, legr, n);
for (k=0; k<n; k++) srra[k]=(legr[k]+prgr[k])/2e0;
if (!nom) { delete[]legr; delete[]prgr; delete[]raspr; delete[]rpr; return srra; }
else if (nom==1) { delete[]legr; delete[]prgr; delete[]srra; delete[]rpr; return raspr; } 
else if (nom==2) { k=1; double m=0.0, *mm=new double[k]; 
if (!mm) { cout << "No memory!" << endl; k=getchar(); exit(1); }
for (k=0; k<n; k++) if (m<prgr[k]) m=prgr[k]; k=0; mm[k]=m; return mm; } } //средний размер каждого из диапазонов
double *rasPorpoRazVer(double poris, int vyve, int vyvyzn, int vysove, int isrp, int vpkf) //0 - старые, 1 - новые значения
{ int n=dmsrpv, k;
legr=new double[n]; prgr=new double[n]; rpr=new double[n];
double rpn=1e0, dp=1e0;
if ((!rpr) || (!prgr) || (!legr)) { cout << "No memory!" << endl; n=getchar(); exit(1); }
prgr=PravGranPoromVer(prgr, n); legr=LevGranPoromVer(prgr, legr, n);
if ((!vysove) || (vysove==1)) {
if (!vyve) rpr=rasPorpoRazmVerI207(rpr,n);
else if (vyve==1) { if (!vpkf) rpr=rasPorpoRazmVerIs84(rpr,n); else if (vpkf==1) rpr=rasPorpoRazmVerIn84(rpr,n); }
else if (vyve==2) rpr=rasPorpoRazmVerI16035(rpr,n); }
if (vysove==2) { if ((!vyve) || (vyve==2)) rpr=rasPorpoRazmVerO16035(rpr,n);
else if (vyve==1) rpr=rasPorpoRazmVerO84(rpr,n); }
double *rp=poisrasprpor0(n, rpn, dp, 0), s=0.0, l=0.0, p=1e-6, h=1e-6, *m=new double[1];
srra=new double[qg]; delete []legr; delete []prgr; delete []rpr;
legr=new double[qg]; 
if (!isrp) { for (k=0; k<qg; k++) { srra[k]=(p+l)/2e0; 
s=s+srra[k]*rp[k]; legr[k]=l; p=p+h; l=l+h; } } 
if (isrp==1) { double sz=0.0, sch=0.0, dv, vpre=0.0, vtek=0.0, rtek=0.0, rpre=0.0; l=0.0; p=1e-6;
for (k=0; k<qg; k++) { rtek=(p+l)/2e0; srra[k]=rtek; 
vtek=(pi/6e0)*pow(rtek,3e0); vpre=(pi/6e0)*pow(rpre,3e0); dv=fabs(vtek-vpre);
sz=sz+rp[k]*dv; sch=sch+rtek*rp[k]*dv; rpre=rtek;
legr[k]=l; p=p+h; l=l+h; } s=sch/sz; } //cout << "Sr raz por = " << s << endl;
if (!vyvyzn) { delete []srra; delete []legr; delete []m; return rp; }
else if (vyvyzn==1) { delete []rp; delete []srra; delete []legr; m[0]=s; return m; }
else if (vyvyzn==2) { delete []rp; delete []srra; delete []legr; m[0]=l*1e6; return m; }
else if (vyvyzn==3) { delete []rp; delete []legr; delete []m; return srra; }
else if (vyvyzn==4) { delete []rp; delete []srra; delete []m; return legr; } }
double *rasppor1(double p1, double p0, double *legr01, double *prgr01, double *ras01, int nn, double rpn, double d, int vlpr)
{ int k, n=nn;
double *ras21=rasrasporporaz1(ras01, rpn, 0, n, ras01, prgr01, legr01, p1, p0, d, vlpr);
return ras21; }
double *rasrasporporaz1(double *raspr0, double rpn, int vyb, int nn, double *ras11, double *prgr01, double *legr01, double p1, double p0, double de, int vlpr)
{ 	int k, q, j, n=nn; 
	double *legr11=new double[n], *prgr11=new double[n], kot, dob1, dob2, dob3, xp, kop, xt, s;
	if ((!legr11) || (!prgr11)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) prgr11[k]=prgr01[k]*p1/p0; legr11[0]=0.0; for (k=1; k<n; k++) legr11[k]=prgr11[k-1]; //for (k=0; k<n; k++) if (k<50) cout << "k = " << k << "\tr = " << ras11[k] << endl;
	struct raspporporazmm *rpr21, *rprp21, *roo=new struct raspporporazmm;
	rpr21=roo; q=0; xp=0.0; xt=0.0;
	for (k=0; k<n; k++) {
    xt=prgr01[k];
    for (j=0; j<n; j++) {
        if ((legr11[j]<xt) && (prgr11[j]>xt)) {
        kot=ras11[j]/fabs(prgr11[j]-legr11[j]);
        dob1=kot*(xt-legr11[j]);
		rpr21->rprm=dob1;
		        if (xp<legr11[j]) {
					if (j) {
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
double *ras21=new double[q], *prgrm=new double[q], *legrm=new double[q], prgr001=rpn, t, *srra21=new double[q], *tm=new double[1];
if ((!ras21) || (!prgrm) || (!legrm) || (!srra21)) { cout << "No memory!" << endl; k=getchar(); exit(1); }
rpr21=roo; rprp21=roo; k=0; legrm[0]=0.0; 
while ((rpr21) && (k<q)) { 
	prgrm[k]=prgr001;
	if (k) legrm[k]=prgrm[k-1];
	ras21[k]=rprp21->rprm; 
	rpr21=rprp21->slel; 
	delete rprp21;
	rprp21=rpr21;
	prgr001=prgr001+de; k++; }
q=k; s=0.0; for (k=0; k<q; k++) { s=s+ras21[k]; srra21[k]=(legrm[k]+prgrm[k])/2e0; }
t=1e2; for (k=0; k<q; k++) ras21[k]=ras21[k]*t/s;
s=1e0; t=0.0; for (k=0; k<q; k++) t=t+s; //for (k=0; k<q; k++) if (k<nk) printf("r = %0.4lf\tpr = %0.0lf\n",ras21[k],prgrm[k]);
delete []prgr11; delete []legr11; delete []legrm; delete []prgrm;
if (!vlpr) { delete []srra21; delete []tm; return ras21; } 
if (vlpr==1) { delete []ras21; delete []tm; return srra21; }
if (vlpr==2) { delete []ras21; delete []srra21; tm[0]=t; return tm; } }
double *poisrasprpor0(int nn, double ng, double shag, int v)
{ int n=nn, k, q, j, p; 
double *ras01=rasrasporporaz0(prgr, legr, rpr, ng, v, n, shag);
return ras01; }
double *poisrasprpor1(double p1, double p0, int nn, double rpn, double dp, int vlpr)
{ int n=nn, k, nk;
double *ras01=poisrasprpor0(n, rpn, dp, 0);
double *prgr01=poisrasprpor0(n, rpn, dp, 1);
double *legr01=poisrasprpor0(n, rpn, dp, 2);
double *ras21=rasppor1(p1, p0, legr01, prgr01, ras01, qg, rpn, dp, 0);
double *srra21=rasppor1(p1, p0, legr01, prgr01, ras01, qg, rpn, dp, 1);
k=1; double *l=rasppor1(p1, p0, legr01, prgr01, ras01, qg, rpn, dp, 2), d=l[0], p=d, hf=1e0, *m=new double[k]; delete []l;
k=0; while (d>0.0) { k++; d=d-hf; } nk=k; 
d=0.0; for (k=0; k<nk; k++) { srra21[k]=srra21[k]*1e-6; ras21[k]=ras21[k]*1e-2; legr01[k]=legr01[k]*1e-6; d=d+ras21[k]*srra21[k]; } //cout << "p = " << p << "\td = " << d << endl;
delete []ras01; delete []prgr01;
if (!vlpr) { delete []srra21; delete []legr01; delete []m; return ras21; }
else if (vlpr==1) { delete []ras21; delete []legr01; delete []m; return srra21; }
else if (vlpr==2) { delete []ras21; delete []srra21; delete []legr01; m[0]=d; return m; }
else if (vlpr==3) { delete []ras21; delete []srra21; delete []m; return legr01; } 
else if (vlpr==4) { delete []ras21; delete []srra21; delete []legr01; m[0]=p; return m; } }
double *rasppor2(double p2, double p0, double *legr01, double *prgr01, double *ras21, int nn, double rpn, double de, int vlpr)
{ int k, n=nn; 
double *legr21=new double[n], *prgr21=new double[n];
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
            if (xp<=legr21[j]) { if (j) {
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
rpr22=roo; rprp22=roo; k=0; legrm[0]=0.0; n=q; 
while ((rpr22) && (k<q)) { 
	prgrm[k]=prgr002;
	if (k>0) legrm[k]=prgrm[k-1];
	ras22[k]=rprp22->rprm; 
	rpr22=rprp22->slel; 
	delete rprp22;
	rprp22=rpr22;
	prgr002=prgr002+de; k++; } 
s=0.0; for (k=0; k<n; k++) { s=s+ras22[k]; srra22[k]=(legrm[k]+prgrm[k])/2.0; } 
t=1e2; for (k=0; k<n; k++) ras22[k]=ras22[k]*t/s; 
t=0.0; s=1e0; for (k=0; k<n; k++) t=t+s; //for (k=0; k<q; k++) printf("r22 = %0.4lf\tpr = %0.1lf\tr21 = %0.4lf\n",ras22[k],srra22[k],ras21[k]);
delete []legr21; delete []prgr21; delete []legrm; delete []prgrm;
if (!vlpr) { delete []srra22; return ras22; } 
else if (vlpr==1) { delete []ras22; return srra22; }
else if (vlpr==2) { delete []ras22; delete []srra22; return &t; } }
double *poisrasprpor2(double p2, double p0, int nn, double rpn, double h, int vlpr)
{ int n=nn, k, nk;
double *ras01=poisrasprpor0(n, rpn, h, 0);
double *prgr01=poisrasprpor0(n, rpn, h, 1);
double *legr01=poisrasprpor0(n, rpn, h, 2);
double *ras22=rasppor2(p2,p0,legr01,prgr01,ras01,qg,rpn,h,0);
double *srra2=rasppor2(p2,p0,legr01,prgr01,ras01,qg,rpn,h,1);
double *l=rasppor2(p2,p0,legr01,prgr01,ras01,qg,rpn,h,2), d=l[0], p=d, kor=1e-2, kos=1e-6, hf=1e0, *m=new double[1];
k=0; while (d>0.0) { k++; d=d-hf; } nk=k; 
for (k=0; k<nk; k++) { ras22[k]=ras22[k]*kor; srra2[k]=srra2[k]*kos; legr01[k]=legr01[k]*kos; }
d=0.0; for (k=0; k<nk; k++) d=d+ras22[k]*srra2[k]; //cout << "p = " << p << "\td = " << d << endl;
delete []ras01; delete []prgr01;
if (!vlpr) { delete []srra2; delete []legr01; return ras22; }
else if (vlpr==1) { delete []ras22; delete []legr01; return srra2; }
else if (vlpr==2) { delete []ras22; delete []srra2; delete []legr01; m[0]=d; return m; }
else if (vlpr==3) { delete []ras22; delete []srra2; return legr01; } 
else if (vlpr==4) { delete []ras22; delete []srra2; delete []legr01; m[0]=p; return m; } }
double *rasrasporporaz0(double *prgr0, double *legr0, double *raspr0, double rpn, int vyb, int n, double de)
{ struct raspporporazmm *rpr01, *rprp01, *roo=new struct raspporporazmm;
double s, e=1e-1, r, koef=1e6, *prgr00=new double[n], *legr00=new double[n], m, t, ht=1e0; 
int k, q, p, j, jk;
for (k=0; k<n; k++) prgr00[k]=prgr0[k]*koef;
legr00[0]=0.0; for (k=1; k<n; k++) legr00[k]=prgr00[k-1];
s=0.0; for (k=0; k<n; k++) if (prgr00[k]<(ht+e)) { p=k; s=s+raspr0[k]; } //размер пор до 1 мкм
roo->rprm=s; roo->slel=NULL; rpr01=roo; q=0; rprp01=roo; 
for (k=p; k<n; k++) { 
    r=prgr00[k]-legr00[k];
    if ((prgr00[k]>(ht+e)) && (r<(ht+e))) {
		rpr01=new raspporporazmm;
		rpr01->rprm=raspr0[k];
		rprp01->slel=rpr01;
		rpr01->slel=NULL;
		rprp01=rpr01;
		q++; }
    if (r>(1.0+e)) {
        m=raspr0[k]/r;
		t=r; j=0; while (t>e) { j++; t=t-ht; } jk=j;
        for (j=0; j<jk; j++) {
		rpr01=new raspporporazmm;
		rpr01->rprm=m; 
		rprp01->slel=rpr01; 
		rpr01->slel=NULL;
		rprp01=rpr01;
		q++; } } }
qg=q; double *ras01=new double[q], *prgrm=new double[q], *legrm=new double[q], prgr01=rpn; rpr01=roo;
k=0; legrm[0]=0.0; while ((rpr01) && (k<q)) 
{ prgrm[k]=prgr01;
if (k>0) legrm[k]=prgrm[k-1];
ras01[k]=rpr01->rprm;
rprp01=rpr01;
rpr01=rpr01->slel; 
delete rprp01;
k++; 
prgr01=prgr01+de; } //for (k=0; k<q; k++) cout << "rp " << k << " = " << ras01[k] << "\t";
delete []prgr00; delete []legr00;
if (!vyb) { delete []prgrm; delete []legrm; return ras01; }
else if (vyb==1) { delete []ras01; delete []legrm; return prgrm; }
else if (vyb==2) { delete []ras01; delete []prgrm; return legrm; } }