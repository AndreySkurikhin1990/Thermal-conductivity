#define _CRT_SECURE_NO_WARNINGS
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <windows.h>
#include <dos.h>
#include <iostream> 
#include <time.h>
using namespace std;
//-----
double *opredKTPTverKarkSha(double *, double *, double, double, double, double, int, double *, double *, double *, int, int, int, double *, double *, double *, int, char *);
double epsisred(double, double *, double *, int, double *, double *, int, int, int, double);
double **poisMasKoefSha(double *, int, int, int, int, int, char *, double);
double **napMasEKTPShaNac(double, double, double, int, int, int, int, char *, double, int, double, double, double, double, double *, int, double *, double *, int, double *, double *, int, double *, double *, double *, int, double);
double **vydelPol(int, int, double **, double **, int, int, char *, int, int);
double NovNapMas(int, int, int, int, int, int, int, int);
double **opredTempHolGor(double *, double *, int, double, int, double **, int, int, double *, int, double, char *);
//----------------------------------------
double **poisMasKoefSha(double *kektp, int n, int vybves, int vybsha, int vystsha, int vrsh, char *snm, double porsha)
{
	int k=0, j=0, m=3, w=4, q=3;
	double ko=1e-2, p20=2e1*ko, p24=24.0*ko, p30=3e1*ko, p33=33.0*ko, p16=16.0*ko, p10=1e1*ko, hf=1e0;
	double s28=28e-2, s38=38e-2, s45=45e-2, **mu=new double*[m], *soox=new double[w], *vsm=new double[q]; 
	if ((!mu) || (!soox) || (!vsm)) { cout << snm << endl; k=getchar(); exit(1); }
	double sfeosha=1.24e-2, smgosha=0.29e-2, salosha=41.9e-2, ssiosha=54e-2, e=1e-6, vs=e, vss=e; 
	if (!vrsh) { sfeosha=1.24*ko; smgosha=0.29*ko, salosha=41.9*ko, ssiosha=54.0*ko; porsha=21.8*ko; } //Роучка
	else if (vrsh==1) { sfeosha=1.21*ko; smgosha=0.29*ko; salosha=42.6*ko; ssiosha=hf-salosha-sfeosha-smgosha; porsha=11.0144*ko; } //ШПД-41
	else if (vrsh==2) { sfeosha=1.64*ko; smgosha=0.36*ko; salosha=35.9*ko; ssiosha=59.1e-2; porsha=25.2*ko; } //ШБ-1 2-1
	else if (vrsh==3) { sfeosha=1.66*ko; smgosha=0.4*ko; salosha=37.3*ko; ssiosha=57.4*ko; porsha=26.5*ko; } //ШВ-1 1-1
	else if (vrsh==4) { sfeosha=1.24*ko; smgosha=0.29*ko; salosha=41.9*ko; ssiosha=54*ko; porsha=11.5*ko; } //ШПД
	else if (vrsh==5) { sfeosha=1.54*ko; smgosha=0.3*ko; salosha=38.6*ko; ssiosha=56.5*ko; porsha=16.5*ko; } //ШКУ-32 3-1
	for (k=0; k<n; k++) kektp[k]=0.0;
	if ((salosha>=s28) && (salosha<=s38)) {
if ((porsha>=p20) && (porsha<p24)) { vybsha=0; vystsha=0; 
k=3; kektp[k]=-0.435e-9; k--; kektp[k]=0.685e-6; k--; kektp[k]=0.134e-3; k--; kektp[k]=0.725; }
else if ((porsha>=p24) && (porsha<p30)) { vybsha=0; vystsha=1; 
k=3; kektp[k]=-0.867e-9; k--; kektp[k]=1.77e-6; k--; kektp[k]=-0.523e-3; k--; kektp[k]=0.806; }  //задание коэффициентов - шамот средней пористости
else if ((porsha>=p16) && (porsha<p20)) { vybsha=1; 
k=3; kektp[k]=-0.397e-9; k--; kektp[k]=0.71e-6; k--; kektp[k]=0.011e-3; k--; kektp[k]=0.851; } //уплотненный шамот
else if ((porsha>=p30) && (porsha<=p33)) { vybsha=2; 
k=3; kektp[k]=-0.377e-9; k--; kektp[k]=0.918e-6; k--; kektp[k]=-0.338e-3; k--; kektp[k]=0.77; }  //низкоплотный шамот
else if ((porsha>=p10) && (porsha<p16)) { vybsha=3; 
k=3; kektp[k]=0.0; k--; kektp[k]=-0.607e-6; k--; kektp[k]=1.14e-3; k--; kektp[k] = 0.641; } } //повышенной плотности шамот
if ((salosha>s38) && (salosha<=s45)) {
if ((porsha>=p20) && (porsha<p24)) { vybsha=0; vystsha=0; 
k=3; kektp[k]=-0.124e-9; k--; kektp[k]=0.215e-6; k--; kektp[k]=0.125e-3; k--; kektp[k]=1.01; }
else if ((porsha>=p24) && (porsha<p30)) { vybsha=0; vystsha=1; 
k=3; kektp[k]=-0.333e-9; k--; kektp[k]=0.805e-6; k--; kektp[k]=-0.289e-3; k--; kektp[k]=0.903; } //задание коэффициентов - шамот средней пористости
else if ((porsha>=p16) && (porsha<p20)) { vybsha=1; 
k=3; kektp[k]=0.0; k--; kektp[k]=-0.154e-6; k--; kektp[k]=0.369e-3; k--; kektp[k]=1.03; } //уплотненный шамот
else if ((porsha>=p30) && (porsha<p33)) { vybsha=2; 
k=3; kektp[k]=-0.377e-9; k--; kektp[k]=0.918e-6; k--; kektp[k]=-0.338e-3; k--; kektp[k]=0.77; }  //низкоплотный шамот
else if ((porsha>=p10) && (porsha<p16)) { vybsha=3; 
k=3; kektp[k]=0.0; k--; kektp[k]=-0.141e-6; k--; kektp[k]=0.437e-3; k--; kektp[k]=1.32; }
} //повышенной плотности шамот
for (k=0; k<vybsha; k++) vs=vs+hf; for (k=0; k<vystsha; k++) vss=vss+hf; 
k=0; vsm[k]=vs; k++; vsm[k]=vss; k++; vsm[k]=porsha;
k=0; soox[k]=smgosha; k++; soox[k]=salosha; k++; soox[k]=ssiosha; k++; soox[k]=sfeosha;
k=0; mu[k]=kektp; k++; mu[k]=soox; k++; mu[k]=vsm;
return mu;
}
double **napMasEKTPShaNac(double wmg, double wsi, double wal, double porsha, int vybves, int vybsha, int vystsha, int vrsh, char *snms, double y0, int dmkos, double tnacs, double detes, double tnoscs, double dtoscs, double *etesha, int cem, double *tkuscs, double *kuscs, int dmkooscs, double *dkoscst, double *dkoscsm, int dkoscsl, double *ktpvozsha, double *temvozsha, double *Prvozsha, int dmvozsha, double dko)
{
	int cemdu=6, k=0, j=0, f=cemdu, cems=cem;
	double ys0=y0, *tgorsha=NULL, *tholsha=NULL, *qobsha=NULL, *vsm=NULL;
	double te0=273.15, tks=22.0+te0, *kektps=new double[dmkos], **mu=NULL;
	double *po=NULL, wfe=0.0, e=1e-6, vs=e, vss=e, t=e, hf=1e0;
	double g=0.0, s=0.0, nf=0.0, tn=0.0, **muv=new double*[f], *ektpsha=NULL;
	if (!kektps) { cout << snms << endl; k=getchar(); exit(1); }
	porsha=NovNapMas(vybves, vybsha, k, k, k, k, k, k); cout << porsha << "\t";
	mu=poisMasKoefSha(kektps, dmkos, vybves, vybsha, vystsha, vrsh, snms, porsha);
	k=0; kektps=mu[k]; k++; po=mu[k]; k++; 
	j=0; wmg=po[j]; j++; wal=po[j]; j++; wsi=po[j]; j++; wfe=po[j]; if (po) delete[]po;
	vsm=mu[k]; j=0; vs=vsm[j]; j++; vss=vsm[j]; j++; porsha=vsm[j]; cout << porsha << "\t";
	j=0; t=e; while (t<vs) { j++; t=t+hf; } vybsha=j; j=0; t=e; while (t<vs) { j++; t=t+hf; } vystsha=j;
	if (vsm) delete[]vsm; if (mu) delete[]mu; mu=new double*[f];
	//-----
	if ((!mu) || (!muv)) { cout << snms << endl; k=getchar(); exit(1); }
	for (k=0; k<f; k++) { po=new double[cems]; if (!po) { cout << snms << endl; k=getchar(); exit(1); } muv[k]=po; } 
	for (k=0; k<cems; k++) { t=etesha[k]-te0; s=0.0; g=0.0; 
	for (j=0; j<dmkos; j++) { s=s+kektps[j]*pow(t, g); g=g+hf; } ektpsha[k]=s; } //for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpi[k] << "\t"; cout << endl;
	muv=opredTempHolGor(ektpsha, etesha, cems, ys0, k, muv, f, cems, etesha, dmkos, tks, snms); //cem - длина массива efte
	k=0; j=1; mu=vydelPol(k, cems, muv, mu, f, j, snms, vybves, vybsha);
	k=0; tholsha=mu[k]; k++; tgorsha=mu[k]; k++; qobsha=mu[k]; k++; 
	ektpsha=mu[k]; k++; etesha=mu[k]; k++; po=mu[k]; 
	k=0; nf=po[k]; t=e; k=0; while (t<nf) { t=t+hf; k++; } cems=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cems; k++) { cout << "tem = " << etesha[k] << "\tktp = " << ektpsha[k] << "\tth = " << tholsha[k] << "\ttg = " << tgorsha[k] << "\tqo = " << qobsha[k] << "\t"; cout << endl; }
	for (k=0; k<f; k++) { po=muv[k]; if (po) delete[]po; } for (k=0; k<f; k++) { po=mu[k]; if (po) delete[]po; } po=NULL;
	if (muv) delete[]muv; if (mu) delete[]mu; if (kektps) delete[]kektps; 
	//-----
	double *stchsrsha=new double[cems]; if (!stchsrsha) { cout << snms << endl; k=getchar(); exit(1); }
	for (k=0; k<cems; k++) { g=epsisred(etesha[k], tkuscs, kuscs, dmkooscs, dkoscst, dkoscsm, dkoscsl, vybves, vybsha, dko);
		stchsrsha[k]=g; } 
	for (k=0; k<cems; k++) cout << "tem = " << etesha[k] << "\tst_ch = " << stchsrsha[k] << "\t"; cout << endl;
	double *kttks=opredKTPTverKarkSha(etesha, ektpsha, porsha, wsi, wal, wmg, dmkooscs, kuscs, tkuscs, stchsrsha, cems, vrsh, vystsha, ktpvozsha, temvozsha, Prvozsha, dmvozsha, snms);
	//-----
	int c=7, q=4; mu=new double*[c]; 
	vsm=new double[q]; if ((!mu) || (!vsm)) { cout << snms << endl; k=getchar(); exit(1); }
	k=0; vsm[k]=nf; k++; vsm[k]=vs; k++; vsm[k]=vss; k++; vsm[k]=porsha;
	k=0; mu[k]=kttks; k++; mu[k]=etesha; k++; mu[k]=vsm; k++; mu[k]=ektpsha; k++; 
	mu[k]=tgorsha; k++; mu[k]=tholsha; k++; mu[k]=qobsha;
	if (stchsrsha) delete[]stchsrsha;
	return mu;
}