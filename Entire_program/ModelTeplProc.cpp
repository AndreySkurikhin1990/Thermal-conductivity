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
#include <time.h>
using namespace std;
//----------------------------------------
double **oprkoefKTPiskhchao(int, int, double *, double, int, double **, int, int, int, char *);
double **vydelPol(int, int, double **, double **, int, int, char *, int, int);
double **opredTempHolGor(double *, double *, int, double, int, double **, int, int, double *, int, double, char *);
double **napMasEKTPVer(int, int, int, double *, int, double, int, int, int, double **, int, int, char *);
double **arrTem_Netzsch(double **, char *); //массив температур - экспериментальные данные на Netzsch - хаотичная засыпка фракции 2-0,7 мм (исходный)
double *arrKTP_Netzsch(char *); //массив КТП вермикулита
double **PoiskZavVelTem(int, double **, int, int, double *, double, char *, int);
double *koefPrib(double *, double *, int, char *);
double opredKTPTKToch(double *, double *, double, int);
double *koefPribl(double *, double *, int, char *);
//----------------
double **PoiskZavVelTem(int v, double **mu, int rmu, int n, double *efte, double h, char *snmv, int dmkov)
{
	int k=1, f=rmu, j=0, q=j, w=j, b=10; double *poin=new double[k], cf=0.0, hf=1e0, nf=cf, tf=cf, e=1e-1, cfp=cf, cft=cf;
	double *temvht=NULL, *temvct=NULL, *tepvt=NULL, *po=NULL, *temhq=NULL, *ts=NULL, ***muu=NULL;
	double *temcq=NULL, *temvs=NULL, *tepv=NULL, *cemt=NULL, *ktpq=NULL, *ktp=NULL, **muv=NULL;
	int vyfrve=0, vysove=0, vymivmf=1, vyukve=1, c=0, nnvyfv=2, nnvysv=3, nnvmivmf=2;
	int nnvyuv=2, kvf=0, jvsv=0, qvmi=0, qvuv=0, d=0;
	double *nvyfv=new double[nnvyfv], *nvysv=new double[nnvysv], *nvmivmf=new double[nnvmivmf];
	double *nvyuv=new double[nnvyuv]; temhq=new double[n]; temcq=new double[n]; temvs=new double[n]; 
	tepv=new double[n]; cemt=new double[n]; ktpq=new double[n]; ktp=new double[n]; muu=new double**[b]; 
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv) || (!temhq) || (!temcq) || (!temvs) || (!tepv) || (!cemt) || (!ktpq) || (!ktp) || (!muu)) 
	{ cout << snmv << endl; k=getchar(); exit(1); } 
	for (k=0; k<n; k++) { temhq[k]=0.0; temcq[k]=0.0; temvs[k]=0.0; tepv[k]=0.0; cemt[k]=0.0; ktpq[k]=0.0; }
	k=0; nvyfv[k]=0; k++; nvyfv[k]=1; //фракции вермикулита
	k=0; nvysv[k]=0; k++; nvysv[k]=1; k++; nvysv[k]=2; //состояния вермикулита
	k=0; nvmivmf[k]=1; k++; nvmivmf[k]=2; //стационарные методы измерений - 2019 и 2020
	k=0; nvyuv[k]=1; k++; nvyuv[k]=2; k=0; j=-1; b=0; //укладка вермикулита
	for (kvf=0; kvf<nnvyfv; kvf++) {
		vyfrve=nvyfv[kvf];
		for (jvsv=0; jvsv<nnvysv; jvsv++) {
			vysove=nvysv[jvsv];
			if ((!vyfrve) || (vyfrve==2)) {
				for (qvmi=0; qvmi<nnvmivmf; qvmi++) {
					vymivmf=nvmivmf[qvmi]; 
					muv=new double*[f]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
					muv=napMasEKTPVer(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmkov, snmv); //d=0; temvct=muv[d]; d=5; po=muv[d]; k=0; nf=po[k]; //k=0; cf=e; while (cf<nf) { cf=cf+hf; k++; } q=k; //cout << "vmi = " << vymivmf << "\tvsv = " << jvsv << "\tvfv = " << kvf << endl; //for (k=0; k<q; k++) cout << "tc ( " << k << " ) = " << temvct[k] << "\t"; cout << endl; //
					muu[b]=muv; b++; 
					if (vysove) break; 
				} }
			else if (vyfrve==1) {
				for (qvuv=0; qvuv<nnvyuv; qvuv++) {
					vyukve=nvyuv[qvuv];
					muv=new double*[f]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
					muv=napMasEKTPVer(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmkov, snmv); // //d=0; temvct=muv[d]; d=5; po=muv[d]; k=0; nf=po[k]; //k=0; cf=e; while (cf<nf) { cf=cf+hf; k++; } q=k; //cout << "vyuv = " << vyukve << "\tvsv = " << jvsv << "\tvfv = " << kvf << endl; //for (k=0; k<q; k++) cout << "tc ( " << k << " ) = " << temvct[k] << "\t"; cout << endl; //
					muu[b]=muv; b++; } } } } 
	for (j=0; j<b; j++) { muv=muu[j];
	d=0; temvct=muv[d]; d++; temvht=muv[d]; d++; tepvt=muv[d]; d++; ktp=muv[d]; d++; ts=muv[d]; d++; po=muv[d]; k=0; nf=po[k]; 
					k=0; cf=e; while (cf<nf) { cf=cf+hf; k++; } q=k; //cout << "vmi = " << vymivmf << endl; //for (k=0; k<q; k++) cout << "tc ( " << k << " ) = " << temvct[k] << "\t"; cout << endl;
					d=1; k=0; cfp=tepvt[k]; for (k=d; k<q; k++) { cft=tepvt[k]; if ((cft<=cfp) && (d>0)) { d=-1; break; } cfp=cft; } //for (k=0; k<q; k++) cout << "c = " << c << "\tqo = " << tepvt[k] << "\tts = " << ts[k] << "\t"; 
						if (d>0) { for (w=0; w<q; w++) { 
							tf=ts[w]; cf=tepvt[w]*tf;
							for (k=0; k<n; k++) {
							if (fabs(tf-efte[k])<=hf) {
								temvs[k]=temvs[k]+tf;
								temhq[k]=temhq[k]+tf*temvht[k]; temcq[k]=temcq[k]+tf*temvct[k]; ktpq[k]=ktpq[k]+tf*ktp[k];
								tepv[k]=tepv[k]+cf; cemt[k]=cemt[k]+hf; //cout << "tf = " << tf << "\tptp = " << tepvt[k] << "\t";
								break; } } } } c++; }
	b=c; for (k=0; k<b; k++) { muv=muu[k]; for (j=0; j<f; j++) { po=muv[j]; delete[]po; } if (muv) delete[]muv; } //for (k=0; k<n; k++) cout << "cemt = " << cemt[k] << "\t";
		for (k=0; k<n; k++) {
			cf=temvs[k]; if (cf>e) {
				temhq[k]=temhq[k]/cf; temcq[k]=temcq[k]/cf; ktpq[k]=ktpq[k]/cf; tepv[k]=tepv[k]/cf; }
			cf=cemt[k]; if (cf>e) temvs[k]=temvs[k]/cf; else temvs[k]=0.0; }
	cf=0.0; for (k=0; k<n; k++) cf=cf+hf;
	k=0; poin[k]=cf; 
	k=0; mu[k]=temcq; k++; mu[k]=temhq; k++; mu[k]=tepv;
	k++; mu[k]=ktpq; k++; mu[k]=temvs; k++;  mu[k]=poin; 
	if (nvyfv) delete[]nvyfv; if (nvysv) delete[]nvysv; if (nvmivmf) delete[]nvmivmf; 
	if (nvyuv) delete[]nvyuv; if (muu) delete[]muu; 
	if (cemt) delete[]cemt;
	return mu;
}
double **opredTempHolGor(double *ektp, double *ete, int n, double h0, int v, double **mu, int rmu, int ni, double *efte, 
	int dmkoef, double tn, char *snmv) //моделирование процесса теплообмена в образце
{ //n - длина массива ektp, ni - длина efte, qob - плотность теплового потока, которую может создать лабораторная установка 
	int k=0, j=1, w=0, rm=k, nn=0, vybves=1, kem=3;
	double g=0.0, p=0.0, *qon=NULL, *tena=NULL, nit=1e10, hf=1e0, *po=new double[j], *kt=NULL;
	double *koeq=new double[dmkoef], ts=0.0, laef=0.0, e=1e-10;
	double tgor=0.0, thol=0.0, tsred=0.0, etem=0.0, qo=0.0, dt=hf, **muv=new double*[rmu], nf=0.0;
	double *qob=new double[ni], *temvs=new double[ni], *temvc=new double[ni], *temvh=new double[ni], *ktp=new double[ni];
	if ((!koeq) || (!po) || (!muv) || (!qob) || (!temvc) || (!temvh) || (!temvs)) { cout << snmv << endl; k=getchar(); exit(1); }
	g=0.0; for (k=0; k<ni; k++) { qob[k]=g; temvs[k]=g; temvc[k]=g; temvh[k]=g; ktp[k]=g; }
	muv=PoiskZavVelTem(k, muv, rmu, ni, efte, h0, snmv, dmkoef);
	k=2; qon=muv[k]; k=4; tena=muv[k]; k++; po=muv[k]; 
	k=0; nf=po[k]; w=k; g=e; while (g<nf) { g=g+hf; w++; } 
	for (k=0; k<dmkoef; k++) koeq[k]=0.0; 
	kt=koefPrib(qon, tena, ni, snmv); for (k=0; k<kem; k++) koeq[k]=kt[k]; if (kt) delete[]kt; 
		for (k=0; k<ni; k++) {
			ts=efte[k]; g=0.0; p=g;
			for (j=0; j<dmkoef; j++) { g=g+pow(ts, p)*koeq[j]; p=p+hf; } 
			qob[k]=g; //cout << "q ( " << k << " ) = " << g << "\tt = " << ts << endl;
		} if (koeq) delete[]koeq; for (k=0; k<rmu; k++) { po=muv[k]; if (po) delete[]po; } if (muv) delete[]muv;
		g=0.0; for (k=0; k<ni; k++) { qo=qob[k]; laef=0.0;
			if (qo>e) {
				laef=opredKTPTKToch(ektp, ete, efte[k], n); etem=efte[k]; 
				dt=qo*h0/laef;
                thol=etem-dt/2e0;
				tgor=thol+dt; }
			else { thol=0.0; tgor=0.0; laef=0.0; } 
			tsred=(tgor+thol)/2e0; 
			temvc[k]=thol; temvh[k]=tgor; temvs[k]=tsred; ktp[k]=laef; qob[k]=qo; 
		} g=0.0; for (k=0; k<ni; k++) g=g+hf; 
	k=1; po=new double[k]; muv=new double*[rmu]; if ((!po) || (!muv)) { cout << snmv << endl; k=getchar(); exit(1); } k=0; po[k]=g;
	muv[k]=temvc; k++; muv[k]=temvh; k++; muv[k]=qob; k++; muv[k]=ktp; k++; muv[k]=temvs; k++; muv[k]=po; 
	k=0; j=-1; mu=vydelPol(k, ni, muv, mu, rmu, j, snmv, vybves, k); 
	for (k=0; k<rmu; k++) { po=muv[k]; if (po) delete[]po; } if (muv) delete[]muv; 
	return mu;
}
double **oprkoefKTPiskhchao(int vmiv, int v, double *efte, double h, int n207, double **mu, int rmu, int cem, int dmk, char *snmv) //vmiv - выбор метода измерений
{
	int cemdum=2, nk=0, k=0, j=0, q=0, qn=0, nn=n207, cedumi=cemdum, w=rmu-1, f=rmu, kem=3; 
	double hf=1e0, nf207=0.0, *po=NULL, *koeft=NULL, *koefh=NULL, s=0.0, t=0.0;
	double *koefc=NULL, *koefs=NULL, *koefq=NULL, **muv=NULL, **muvv=NULL, r=0.0, e=1e-6;
	double *tsv=NULL, *tgv=NULL, *thv=NULL, *qov=NULL, *ktpv=NULL, *mt=NULL, *ts=NULL, *kt=NULL; 
	double te0=273.15, tkv=22.0+te0;
			if (!vmiv) {  //0 - установка Netzsch - нестационарный метод
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snmv << endl; k = getchar(); exit(1); }
		muv=arrTem_Netzsch(muv, snmv);
		k=0; po=muv[k]; k++; mt=muv[k];
		k=0; s=e; nf207=po[k]; while (s<nf207) { s=s+hf; k++; } nk=k; //cout << "nk = " << nk << "\t"; //nk=8 - длина массива ktpv
		ktpv=arrKTP_Netzsch(snmv); 
		if (muv) delete[]muv; if (po) delete[]po;
		muv=new double*[rmu]; k=0; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
		muv=opredTempHolGor(ktpv, mt, nk, h, k, muv, rmu, cem, efte, dmk, tkv, snmv); //cem - длина массива efte //for (k=0; k<f; k++) { po=muv[k]; delete[]po; }
		if (ktpv) delete[]ktpv; 
		if (mt) delete[]mt; 
		k=0; thv=muv[k];  k++; tgv=muv[k]; k++; qov=muv[k]; 
		k++; ktpv=muv[k]; k++; tsv=muv[k]; k++; po=muv[k]; 
		q=0; nf207=po[q]; s=e; while (s<nf207) { s=s+hf; q++; } qn=q; 
		koefh=new double[dmk]; koefc=new double[dmk]; koefq=new double[dmk]; koeft=new double[dmk]; koefs=new double[dmk]; cout << "dmk = " << "\t";
		if ((!koefh) || (!koefc) || (!koefq) || (!koeft) || (!koefs)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; }
		kt=koefPrib(thv, tsv, qn, snmv); for (k=0; k<kem; k++) koefc[k]=kt[k]; if (kt) delete[]kt; //kt=koefPribl(thv, tsv, qn, snmv); 
		kt=koefPrib(tgv, tsv, qn, snmv); for (k=0; k<kem; k++) koefh[k]=kt[k]; if (kt) delete[]kt; 
		kt=koefPrib(qov, tsv, qn, snmv); for (k=0; k<kem; k++) koefq[k]=kt[k]; if (kt) delete[]kt;  
		kt=koefPrib(ktpv, tsv, qn, snmv); for (k=0; k<kem; k++) koeft[k]=kt[k]; if (kt) delete[]kt; 
		kt=koefPrib(tsv, tsv, qn, snmv); for (k=0; k<kem; k++) koefs[k]=kt[k]; if (kt) delete[]kt; 
		for (k=0; k<rmu; k++) { po=muv[k]; if (po) delete[]po; } if (muv) delete[]muv; 
	}
	k=1; po=new double[k]; if (!po) { cout << snmv << endl; k=getchar(); exit(1); }
	s=0.0; for (k=0; k<rmu; k++) s=s+hf; k=0; po[k]=s; //cout << "rmu = " << s << "\t";
	k=0; mu[k]=koefc; k++; mu[k]=koefh; k++; mu[k]=koefq; k++; mu[k]=koeft; k++; mu[k]=koefs; k++; mu[k]=po;
	return mu; 
}
//-----
double **vydelPol(int v, int n, double **mu, double **muv, int f, int fl, char *snm, int vybves, int vybmar)
{ 
	int q=0, qn=n, k=0, m=0, w=f-1, j=1, nk=n, *mui=NULL, x=0, qk=0;
	double *po=NULL, nf=0.0, t=0.0, hf=1e0, *vm=NULL, e=1e-7, r=0.0, s=0.0;
	double te0=273.15, matepr=0.0, ko=1e2;
	if (!vybves) { if ((vybmar>=4) && (vybmar<=6)) matepr=11.5*ko;
	else if ((vybmar>=7) && (vybmar<=8)) matepr=13.0*ko;
	else if ((vybmar>=9) && (vybmar<=10)) matepr=12.7*ko;
	else if ((vybmar>=11) && (vybmar<=13)) matepr=13.0*ko; }
	if (vybves==1) matepr=11.0*ko; 
	if (vybves==2) matepr=11.0*ko; 
	if (vybves==3) { if (vybmar==4) matepr=1e1*ko; 
	else if (vybmar==5) matepr=10.5*ko; 
	else if ((vybmar>=6) && (vybmar<=7)) matepr=11.0*ko; 
	else if ((vybmar>=8) && (vybmar<=10)) matepr=11.5*ko; }
	matepr=matepr+te0;
	k=w; po=mu[k]; k=0; nf=po[k]; 
	t=e; while (t<nf) { k++; t=t+hf; } nk=k; mui=new int[nk]; if (!mui) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<nk; k++) {
		x=1; for (m=0; m<w; m++) { 
			po=mu[m]; if ((po[k]<e) && (x>0)) { x=-1; break; } } 
		if (x>0) mui[k]=k; else mui[k]=-1; }
		if (fl>0) { qn=1; for (m=0; m<=qn; m++) { po=mu[m]; for (k=0; k<nk; k++) if (po[k]>matepr) mui[k]=-1; } } 
		q=0; qn=-1; for (k=0; k<nk; k++) if (mui[k]>qn) q++; qn=q; 
	for (m=0; m<w; m++) {
	vm=new double[qn]; if (!vm) { cout << snm << endl; k=getchar(); exit(1); } 
	q=0; po=mu[m]; for (k=0; k<nk; k++) {
		x=mui[k]; if (x>=0) { vm[q]=po[x]; q++; } } muv[m]=vm; } qk=q; if (mui) delete[]mui; 
	nf=0.0; for (k=0; k<qn; k++) nf=nf+hf; 
	k=1; po=new double[k]; if (!po) { cout << snm << endl; k=getchar(); exit(1); } k=0; po[k]=nf; muv[w]=po; 
	return muv;
}