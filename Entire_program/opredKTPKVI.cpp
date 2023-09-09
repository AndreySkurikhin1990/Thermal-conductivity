#define _CRT_SECURE_NO_WARNINGS
#include <string>
#include <cstring>
#include <stdio.h>
#include <conio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
using namespace std;
//----
const int cemnk=11, cemdu=6;
double epsisred(double, double *, double *, int, double *, double *, int, int, int, double);
double *opredKTPTverKarkkvi(double *, double *, double, double, double, double, int, int, double *, double *, int, double *, int, double *, double *, double *, int, char *);
double **opredTempHolGor(double *, double *, int, double, int, double **, int, int, double *, int, double, char *);
double **napMasEKTPkviNac(double, double, double, double, int, int, char *, double, int, double, double, double *, int, double *, double *, int, double *, double *, int, double *, double *, double *, int, double);
double *poisMasKoefkvi(int, double *, int);
double *oprsodoxkvi(int, char *);
double **vydelPol(int, int, double **, double **, int, int, char *, int, int);
//-----
double *poisMasKoefkvi(int vyb, double *kktp, int n)
{
	int k=0, j=0, w=0, q=1; double t1=25e0, t2=5e2, dt=t2-t1, kn=0.0;
	for (k=0; k<n; k++) kktp[k]=0.0;
	if (vyb==3) { kktp[q]=0.00015; kktp[w]=0.068; } //350
	else if (vyb==4) { kktp[q]=0.000125; kktp[w]=0.082; //1
	kn=(0.156-0.087)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.156-kn; //2
	kn=(0.155-0.087)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.155-kn; //3 
	kktp[q]=0.00016; kktp[w] = 0.14; //4  //Spirina //выбор //kn=(0.146-0.087)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.146-kn; //к вопросу о стандартизации КВИ
	} //400
	else if (vyb==5) { kktp[q]=0.0001; kktp[w]=0.103; //из Ахтямова, 1991
	kn=(0.165-0.105)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.165-kn; //выбор //kn=(0.178-0.105)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.178-kn; //Двухслойные //kn=(0.152-0.105)/dt; kktp[q]=kn; kn=kn*t2; kktp[w] = 0.152-kn; //к вопросу о стандартизации КВИ
	} //500
	else if (vyb==6) { kktp[q]=0.00015; kktp[w]=0.116; //выбор
	kn=(0.201-0.12)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.201-kn;
	kn=(0.195-0.12)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.195-kn;
	kktp[q]=0.00015; kktp[w]=0.17; //Spirina
	kn=(0.196-0.12)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.196-kn; //Двухслойные //выбор //kn=(0.178-0.12)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.178-kn; //к вопросу о стандартизации КВИ
	} //600
	else if (vyb==7) { kktp[q]=0.00017; kktp[w]=0.146; 
	kn=(0.216-0.15)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.216-kn; 
	kn=(0.235-0.15)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.235-kn; //к вопросу о стандартизации КВИ
	kn=(0.251-0.16)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.251-kn; //Двухслойные //выбор
	} //700
	else if (vyb==8) { kktp[q]=0.00018;  kktp[w]=0.156; 
	kn=(0.226-0.16)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.226-kn; //выбор //kn=(0.25-0.16)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.25-kn; //kktp[q]=0.00014; kktp[w]=0.21; //Spirina //kn=(0.23-0.16)/dt; kktp[q]=kn; kn=kn*t2; kktp[w] = 0.23-kn; //к вопросу о стандартизации КВИ
	} //800
	else if (vyb == 9) { kktp[q]=0.00019;  kktp[w]=0.185; 
	kn=(0.23-0.195)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.23-kn; //выбор //kn=(0.29-0.195)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.29-kn;
	} //900
	else if (vyb==10) { kktp[q]=0.00025; kktp[w]=0.246; 
	kn=(0.287-0.25)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.287-kn; //выбор //kn=(0.36-0.25)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.36-kn; //kn=(0.35-0.25)/dt; kktp[q]=kn; kn=kn*t2; kktp[w]=0.35-kn; //к вопросу о стандартизации КВИ
	} //1000
	return kktp;
}
double **napMasEKTPkviNac(double wmg, double wsi, double wal, double porkvi, int vybves, int vybkvi, char *snmk, double ys0, int dmkok, double tnosck, double dtosck, double *etek, int cem, double *tkusck, double *kusck, int dmkoosck, double *dkosckt, double *dkosckm, int dkosckl, double *ktpvoz, double *temvoz, double *Prvoz, int dmkov, double dko)
{
	int vpk=vybkvi, k=0, j=0, f=cemdu, cemk=cem;
	double *kektpk=new double[dmkok], *kttkk=NULL, *stchsrkvi=NULL; if (!kektpk) { cout << snmk << endl; k=getchar(); exit(1); }
	double *sodoxkvi=NULL, **mu=NULL, **muv=NULL;
	double *ektpk=NULL, tek0=273.15, tkk=22.0+tek0; 
	kektpk=poisMasKoefkvi(vpk, kektpk, dmkok);
	//-----
	double hf=1e0, g=0.0, s=0.0, t=0.0, *po=NULL, nf = 0.0, tn=0.0;
	mu=new double*[f]; muv=new double*[f];
	double *qobk=NULL, *tholk=NULL, *tgork=NULL, *tsredk=NULL;
	if ((!ektpk) || (!mu) || (!muv)) { cout << snmk << endl; k=getchar(); exit(1); }
	for (k=0; k<f; k++) { po=new double[cemk]; if (!po) { cout << snmk << endl; k=getchar(); exit(1); } muv[k]=po; } 
	for (k=0; k<cemk; k++) { t=etek[k]; s=0.0; g=0.0; 
	for (j=0; j<dmkok; j++) { s=s+kektpk[j]*pow(t, g); g=g+hf; } ektpk[k]=s; } //for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpi[k] << "\t"; cout << endl;
	muv=opredTempHolGor(ektpk, etek, cemk, ys0, k, muv, f, cemk, etek, dmkok, tkk, snmk); //cem - длина массива efte
	k=0; j=1; mu=vydelPol(k, cemk, muv, mu, f, j, snmk, vybves, vpk);
	k=0; tholk=mu[k]; k++; tgork=mu[k]; k++; qobk=mu[k]; k++; 
	ektpk=mu[k]; k++; tsredk=mu[k]; etek=tsredk; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cemk=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t"; //for (k=0; k<cemk; k++) cout << "tem = " << etek[k] << "\tktp = " << ektpk[k] << "\tth = " << tholk[k] << "\ttg = " << tgork[k] << "\tqo = " << qobk[k] << "\t"; cout << endl;
	for (k=0; k<f; k++) { po=muv[k]; if (po) delete[]po; } 
	for (k=0; k<f; k++) { po=mu[k]; if (po) delete[]po; } po=NULL;
	if (muv) delete[]muv; if (mu) delete[]mu; if (kektpk) delete[]kektpk; 
	//-----
	stchsrkvi=new double[cemk]; if (!stchsrkvi) { cout << snmk << endl; k=getchar(); exit(1); }
	for (k=0; k<cemk; k++) { g=epsisred(etek[k], tkusck, kusck, dmkoosck, dkosckt, dkosckm, dkosckl, vybves, vybkvi, dko); 
	stchsrkvi[k]=g; } //for (k=0; k<cemk; k++) { cout << "tem = " << etek[k] << "\tst_ch = " << stchsrkvi[k] << "\t"; cout << endl; } 
	//-----
	sodoxkvi=oprsodoxkvi(vybkvi, snmk);
	if ((!dkosckm) || (!dkosckt) || (!sodoxkvi)) { cout << snmk << endl; k=getchar(); exit(1); }
	k=0; wal=sodoxkvi[k]; k++; wsi=sodoxkvi[k]; k++; wmg=sodoxkvi[k]; if (sodoxkvi) delete[]sodoxkvi; 
	kttkk=opredKTPTverKarkkvi(etek, ektpk, porkvi, wsi, wal, wmg, vybves, vybkvi, kusck, tkusck, dmkoosck, stchsrkvi, cemk, ktpvoz, temvoz, Prvoz, dmkov, snmk);
	int c=8, q=1;
	mu=new double*[c]; po=new double[q]; if ((!mu) || (!po)) { cout << snmk << endl; k=getchar(); exit(1); }
	k=0; po[k]=nf; mu[k]=kttkk; k++; mu[k]=etek; k++; mu[k]=po; k++; mu[k]=ektpk; k++;
	mu[k]=tgork; k++; mu[k]=tholk; k++; mu[k]=qobk; k++; mu[k]=stchsrkvi;
	return mu;
}
double *oprsodoxkvi(int vybkvi, char *snmk)
{ int k=0, j=6; double *sodoxkvi=new double[j];
if (!sodoxkvi) { cout << snmk << endl; k=getchar(); exit(1); }
double wal=0.0, wsi=0.0, wmg=0.0, salok=0.0, smgok=0.0, ssiok=0.0, ko=1e-2;
if (vybkvi == 4) {
	salok=33e0; smgok=15e0; ssiok=52e0; wal=25e0; wsi=11e0; wmg=4e1; } //КВИ-400
else if (vybkvi == 5) {
	salok=34e0; smgok=11e0; ssiok=54e0; wal=28e0; wmg=8e0; wsi=44e0; } //КВИ-500
else if (vybkvi == 6) {
	salok=36e0; smgok=9e0; ssiok=55e0; wal=3e1; wsi=7e0; wmg=45e0; } //КВИ-600
else if (vybkvi == 7) {
	salok=37e0; smgok=8e0; ssiok=55e0; wal=31e0; wmg=6e0; wsi=45e0; } //КВИ-700
else if (vybkvi == 8) {
	salok=38e0; smgok=7e0; ssiok=55e0; wal=3e1; wmg=5e0; wsi=45e0; } //КВИ-800
else if (vybkvi == 9) {
	salok=39e0; smgok=6e0; ssiok=55e0; wal=32e0; wmg=5e0; wsi=45e0; } //КВИ-900
else if (vybkvi == 10) {
	salok=39e0; smgok=6e0; ssiok=55e0; wal=32e0; wmg=4e0; wsi=45e0; } //КВИ-1000
else { cout << "Net takoy marki KVI!" << endl; k=getchar(); exit(1); }
	salok=salok*ko; smgok=smgok*ko; ssiok=ssiok*ko; wal=wal*ko; wsi=wsi*ko; wmg=wmg*ko;
	k=0; sodoxkvi[k]=wal;   k++;   sodoxkvi[k]=wsi; k++; sodoxkvi[k]=wmg;  k++; 
		 sodoxkvi[k]=salok; k++; sodoxkvi[k]=ssiok; k++; sodoxkvi[k]=smgok; 
return sodoxkvi;
}