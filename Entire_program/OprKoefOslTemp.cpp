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
double *koefoslab(double, double, double, double *, int, double *, char *);
double *reshMetKram(double **, double *, int, char *);
double **oslaintestepcher(int, char *);
double *RasKorOV(double, double, double, double, double, double, char *);
double UsredMasOV(double **, int, int);
//----
double *koefoslab(double wmg, double wsi, double wal, double *tere, int n, double *kuo, char *snm)
{ int lt=n, k=0; 
double wo=wmg+wsi+wal, kmg=0.0, kal=0.0, ksi=0.0, ht=1e0, eps=1e-6;
double *mgo=NULL, *alo=NULL, *sio=NULL, **mu=oslaintestepcher(lt, snm);
k=0; sio=mu[k]; k++; alo=mu[k]; k++; mgo=mu[k]; //for (k=0; k<n; k++) cout << "sio = " << sio[k] << "\talo = " << alo[k] << "\tmgo = " << mgo[k] << endl; //cout << "wmg = " << wmg << "\twsio = " << wsi << "\twal = " << wal << "\two = " << wo << endl;
for (k=0; k<n; k++) {
kmg=mgo[k]; kal=alo[k]; ksi=sio[k];
if (fabs(wo)>eps) 
kuo[k]=(kmg*wmg+kal*wal+ksi*wsi)/wo; 
else kuo[k]=0.0; } //for (k=0; k<n; k++) cout << "kuo = " << kuo[k] << "\ttere = " << tere[k] << endl;
if (mgo) delete []mgo; if (alo) delete []alo; if (sio) delete []sio; 
return kuo; }
double **oslaintestepcher(int dm, char *snm) //рассчитывает ослабление интегральной степени черноты для SiO2, Al2O3, MgO
{	int n=2, k=0, j=0, w=0, p=6, q=10, m=0, l=0, r=1, v=0;
	double **scv=new double*[p], **hsv=new double*[p], *sc=new double[dm], hf=1e0;
	double *hs=new double[n], **kor=new double*[q], t0=0.0, t1=t0;
	double *vscs=new double[dm], *vscm=new double[dm], *vsca=new double[dm], *po=NULL;
if ((!scv) || (!sc) || (!hsv) || (!hs) || (!kor)) { cout << snm << endl; k=getchar(); exit(1);}
	k=0; sc[k]=hf;    k++; sc[k]=0.94;  k++; sc[k]=0.87;  k++; sc[k]=0.801; k++; sc[k]=0.736;
	k++; sc[k]=0.676; k++; sc[k]=0.635; k++; sc[k]=0.590; k++; sc[k]=0.567; k++; sc[k]=0.543;
	k++; sc[k]=0.53;  k++; sc[k]=0.525; k++; sc[k]=0.515; k++; sc[k]=0.507; //степень черноты магнезита
	w=0; scv[w]=sc; sc=new double[dm];
	k=0; hs[k]=98e-2; k++; hs[k]=2e-2;
	k=0; t0=hs[k]; k++; t1=hs[k];
	k=0; hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1);
	hsv[w]=hs; w++; hs=new double[n]; //Магнезит: 98 % - MgO, 2 % - SiO2
	if ((!sc) || (!hs)) { cout << snm << endl; k=getchar(); exit(1); }
k=0; sc[k]=hf;    k++; sc[k]=0.976; k++; sc[k]=0.949; k++; sc[k]=0.905; k++; sc[k]=0.859;  
k++; sc[k]=0.812; k++; sc[k]=0.774; k++; sc[k]=0.737; k++; sc[k]=0.709; k++; sc[k]=0.681; 
k++; sc[k]=0.661; k++; sc[k]=0.639; k++; sc[k]=0.626; k++; sc[k]=0.620; //степень черноты шамота
scv[w]=sc; sc=new double[dm]; 
k=0; hs[k]=56e-2; k++; hs[k]=396e-3; 
k=0; t0=hs[k]; k++; t1=hs[k]; 
k=0; hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Шамот: 56 % - SiO2, 39,6 % - Al2O3 - шамот, считаем, что магния нет в шамоте
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << snm << endl; k=getchar(); exit(1);}
k=0; sc[k]=hf;     k++; sc[k]=98e-2;  k++; sc[k]=951e-3; k++; sc[k]=92e-2;  k++; sc[k]=883e-3;
k++; sc[k]=853e-3; k++; sc[k]=821e-3; k++; sc[k]=79e-2;  k++; sc[k]=767e-3; k++; sc[k]=746e-3;
k++; sc[k]=73e-2;  k++; sc[k]=715e-3; k++; sc[k]=705e-3; k++; sc[k]=692e-3; //степень черноты корундошамота
scv[w]=sc; sc=new double[dm]; 
k=0; hs[k]=28e-2; k++; hs[k]=7e-1; 
k=0; t0=hs[k]; k++; t1=hs[k]; 
k=0; hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Корундошамот: 28 % - SiO2, 70 % - Al2O3 - корундошамот
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << snm << endl; k=getchar(); exit(1);}
k=0; sc[k]=hf;     k++; sc[k]=983e-3; k++; sc[k]=936e-3; k++; sc[k]=867e-3; k++; sc[k]=819e-3; 
k++; sc[k]=721e-3; k++; sc[k]=659e-3; k++; sc[k]=593e-3; k++; sc[k]=541e-3; k++; sc[k]=49e-2;
k++; sc[k]=453e-3; k++; sc[k]=429e-3; k++; sc[k]=403e-3; k++; sc[k]=384e-3; //степень черноты каолинового теплоизоляционного кирпича (КТК)
scv[w]=sc; sc=new double[dm];
k=0; hs[k]=57e-2; k++; hs[k]=4e-1; 
k=0; t0=hs[k]; k++; t1=hs[k]; 
k=0; hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //КТК: 57 % - SiO2, 40 % - Al2O3
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << snm << endl; k=getchar(); exit(1); }
k=0; sc[k]=hf;     k++; sc[k]=984e-3; k++; sc[k]=941e-3; k++; sc[k]=882e-3; k++; sc[k]=813e-3;
k++; sc[k]=751e-3; k++; sc[k]=695e-3; k++; sc[k]=641e-3; k++; sc[k]=594e-3; k++; sc[k]=558e-3;
k++; sc[k]=53e-2;  k++; sc[k]=499e-3; k++; sc[k]=479e-3; k++; sc[k]=462e-3; //степень черноты муллита
scv[w]=sc; sc=new double[dm];
k=0; hs[k]=28e-2; k++; hs[k]=72e-2; 
k=0; t0=hs[k]; k++; t1=hs[k]; 
k=0; hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Муллит: 28 % - SiO2, 72 % - Al2O3
hsv[w]=hs; w++; hs=new double[n];
if ((!sc) || (!hs)) { cout << "No memory" << endl; k=getchar(); exit(1);}
k=0; sc[k]=1e0;    k++; sc[k]=984e-3; k++; sc[k]=953e-3; k++; sc[k]=917e-3; k++; sc[k]=854e-3;
k++; sc[k]=808e-3; k++; sc[k]=756e-3; k++; sc[k]=711e-3; k++; sc[k]=578e-3; k++; sc[k]=523e-3;
k++; sc[k]=495e-3; k++; sc[k]=468e-3; k++; sc[k]=448e-3; k++; sc[k]=429e-3; //степень черноты кремнезема
scv[w]=sc; 
k=0; hs[k]=985e-3; k++; hs[k]=1e-2; 
k=0; t0=hs[k]; k++; t1=hs[k]; 
k=0; hs[k]=t0/(t0+t1); k++; hs[k]=t1/(t0+t1); //Кремнезем: 98,5 % - SiO2, 1 % - Al2O3
hsv[w]=hs; w++; r=1; v=0;
for (j=0; j<dm; j++) { k=0; 
for (m=1; m<p; m++)
	for (l=m+1; l<p; l++) {
		kor[k]=RasKorOV(hsv[m][v], hsv[m][r], hsv[l][v], hsv[l][r], scv[m][j], scv[l][j], snm); k++; }
		vscs[j]=UsredMasOV(kor, q, v); vsca[j]=UsredMasOV(kor, q, r); } 
k=0; hs=hsv[k]; sc=scv[k]; 
for (k=0; k<dm; k++) vscm[k]=(sc[k]-vscs[k]*hs[r])/hs[v];
for (k=0; k<q; k++) { sc=kor[k]; if (sc) delete[]sc; }
for (k=0; k<p; k++) { hs=hsv[k]; if (hs) delete[]hs; } if (hsv) delete[]hsv; //0 - SiO2, 1 - Al2O3, 2 - MgO
for (k=0; k<p; k++) { po=scv[k]; if (po) delete[]po; } if (scv) delete[]scv;
double kmg=hf, kal=hf, ksi=hf;
for (k=0; k<dm; k++) {
kmg=vscm[k]; kal=vsca[k]; ksi=vscs[k];
if (kmg<0.0) kmg=0.0; if (kmg>hf) kmg=hf; 
if (kal<0.0) kal=0.0; if (kal>hf) kal=hf; 
if (ksi<0.0) ksi=0.0; if (ksi>hf) ksi=hf; 
vscm[k]=kmg; vsca[k]=kal; vscs[k]=ksi; }
k=3; double **mu=new double*[k];
k=0; mu[k]=vscs; k++; mu[k]=vsca; k++; mu[k]=vscm;
return mu; }
double *RasKorOV(double a11, double a12, double a21, double a22, double b1, double b2, char *snmov)
{ 
	int l=2, k=0, j=k+1, i=0, f=1;
	double **mas=new double*[l], *st=new double[l];
	double *bs=new double[l], *x=NULL, *kor=new double[l], n=0.0, e=1e0;
if ((!mas) || (!st) || (!bs)) { cout << snmov << endl; k=getchar(); exit(1);}
st[k]=a11; st[j]=a12; mas[k]=st; st=new double[l];
if (st) { st[k]=a21; st[j]=a22; mas[j]=st; } else { cout << snmov << endl; k=getchar(); exit(1);}
bs[i]=b1; bs[f]=b2;
x=reshMetKram(mas,bs,l,snmov);
if ((x[i]>=n) && (x[f]>=n) && (x[i]<=e) && (x[f]<=e)) {
	kor[i]=x[i]; kor[f]=x[f]; }
else { kor[i]=n; kor[f]=n; }
for (k=0; k<l; k++) { st=mas[k]; if (st) delete []st; } if (bs) delete []bs;
return kor; }
double UsredMasOV(double **kor, int q, int vy)
{ 
	double s1=0.0, s2=0.0, p=0.0, *s=NULL, hf=1e0; 
	int j=0, k=0, l=1;
for (j=0; j<q; j++) { s=kor[j];
    if ((s[k]>0.0) && (s[l]>0.0)) {
        s1=s1+s[k]; s2=s2+s[l]; p=p+hf; } }
if (fabs(p)>0.0) { s1=s1/p; s2=s2/p; }
if (!vy) return s1; 
else return s2; }