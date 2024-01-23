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
//-----------
double *koefPrib(double *, double *, int, char *);
double *opredKTPTverKarkNach(double *, double *, double, double, double, double, int, int, int, int, double *, double *, int, 
	double *, int, int, double *, double *, double *, int, char *, int, int, int, int);
double epsisred(double, double *, double *, int, double *, double *, int, int, int, double);
double **oprkoefKTPiskhchao(int, int, double *, double, int, double **, int, int, int, char *);
double **napMasEKTPVer(int, int, int, double *, int, double, int, int, int, double **, int, int, char *);
double *danPoTemTepl840(double *, double *, int, char *);
double *danPoTemTepl841(double *, double *, int, char *);
double *danPoTemTepl842(double *, double *, int, char *);
double *danPoTemTepl843(double *, double *, int, char *);
double *danPoTemTepl844(double *, double *, int, char *);
double *danPoTemTepl845(double *, double *, int, char *);
double *danPoTemTepl2071(double *, double *, int, char *);
double *danPoTemTepl2072(double *, double *, int, char *);
double *danPoTemTepl2073(double *, double *, int, char *);
double *danPoTemTepl2074(double *, double *, int, char *);
double **arrTemCold84(double **, char *);
double *arrTemHigh84(char *);
double *arrTepPot84(char *);
double **arrTemCold207(double **, char *);
double *arrTemHigh207(char *);
double *arrTepPot207(char *);
double *danIskh207(double *, double *, int, int, int, char *);
double *opredKTPTverKarkNach(double *, double *, double, double, double, double, int, int, int, int, double *, double *, int, 
	double *, int, int, double *, double *, double *, int, char *, int, int, int, int);
double *arrKTP_2020(char *);
double *arrKTP_Netzsch(char *);
double **arrTem_Netzsch(double **, char *);
double *arrTem1_2020(char *);
double **arrTem2_2020(double **, char *);
double *arrTem3_2020(char *);
double **gotSredKTPTK(char *);
double **vydelPol(int, int, double **, double **, int, int, char *, int, int);
double **gotZnachExperVel(int, int, int, int, int, int);
//----------------
double **napMasEKTPVerNac(double wmg, double wsi, double wal, double porver, int vybves, int vyfv, int vyuv, int vmi, int vysv, 
	int vpmf, int vpkf, int isrp, char *snmv, int dmkov, double y0, double tnoscv, double dtoscv, double *etev, int cem, 
	double *tkuscv, double *kuscv, int dmkooscv, double *dkoscvt, double *dkoscvm, int dkoscvl, int fl, double *ktpvoz, 
	double *tevoz, double *Prvoz, int dmkvoz, double dko)
{
	int cemdu=6, k=1, nk=0, j=0, jk=0, f=cemdu, qn=0, cemv=cem, vybmar=0, z=0, rgz=0, vrsh=0, vystsha=0;
	double *stchsrver=NULL, *ektpv=NULL, *kttkv=NULL, *tsred=NULL, r=0.0;
	double hf=1e0, *po=NULL, s=0.0, ys0=y0, *tholv=NULL, *tgorv=NULL, *qobv=NULL;
	double **mu=new double*[f], **muv=new double*[f]; if ((!mu) || (!muv)) { cout << snmv << endl; k=getchar(); exit(1); } 
	double t=0.0, nf=t, g=t, e=1e-6, tn=t, wo=wmg+wsi+wal, wmgo=t, wsio=t, walo=t;
	if (fabs(wo)>e) { wmgo=wmg/wo; walo=wal/wo; wsio=wsi/wo; } //cout << "por = " << porver << "\twmg = " << wmgo << "\twal = " << walo << "\twsi = " << wsio << "\tce = " << cemv << "\t";
	muv=napMasEKTPVer(vyfv, vysv, vmi, etev, k, ys0, vyuv, cemv, k, muv, f, dmkov, snmv); //dmkov - длина массива коэффициентов //k=0; tholv=muv[k]; k++; tgorv=muv[k]; k++; qobv=muv[k]; k++; ektpv=muv[k]; k++; etev=muv[k];  k++; po=muv[k]; jk=cemv; k=0; nf=po[k]; t=e; while (t<nf) { t=t+hf; k++; } cemv=k; cout << "cemv = " << cemv << endl; for (k=0; k<cemv; k++) cout << "qob = " << qobv[k] << "\tektp = " << ektpv[k] << "\tte = " << etev[k] << "\ttc = " << tholv[k] << "\ttg = " << tgorv[k] << "\tts = " << etev[k] << endl; cout << endl << endl;
	k=0; j=1; mu=vydelPol(k, cemv, muv, mu, f, j, snmv, vybves, vybmar); 
	k=0; tholv=mu[k]; k++; tgorv=mu[k]; k++; qobv=mu[k];
	k++; ektpv=mu[k]; k++; etev=mu[k];  k++; po=mu[k]; 
	if (muv) delete[]muv; if (mu) delete[]mu; 
	jk=cemv; k=0; nf=po[k]; if (po) delete[]po; 
	t=e; while (t<nf) { t=t+hf; k++; } cemv=k; cout << "cemv = " << cemv << endl; 
	for (k=0; k<cemv; k++) cout << "qob = " << qobv[k] << "\tektp = " << ektpv[k] << "\tte = " << etev[k] << "\ttc = " << tholv[k] << "\ttg = " << tgorv[k] << "\tts = " << etev[k] << endl;
	if (!rgz) { stchsrver=new double[cemv]; if (!stchsrver) { cout << snmv << endl; k=getchar(); exit(1); } //степени черноты рассчитываются
	for (k=0; k<cemv; k++) { g=epsisred(etev[k], tkuscv, kuscv, dmkooscv, dkoscvt, dkoscvm, dkoscvl, vybves, vybmar, dko); 
	stchsrver[k]=g; } for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tst_ch = " << stchsrver[k] << "\t"; //cout << endl; //cout << "n = " << cemv << "\t";
	kttkv=opredKTPTverKarkNach(etev, ektpv, porver, wsio, walo, wmgo, vyfv, vysv, vpkf, isrp, kuscv, tkuscv, dmkooscv, stchsrver, 
		cemv, fl, ktpvoz, tevoz, Prvoz, dmkvoz, snmv, vybves, vrsh, vystsha, vybmar); }
	else { 
		mu=gotZnachExperVel(vyfv, vysv, vmi, vyuv, vybves, vybmar);
		if (tholv) delete[]tholv; if (tgorv) delete[]tgorv;
		k=0; tsred=mu[k];   k++; ektpv=mu[k];  k++; tgorv=mu[k]; k++; tholv=mu[k]; 
		k++; kttkv=mu[k]; k++;   stchsrver=mu[k]; k++; po=mu[k]; if (tsred) delete[]tsred; 
		k=0; t=e; r=po[k]; while (t<r) { k++; t=t+hf; } cem=k; }
	for (k=0; k<cemv; k++) cout << "tem = " << etev[k] << "\tkttkv = " << kttkv[k] << "\tektp = " << ektpv[k] << endl;
	if (z) { 
		FILE *fo=fopen("Stepen chernoty.txt", "at"); 
	if (!fo) { cout << "File is not open" << endl; k=getchar(); exit(1); }
	for (k=0; k<cemv; k++) fprintf(fo, "T = %0.2lf\tstch = %0.15lf\tvyfv = %d\tvysv = %d\tvyuv = %d\tvmimf = %d\tvybves = %d\n", 
		etev[k], stchsrver[k], vyfv, vysv, vyuv, vmi, vybves); fprintf(fo,"\n"); 
		fclose(fo); }
	int q=1, c=8; po=new double[q]; 
	mu=new double*[c]; 
	if ((!mu) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	k=0; po[k]=nf; 
	k=0; mu[k]=etev; k++; mu[k]=ektpv; k++; mu[k]=tgorv;     k++; mu[k]=tholv; 
	k++; mu[k]=qobv; k++; mu[k]=kttkv; k++; mu[k]=stchsrver; k++; mu[k]=po;
	return mu;
}
double **napMasEKTPVer(int vyfrve, int vysove, int vymeizvemafr, double *te, int v, double h, int vyukve, int nt, int wa, double **mu, 
	int rmu, int dmk, char *snmv)
{
	int cemdum=2, k=0, nvm=dmk, j=0, vy=0, n84=0, n207=0, cedu=rmu, q=0, cedumi=cemdum, nn=0, vyve=1, kem=3;
	double F=13.85*1e-4, nf84=0.0, nf207=0.0, *po=NULL, *tc84=NULL, s=0.0, t=0.0, e=1e-7; 
	double *th84=NULL, hf=1e0, *tp84=NULL, *th207=NULL, *tc207=NULL, *tp207=NULL;
	double *koeft=NULL, *koefh=NULL, *koefc=NULL, *koefq=NULL, *koefs=NULL, **muv=NULL, r=0.0; //cout << "vf = " << vyfrve << "\tvs = " << vysove << "\tvmi = " << vymeizvemafr << "\tvu = " << vyukve << endl; //for (k=0; k<nt; k++) cout << "te ( " << k << " ) = " << te[k] << "\t";
	double *ktp=NULL, *temvs=NULL, *temvh=NULL, *temvc=NULL, *tepv=NULL, *vm=NULL, *ma=NULL, *kt=NULL; 
	//-------
	if ((!vyfrve) || (vyfrve==2)) { //фракция 2-0,7 мм или фракция 1,6-0,35 мм
		th207=arrTemHigh207(snmv); 
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snmv << endl; k = getchar(); exit(1); }
		muv=arrTemCold207(muv, snmv); k=0; po=muv[k]; k++; tc207=muv[k]; 
		tp207=arrTepPot207(snmv);
		k=0; t=0.0; nf207=po[k]; while (t<nf207) { t=t+hf; k++; } n207=k; nn=k; //cout << "n207 = " << n207 << "\t";
		temvs=new double[n207]; temvh=new double[n207]; temvc=new double[n207]; 
		tepv=new double[n207]; ktp=new double[n207]; 
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snmv << endl; k=getchar(); exit(1); } //for (k=0; k<n207; k++) cout << "th = " << th207[k] << "\ttc = " << tc207[k] << "\tQ = " << tp207[k] << endl;
		for (k=0; k<n207; k++) {
		temvh[k]=th207[k]; temvc[k]=tc207[k]; temvs[k]=(temvh[k]+temvc[k])/2e0; 
		tepv[k]=tp207[k]/F; ktp[k]=fabs(temvh[k]-temvc[k])/h; ktp[k]=tepv[k]/ktp[k]; }
		if (tc207) delete[]tc207; if (th207) delete[]th207; if (tp207) delete[]tp207; if (po) delete[]po; 
				//исходный
				if (!vysove) {
					//установка Netzsch
					if (!vymeizvemafr) { 
		if (muv) delete[]muv; muv=new double*[cedu]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
		muv=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, n207, muv, cedu, nt, dmk, snmv); 
		k=0; koefc=muv[k]; k++; koefh=muv[k]; k++; koefq=muv[k]; 
		k++; koeft=muv[k]; k++; koefs=muv[k]; k++; po=muv[k]; }
		//данные 2020 года - ГОСТ 12170 - стационарный метод
					else if (vymeizvemafr==1) { 
		if (ktp) delete[]ktp; ktp=arrKTP_2020(snmv); 
		if (temvh) delete[]temvh; temvh=arrTem1_2020(snmv); 
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
		muv=arrTem2_2020(muv, snmv);
		k=0; po=muv[k]; nf207=po[k]; k++; if (temvc) delete[]temvc; temvc=muv[k]; 
		if (temvs) delete[]temvs; temvs=arrTem3_2020(snmv); 
		k=0; s=e; while (s<nf207) { s=s+hf; k++; } n207=k; k=0; q=kem; //cout << "nk 2020 = " << n207 << endl; for (k=0; k<n207; k++) cout << "ts ( " << k << " ) = " << temvs[k] << "\tktp = " << ktp[k] << endl;
		double *ts=new double[n207], *tst=new double[n207]; if ((!ts) || (!tst)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<n207; k++) { tst[k]=(temvh[k]+temvc[k])/2e0; ts[k]=tst[k]; } for (k=0; k<n207; k++) temvs[k]=ts[k];
		ts=danIskh207(ts, tst, n207, q, q, snmv); 
		ktp=danIskh207(ktp, tst, n207, q, q, snmv); 
		temvh=danIskh207(temvh, tst, n207, q, q, snmv); 
		temvc=danIskh207(temvc, tst, n207, q, q, snmv); 
		temvs=danIskh207(temvs, tst, n207, q, q, snmv); if (tst) delete[]tst; //for (k=0; k<q; k++) cout << "ts ( " << k << " ) = " << ts[k] << "\tktp = " << ktp[k] << "\ttg = " << temvh[k] << "\tth = " << temvc[k] << endl;
		koeft=new double[dmk]; koefh=new double[dmk];
		koefc=new double[dmk]; koefs=new double[dmk]; koefq=new double[dmk];
		if ((!koeft) || (!koefh) || (!koefc) || (!koefs) || (!koefq)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; } //cout << "\tktp\t";
		if (ktp) { kt=koefPrib(ktp, ts, q, snmv); for (k=0; k<kem; k++) koeft[k]=kt[k]; if (kt) delete[]kt; } //cout << "\ttemvh\t";
		if (temvh) { kt=koefPrib(temvh, ts, q, snmv); for (k=0; k<kem; k++) koefh[k]=kt[k]; if (kt) delete[]kt; } //cout << "\ttemvc\t";
		if (temvc) { kt=koefPrib(temvc, ts, q, snmv); for (k=0; k<kem; k++) koefc[k]=kt[k]; if (kt) delete[]kt; } //cout << "\ttemvs\t";
		if (temvs) { kt=koefPrib(temvs, ts, q, snmv); for (k=0; k<kem; k++) koefs[k]=kt[k]; if (kt) delete[]kt; } 
		tepv=new double[q]; if (!tepv) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<q; k++) {
			s=0.0; r=0.0; for (j=0; j<dmk; j++) { s=s+koeft[j]*pow(ts[k], r); r=r+hf; }
			r=temvh[k]-temvc[k]; tepv[k]=s*fabs(r)/h; }
			kt=koefPrib(tepv, ts, q, snmv); for (k=0; k<kem; k++) koefq[k]=kt[k]; 
			if (kt) delete[]kt; if (ts) delete[]ts; }
		//данные 2019 года - ГОСТ 12170
					else if (vymeizvemafr==2) { 
		koefq=danPoTemTepl2071(temvs, tepv, nvm, snmv); koefh=danPoTemTepl2071(temvs, temvh, nvm, snmv);
		koefc=danPoTemTepl2071(temvs, temvc, nvm, snmv); koeft=danPoTemTepl2071(temvs, ktp, nvm, snmv);
		koefs=danPoTemTepl2071(temvs, temvs, nvm, snmv); } }
				//после повторных измерений
				if (vysove==1) { 
			koefc=danPoTemTepl2072(temvs, temvc, nvm, snmv); koefh=danPoTemTepl2072(temvs, temvh, nvm, snmv); 
			koefq=danPoTemTepl2072(temvs, tepv, nvm, snmv); koeft=danPoTemTepl2072(temvs, ktp, nvm, snmv); 
			koefs=danPoTemTepl2072(temvs, temvs, nvm, snmv); }
				//после обжига при 1000 °С
				if (vysove == 2) { 
			koefc=danPoTemTepl2073(temvs, temvc, nvm, snmv); koefh=danPoTemTepl2073(temvs, temvh, nvm, snmv); 
			koefq=danPoTemTepl2073(temvs, tepv, nvm, snmv); koeft=danPoTemTepl2073(temvs, ktp, nvm, snmv); 
			koefs=danPoTemTepl2073(temvs, temvs, nvm, snmv); }
				if (vysove == 3) 
				//после повторного обжига при 1000 °С
				{ 
			koefc=danPoTemTepl2074(temvs, temvc, nvm, snmv); koefh=danPoTemTepl2074(temvs, temvh, nvm, snmv); 
			koefq=danPoTemTepl2074(temvs, tepv, nvm, snmv); koeft=danPoTemTepl2074(temvs, ktp, nvm, snmv); 
			koefs=danPoTemTepl2074(temvs, temvs, nvm, snmv); } }
	//фракция 8-4 мм
	else if (vyfrve==1) { 
		k=cedumi; muv=new double*[k]; muv=arrTemCold84(muv, snmv); 
		k=0; po=muv[k]; k++; tc84=muv[k]; 
		th84=arrTemHigh84(snmv); tp84=arrTepPot84(snmv);
		k=0; t=0.0; nf84=po[k]; while (t<nf84) { t=t+hf; k++; } n84=k; nn=k; //for (k=0; k<n84; k++) cout << "th = " << th84[k] << "\ttc = " << tc84[k] << "\ttp = " << tp84[k] << endl;
		temvs=new double[n84]; temvh=new double[n84]; temvc=new double[n84]; tepv=new double[n84]; ktp=new double[n84];
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<n84; k++) { 
		temvs[k]=(th84[k]+tc84[k])/2e0; 
		temvh[k]=th84[k]; temvc[k]=tc84[k]; tepv[k]=tp84[k]; 
		tepv[k]=tepv[k]/F; ktp[k]=fabs(th84[k]-tc84[k])/h; 
		ktp[k]=tepv[k]/ktp[k]; }
		if (po) delete[]po; if (th84) delete[]th84; if (tp84) delete[]tp84; if (tc84) delete[]tc84; 
		if (vyukve==1) { //плоско-параллельная засыпка
			if (!vysove) { //исходный
				koefc=danPoTemTepl840(temvs, temvc, nvm, snmv); koefh=danPoTemTepl840(temvs, temvh, nvm, snmv); 
				koefq=danPoTemTepl840(temvs, tepv, nvm, snmv); koeft=danPoTemTepl840(temvs, ktp, nvm, snmv); 
				koefs=danPoTemTepl840(temvs, temvs, nvm, snmv); }
			else if (vysove==1) { //после повторных измерений
				koefc=danPoTemTepl842(temvs, temvc, nvm, snmv); koefh=danPoTemTepl842(temvs, temvh, nvm, snmv); 
				koefq=danPoTemTepl842(temvs, tepv, nvm, snmv); koeft=danPoTemTepl842(temvs, ktp, nvm, snmv); 
				koefs=danPoTemTepl842(temvs, temvs, nvm, snmv); }
			else if (vysove==2) { //после обжига
				koefc=danPoTemTepl845(temvs, temvc, nvm, snmv); koefh=danPoTemTepl845(temvs, temvh, nvm, snmv); 
				koefq=danPoTemTepl845(temvs, tepv, nvm, snmv); koeft=danPoTemTepl845(temvs, ktp, nvm, snmv); 
				koefs=danPoTemTepl845(temvs, temvs, nvm, snmv); } }
		else if (vyukve==2) { //вертикальная засыпка
			if (!vysove) { //исходный
				koefc=danPoTemTepl841(temvs, temvc, nvm, snmv); koefh=danPoTemTepl841(temvs, temvh, nvm, snmv); 
				koefq=danPoTemTepl841(temvs, tepv, nvm, snmv); koeft=danPoTemTepl841(temvs, ktp, nvm, snmv); 
				koefs=danPoTemTepl841(temvs, temvs, nvm, snmv); }
			else if (vysove == 1) { //после повторных измерений
				koefc = danPoTemTepl844(temvs, temvc, nvm, snmv); koefh = danPoTemTepl844(temvs, temvh, nvm, snmv); 
				koefq = danPoTemTepl844(temvs, tepv, nvm, snmv); koeft = danPoTemTepl844(temvs, ktp, nvm, snmv); 
				koefs = danPoTemTepl844(temvs, temvs, nvm, snmv); }
			else if (vysove == 2) { //после обжига
				koefc = danPoTemTepl843(temvs, temvc, nvm, snmv); koefh = danPoTemTepl843(temvs, temvh, nvm, snmv); 
				koefq = danPoTemTepl843(temvs, tepv, nvm, snmv); koeft = danPoTemTepl843(temvs, ktp, nvm, snmv); 
				koefs = danPoTemTepl843(temvs, temvs, nvm, snmv); } } }
	if (temvs) delete[]temvs; if (temvh) delete[]temvh; 
	if (temvc) delete[]temvc; if (tepv) delete[]tepv; 
	if (ktp) delete[]ktp; if (muv) delete[]muv; 
	k=1; muv=new double*[cedu]; po=new double[k]; if ((!muv) || (!po)) { cout << snmv << endl; k=getchar(); exit(1); }
	t=0.0; for (k=0; k<cedu; k++) t=t+hf; k=0; po[k]=t; 
	k=0; muv[k]=koefc; k++; muv[k]=koefh; k++; muv[k]=koefq; 
	k++; muv[k]=koeft; k++; muv[k]=koefs; k++; muv[k]=po; //vy=cedu-1; for (k=0; k<vy; k++) { ma=muv[k]; for (j=0; j<dmk; j++) cout << "k = " << k << "\tkoef ( " << j << " ) = " << ma[j] << "\t"; cout << endl; }
	q=cedu-1; for (vy=0; vy<q; vy++) { 
		vm=muv[vy]; ma=new double[nt]; if (!ma) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<nt; k++) { t=te[k]; s=0.0; r=s;
	for (j=0; j<nvm; j++) { s=s+vm[j]*pow(t, r); r=r+hf; } 
	ma[k]=s; //if (!vy) cout << "te_s = " << t << "\tt_c = " << s << "\t"; if (vy==1) cout << "te_s = " << t << "\tt_h = " << s << "\t"; if (vy==2) cout << "te_s = " << t << "\tq = " << s << "\t"; if (vy==3) cout << "te_s = " << t << "\tktp = " << s << "\t";if (vy==4) cout << "te_s = " << t << "\tts = " << s << "\t";
	} //cout << endl;
	mu[vy]=ma; if (vm) delete[]vm; } po=muv[q]; if (po) delete[]po; if (muv) delete[]muv;
	k=1; po=new double[k]; if (!po) { cout << snmv << endl; k=getchar(); exit(1); }
	s=0.0; for (k=0; k<nt; k++) s=s+hf; k=0; po[k]=s; mu[q]=po; //cout << "nt = " << nt << "\t";
	muv=new double*[cedu]; if (!muv) { cout << snmv << endl; k=getchar(); exit(1); }
	k=0; j=-1; muv=vydelPol(k, nt, mu, muv, rmu, j, snmv, vyve, k); 
	for (k=0; k<cedu; k++) { po=muv[k]; mu[k]=po; } //k=5; po=mu[k]; k=0; t=po[k]; //cout << "s = " << s << "\tt = " << t << "\t"; 
	if (muv) delete[]muv; 
	return mu;
}
double *danIskh207(double *ma, double *x, int v, int n, int np, char *snmv)
{
		int k=0, p=0, q=0; 
		double s=0.0, t=0.0, hf=1e0, *m=new double[n], e=1e-10;
		if ((!m) || (!ma)) { cout << snmv << endl; k=getchar(); exit(1); }
		for (k=0; k<np; k++) { s=s+ma[k]*x[k]; t=t+x[k]; } if (t>e) s=s/t; else s=0.0;
		k=0; m[k]=s; k++;
		p=3; m[k]=ma[p]; k++;
		p=4; q=5; t=x[p]+x[q]; 
		if (t>e) s=(ma[p]*x[p]+ma[q]*x[q])/t; else s=0.0; m[k]=s; //600, 800 и 1000 °C
		if (ma) delete[]ma; 
		return m;
}
//---------
double **arrTem_Netzsch(double **mu, char *snmv)
{
	int k=1, le=8; 
	double *tem=new double[le], hf=1e0, nf=0.0, *po=new double[k], tev0=273.15;;
	if ((!tem) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	k=0; tem[k]=27.967;  k++; tem[k]=93.833;  k++; tem[k]=192.5; 
	k++; tem[k]=341.667; k++; tem[k]=491.467; k++; tem[k]=641.2; 
	k++; tem[k]=790.933; k++; tem[k]=991.133; k++;
	nf=0.0; for (k=0; k<le; k++) { tem[k]=tem[k]+tev0; nf=nf+hf; } k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=tem;
	return mu;
}
double *arrKTP_Netzsch(char *snmv)
{
	int k=0, le=8; double *ktp=new double[le], k1=1e-2, k2=1e-3;
	if (!ktp) { cout << snmv << endl; k = getchar(); exit(1); }
	ktp[k] = 7.0*k1;   k++; ktp[k] = 8.0*k1;   k++; ktp[k] = 103.0*k2; k++; 
	ktp[k] = 15.0*k1;  k++; ktp[k] = 207.0*k2; k++; ktp[k] = 283.0*k2; k++; 
	ktp[k] = 373.0*k2; k++; ktp[k] = 477.0*k2; 
	return ktp;
}
double *arrKTP_2020(char *snmv)
{
	int k=0, n=6; double *a = new double[n]; if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 0.175566644058715; k++; a[k] = 0.176801537812368; k++; a[k] = 0.179324717653617; k++;
	a[k] = 0.211768953068592; k++; a[k] = 0.237194543621728; k++; a[k] = 0.237231989760775;
	return a;
}
double *arrTem1_2020(char *snmv)
{
	int k=0, n=6; 
	double *a=new double[n], tev0=273.15; 
	if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 585.0; k++; a[k] = 600.0;  k++; a[k] = 585.0; k++; 
	a[k] = 800.0; k++; a[k] = 1000.0; k++; a[k] = 1000.0;
	for (k = 0; k < n; k++) a[k] = a[k] + tev0;
	return a;
}
double **arrTem2_2020(double **mu, char *snmv)
{
	int k=1, n=6; 
	double *a=new double[n], nf=0.0, hf=1e0, *po=new double[k], tev0=273.15;
	if ((!a) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) nf = nf + hf; k=0; po[k]=nf;
	a[k] = 119.75; k++; a[k] = 138.0; k++; a[k] = 129.5; k++; 
	a[k] = 200.0;  k++; a[k] = 273.0; k++; a[k] = 261.0;
	for (k=0; k<n; k++) a[k]=a[k]+tev0;
	k=0; mu[k]=po; k++; mu[k]=a;
	return mu;
}
double *arrTem3_2020(char *snmv)
{
	int k=0, n=6; 
	double *a=new double[n], tev0=273.15; 
	if (!a) { cout << snmv << endl; k = getchar(); exit(1); }
	a[k] = 377.0; k++; a[k] = 396.0; k++; a[k] = 383.5; k++; 
	a[k] = 548.0; k++; a[k] = 703.0; k++; a[k] = 697.25;
	for (k = 0; k < n; k++) a[k] = a[k] + tev0;
	return a;
}
//---------
double *danPoTemTepl840(double *temvs, double *temvh, int n, char *snmv) //Засыпка плоско-параллельная, исходный
{
	double tho1=0.0, tho2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], e=1e-10; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=14; q=16; tem1=temvs[p]+temvs[q];
	p=15; q=17; tem2=temvs[p]+temvs[q];
	p=14; q=16; if (tem1>e) tho1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1;
	p=15; q=17; if (tem2>e) tho2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0, t=0.0, s=0.0; 
	t=tem2-tem1; s=tho2-tho1; 
	if (fabs(t)>e) k1=s/t; else k1=0.0; 
	k2=tho2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl841(double *temvs, double *temvh, int n, char *snmv) //Засыпка вертикальная, исходный
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], e=1e-10; 
	int k=0, p=0, q=0, r=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=4; q=6; r=10; tem1=temvs[p]+temvs[q]+temvs[r];
	p=5; q=7; r=11; tem2=temvs[p]+temvs[q]+temvs[r];
	p=4; q=6; r=10; 
	if (fabs(tem1)>e) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem1;
	p=5; q=7; r=11; 
	if (fabs(tem2)>e) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem2;
	tem1=tem1/3e0; tem2=tem2/3e0;
	double k1=0.0, k2=0.0, s=0.0, t=0.0;
	s=tc2-tc1; t=tem2-tem1;
	if (fabs(t)>e) k1=s/t; else k1=0.0; 
	k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl842(double *temvs, double *temvh, int n, char *snmv) //Засыпка плоско-параллельная, повторы
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], e=1e-10; 
	int k=0, p=0, q=0, r=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=18; q=20; r=22; tem1=temvs[p]+temvs[q]+temvs[r];
	p=19; q=21; r=23; tem2=temvs[p]+temvs[q]+temvs[r];
	p=18; q=20; r=22; 
	if (fabs(tem1)>e) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem1; else tc1=0.0;
	p=19; q=21; r=23; 
	if (fabs(tem2)>e) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem2; else tc2=0.0;
	tem1=tem1/3e0; tem2=tem2/3e0;
	double k1=0.0, k2=0.0, s=0.0, t=0.0;
	t=tem2-tem1; s=tc2-tc1;
	if (fabs(t)>0.0) k1=s/t; else k1=0.0; 
	k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl843(double *temvs, double *temvh, int n, char *snmv) //Засыпка вертикальная, после обжига при 1000 °С
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], e=1e-10; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=0; q=2; tem1=temvs[p]+temvs[q];
	p=1; q=3; tem2=temvs[p]+temvs[q];
	p=0; q=2; if (fabs(tem1)>e) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1; else tc1=0.0;
	p=1; q=3; if (fabs(tem2)>e) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2; else tc2=0.0;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0, s=0.0, t=0.0; 
	s=tc2-tc1; t=tem2-tem1; 
	if (fabs(t)>e) k1=s/t; else k1=0.0; 
	k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl844(double *temvs, double *temvh, int n, char *snmv) //Засыпка вертикальная, повторы
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n];
	double k1=0.0, k2=0.0, e=1e-10, s=0.0, t=0.0; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=8; q=12; tem1=temvs[p]+temvs[q];
	p=9; q=13; tem2=temvs[p]+temvs[q];
	p=8; q=12; 
	if (fabs(tem1)>e) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1; else tc1=0.0;
	p=9; q=13; 
	if (fabs(tem2)>e) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2; else tc2=0.0;
	tem1=tem1/2e0; tem2=tem2/2e0;
	t=tem2-tem1; s=tc2-tc1;
	if (fabs(t)>e) k1=s/t; else k1=0.0; 
	k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl845(double *temvs, double *temvh, int n, char *snmv) //Засыпка плоско-параллельная, после обжига при 1000 °С
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n];
	double hf=1e0, s=0.0, t=0.0, e=1e-10;
	int k=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	k=24; tem1=temvs[k];
	k=25; tem2=temvs[k];
	k=24; if (fabs(tem1)>e) tc1=(temvh[k]*temvs[k])/tem1; else tc1=0.0;
	k=25; if (fabs(tem2)>e) tc2=(temvh[k]*temvs[k])/tem2; else tc2=0.0;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; 
	t=tem2-tem1; s=tc2-tc1;
	if (fabs(t)>0.0) k1=s/t; else k1=0.0; 
	k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
//---------
double *danPoTemTepl2071(double *temvs, double *tepv, int n, char *snmv) //Засыпка исходная, фракция 2-0,7 мм
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0;
	double *isdan=new double[n], hf=1e0, s=0.0, t=0.0, e=1e-10;
	int k=0;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k] = 0.0;
	k=0; tem1=temvs[k];
	k=1; tem2=temvs[k];
	k=0; if (fabs(tem1)>e) tepv1=(tepv[k]*temvs[k])/tem1; else tepv1=0.0;
	k=1; if (fabs(tem2)>e) tepv2=(tepv[k]*temvs[k])/tem2; else tepv2=0.0;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; 
	t=tem2-tem1; s=tepv2-tepv1;
	if (fabs(t)>0.0) k1=s/t; else k1=0.0; 
	k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1;
	return isdan;
}
double *danPoTemTepl2072(double *temvs, double *tepv, int n, char *snmv) //Фракция 2-0,7 мм (повторные измерения)
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0;
	double *isdan=new double[n], s=0.0, t=0.0; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=2; q=4; tem1=temvs[p]+temvs[q];
	p=3; q=5; tem2=temvs[p]+temvs[q];
	p=2; q=4; 
	if (fabs(tem1)>0.0) tepv1=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem1; else tepv1=0.0;
	p=3; q=5; 
	if (fabs(tem2)>0.0) tepv2=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem2; else tepv2=0.0;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0, e=1e-10; 
	t=tem2-tem1; s=tepv2-tepv1;
	if (fabs(t)>e) k1=s/t; else k1=0.0; 
	k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k] = k1; return isdan;
}
double *danPoTemTepl2073(double *temvs, double *tepv, int n, char *snmv) //Фракция 2-0,7 мм, после обжига при 1000 °С
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, e=1e-10;
	double *isdan=new double[n], t=0.0, s=0.0;
	int k=0, p=0, q=0;
	if (!isdan) { cout << snmv << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=6; q=8; tem1=temvs[p]+temvs[q];
	p=7; q=9; tem2=temvs[p]+temvs[q];
	p=6; q=8; 
	if (fabs(tem1)>e) tepv1=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem1; else tepv1=0.0;
	p=7; q=9; 
	if (fabs(tem2)>e) tepv2=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem2; else tepv2=0.0;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; 
	t=tem2-tem1; s=tepv2-tepv1;
	if (fabs(t)>0.0) k1=s/t; 
	else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl2074(double *temvs, double *tepv, int n, char *snmv) //Фракция 2-0,7 мм, после повторного обжига при 1000 °С
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, e=1e-10;
	double *isdan=new double[n], hf=1e0, s=0.0, t=0.0; 
	int k=0, p=0;
	if (!isdan) { cout << snmv << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=2; tem1=temvs[p];
	p=3; tem2=temvs[p];
	p=2; 
	if (fabs(tem1)>e) tepv1=(tepv[p]*temvs[p])/tem1; else tepv1=0.0;
	p=3; 
	if (fabs(tem2)>e) tepv2=(tepv[p]*temvs[p])/tem2; else tepv2=0.0;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; 
	t=tem2-tem1; s=tepv2-tepv1;
	if (fabs(t)>0.0) k1=s/t; else k1=0.0; 
	k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
//---------
double **arrTemCold84(double **mu, char *snmv)
{
	int k=1, n=26; 
	double *temcol=new double[n], *po=new double[k], nf=0.0, hf=1e0, tev0=273.15;
	if ((!temcol) || (!po)) { cout << snmv << endl; k = getchar(); exit(1); }
	k=0; 
	temcol[k] = 168.0;  k++; temcol[k] = 369.0;  k++; 
	temcol[k] = 168.0;  k++; temcol[k] = 369.0;  k++;
	temcol[k] = 148.6;  k++; temcol[k] = 356.5;  k++; 
	temcol[k] = 184.0;  k++; temcol[k] = 396.0;  k++;
	temcol[k] = 148.75; k++; temcol[k] = 350.0;  k++; 
	temcol[k] = 166.0;  k++; temcol[k] = 375.0;  k++;
	temcol[k] = 171.0;  k++; temcol[k] = 383.5;  k++; 
	temcol[k] = 106.75; k++; temcol[k] = 242.0;  k++;
	temcol[k] = 123.0;  k++; temcol[k] = 294.0;  k++; 
	temcol[k] = 111.0;  k++; temcol[k] = 240.0;  k++;
	temcol[k] = 109.0;  k++; temcol[k] = 232.75; k++; 
	temcol[k] = 127.0;  k++; temcol[k] = 291.0;  k++;
	temcol[k] = 221.0;  k++; temcol[k] = 443.75;
	for (k=0; k<n; k++) temcol[k]=temcol[k] + tev0;
	for (k=0; k<n; k++) nf=nf+hf; k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=temcol;
	return mu;
}
double *arrTemHigh84(char *snmv)
{
	int k=0, n=26; double *temh=new double[n], tev0=273.15;
	if (!temh) { cout << snmv << endl; k=getchar(); exit(1); }
	temh[k] = 545.0; k++; temh[k] = 927.0;  k++; 
	temh[k] = 530.0; k++; temh[k] = 925.0;  k++;
	temh[k] = 560.0; k++; temh[k] = 950.0;  k++; 
	temh[k] = 560.0; k++; temh[k] = 900.0;  k++;
	temh[k] = 558.0; k++; temh[k] = 950.0;  k++; 
	temh[k] = 540.0; k++; temh[k] = 920.0;  k++;
	temh[k] = 540.0; k++; temh[k] = 920.0;  k++; 
	temh[k] = 600.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 580.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 587.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 590.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 580.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 483.0; k++; temh[k] = 850.0;
	for (k=0; k<n; k++) temh[k]=temh[k]+tev0;
	return temh;
}
double *arrTepPot84(char *snmv)
{
	int k=0, n=26; double *tep=new double[n];
	if (!tep) { cout << snmv << endl; k=getchar(); exit(1); }
	tep[k] = 3.9805;  k++; tep[k] = 9.5532;   k++; 
	tep[k] = 3.6447;  k++; tep[k] = 9.4779;   k++;
	tep[k] = 3.55732; k++; tep[k] = 11.4997;  k++; 
	tep[k] = 4.4624;  k++; tep[k] = 11.6021;  k++;
	tep[k] = 4.7766;  k++; tep[k] = 11.5016;  k++; 
	tep[k] = 3.4023;  k++; tep[k] = 10.2068;  k++;
	tep[k] = 3.92812; k++; tep[k] = 11.17333; k++; 
	tep[k] = 3.5144;  k++; tep[k] = 8.6593;   k++;
	tep[k] = 2.977;   k++; tep[k] = 8.0448;   k++; 
	tep[k] = 3.352;   k++; tep[k] = 9.218;    k++;
	tep[k] = 3.0313;  k++; tep[k] = 7.7946;   k++; 
	tep[k] = 3.0671;  k++; tep[k] = 6.1342;   k++;
	tep[k] = 1.73466; k++; tep[k] = 4.32967;
	return tep;
}
double **arrTemCold207(double **mu, char *snmv)
{
	int k=1, n=10; 
	double *temcol=new double[n], *po=new double[k], nf=0.0, hf=1e0, tev0=273.15;
	if ((!temcol) || (!po)) { cout << snmv << endl; k=getchar(); exit(1); } 
	k=0;
	temcol[k] = 109.0;  k++; temcol[k] = 235.0; k++; 
	temcol[k] = 101.0;  k++; temcol[k] = 199.0; k++;
	temcol[k] = 108.75; k++; temcol[k] = 238.0; k++; 
	temcol[k] = 124.0;  k++; temcol[k] = 266.0; k++;
	temcol[k] = 111.0;  k++; temcol[k] = 262.0;
	for (k=0; k<n; k++) temcol[k]=temcol[k]+tev0;
	for (k=0; k<n; k++) nf=nf+hf; k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=temcol;
	return mu;
}
double *arrTemHigh207(char *snmv)
{
	int k=0, n=10; 
	double *temh=new double[n], tev0=273.15;
	if (!temh) { cout << snmv << endl; k=getchar(); exit(1); }
	temh[k] = 585.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 603.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 603.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 571.5; k++; temh[k] = 1000.75; k++;
	temh[k] = 583.0; k++; temh[k] = 1000.0;
	for (k=0; k<n; k++) temh[k]=temh[k]+tev0;
	return temh;
}
double *arrTepPot207(char *snmv)
{
	int k=0, n=10; 
	double *tepot=new double[n];
	if (!tepot) { cout << snmv << endl; k = getchar(); exit(1); }
	tepot[k] = 3.9596;  k++; tepot[k] = 8.6377; k++; 
	tepot[k] = 2.3003;  k++; tepot[k] = 5.3674; k++;
	tepot[k] = 3.56149; k++; tepot[k] = 7.123;  k++; 
	tepot[k] = 2.12992; k++; tepot[k] = 7.6956; k++;
	tepot[k] = 2.3003;  k++; tepot[k] = 6.9009;
	return tepot;
}
double **gotSredKTPTK(char *snmv)
{
int f=3, k=5, j=1; 
double **mu=new double*[f], *te=new double[k], *ktptk=new double[k], nf=0.0, hf=1e0, *po=new double[j], dt=1e2;
if ((!ktptk) || (!te) || (!mu) || (!po)) { cout << snmv << endl; k=getchar(); exit(1); }
j=0; ktptk[j]=1.0; j++; ktptk[j]=1.32; j++; ktptk[j]=1.66; j++; ktptk[j]=2.02; j++; ktptk[j]=2.39;
j=0; te[j]=573.15; for (j=1; j<k; j++) te[j]=te[j-1]+dt;
nf=0.0; for (j=0; j<k; j++) nf=nf+hf; j=0; po[j]=nf;
k=j; mu[k]=ktptk; k++; mu[k]=te; k++; mu[k]=po;
return mu;
}