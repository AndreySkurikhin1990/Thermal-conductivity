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
/*#define vmivmf 0 //выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный
#define vyfv 0 //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм
#define vyuv 1 //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
#define vysv 0 //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град С
using namespace std;
const int dso=12, dmko=4, isrp=0, vpkf=0, vpmf=0, cemdu=5, cemdum=2; //dmko - длина массива коэффициентов
const double ys0=3e1*1e-3, dete=1e2, te0=273.15, tnac=2e2+te0, templa=134.0*1e1+te0; 
const double por207=66.35*1e-2, poris84=55.75*1e-2, porin84=81.53*1e-2;
const double poro84=86.61*1e-2, poro16035=84.36*1e-2, por16035=83.97*1e-2;
const double pori440=(8e1+82e0)*1e-2/2e0, pori620=(75e0+78e0)*1e-2/2e0, pori860=(65e0+68e0)*1e-2/2, pori1000=(62e0+65e0)*1e-2/2e0;
const double porsha0=21.8*1e-2, porsha1=11.014*1e-2, porsha2=25.2*1e-2, porsha3=26.5*1e-2, porsha4=11.5*1e-2, porsha5=16.5*1e-2; //0 - Роучка, 1 - ШПД-41, 2 - ШБ-1 2-1, 3 - ШВ-1 1-1, 4 - ШПД, 5 - ШКУ-32 3-1
const double pork400=51.74*1e-2, pork500=52.07*1e-2, pork600=51.51*1e-2, pork700=39.75*1e-2;
const double pork800=40.85*1e-2, pork900=39.37*1e-2, pork1000=36.07*1e-2;
struct derevo {
	int otre; //1 - отражение или 0 - пропускание
	int ste; //номер стенки
	int vis; //видимость: 1 - виден, 0 - нет
	int lev; //номер уровня
	struct derevo *back; //указатель назад
	struct derevo *next; //указатель вперед
};
class KoePog { public: double alp, tem; KoePog *nex; };
//----------------------------------------
double *qob = NULL, *ete = NULL, tna=0.0, *tgor = NULL, *thol = NULL;
double por = 0.0, *tsred=NULL, *ektp = NULL, *kektp = NULL;
char *snm = NULL;
int cem=12; 
//----------------------------------------
double **oprkoefKTPiskhchao(int, int, double *, double *, double, double *, double *, double *, int, double **, int, int);
double opredKTPTKTochSha(double *, double *, double, int);
double *opredKTPTverKarkVerm(double *, double *, double, double, double, double, int, int, double, double, int, double *, double *, double *, int, int, int, int);
void initarrver(int, double, double, double);
void osvpamver();
double **opredtemphc(double *, double *, double *, double *, double *, int, int, double, double *, double *, double *, double **);
double **vydelPol(double *, double *, double *, double *, double *, int, int, double **, int);
double **opredTempHolGor(double *, double *, int, int, double, int, double **, int, double *, double *, int, double *, double *, double *);
double **napMasEKTP(int, int, int, double *, int, double, int, int, int, double **, int, int);
double **arrTem_Netzsch(double **); //массив температур - экспериментальные данные на Netzsch - хаотичная засыпка фракции 2-0,7 мм (исходный)
double *arrKTP_Netzsch(); //массив КТП вермикулита
double *arrKTP_2020();
double *arrTem1_2020();
double **arrTem2_2020(double **);
double *arrTem3_2020();
double *danPoTemTepl840(double *, double *, int);
double *danPoTemTepl841(double *, double *, int);
double *danPoTemTepl842(double *, double *, int);
double *danPoTemTepl843(double *, double *, int);
double *danPoTemTepl844(double *, double *, int);
double *danPoTemTepl845(double *, double *, int);
double *danPoTemTepl2071(double *, double *, int);
double *danPoTemTepl2072(double *, double *, int);
double *danPoTemTepl2073(double *, double *, int);
double *danPoTemTepl2074(double *, double *, int);
double **arrTemCold84(double **);
double *arrTemHigh84();
double *arrTepPot84();
double **arrTemCold207(double **);
double *arrTemHigh207();
double *arrTepPot207();
double *danIskh207(double *, int, int, int);
double **PoiskZavVelTem(int, double **, int, int);
void napMasEKTPVerNac();
double *koefPribSha(double *, double *, int, double *, char *);
double novNapMas(double, int, int, int, int, int, int, int);
//----------------
double novNapMas(double tnach, int vyve, int vymave, int vmimfv, int vyfrve, int vyukve, int vysove, int cemv)
{
	int k=0;
	snm=new char[dso]; if (!snm) { cout << "No memory"; k=getchar(); exit(1); } 
	k=0; snm[k]='N'; k++; snm[k]='o'; k++; snm[k]='_'; k++; snm[k]='m'; k++; snm[k]='e'; k++; 
	snm[k]='m';  k++; snm[k]='o'; k++; snm[k]='r'; k++; snm[k]='y'; k++; snm[k]='!'; k++; snm[k]='\0';
	if (!vyve) //если выбран шамот
	{
	if (!vymave) por=porsha0; else if (vymave==1) por=porsha1; else if (vymave==2) por=porsha2;
	else if (vymave==3) por=porsha3; else if (vymave==4) por=porsha4; else if (vymave==5) por=porsha5;
	}
	else if (vyve==1) { //если выбран вермикулит
	if ((!vysv) || (vysv==1)) {  //исходный или после повторных измерений
	if (!vyfv) { if (!vpmf) por = por207; else if (vpmf==1) por=por16035; } //выбор пористости мелкой фракции
	else if (vyfv==1) { if (!vpkf) por = poris84; else if (vpkf==1) por=porin84; } //выбор пористости крупной фракции
	else if (vyfv==2) por = poro16035; }
	if (vysv==2) { //после обжига
	if ((!vyfv) || (vyfv==2)) por=poro16035;
	else if (vyfv==1) por=poro84; } 
	}
	else if (vyve==2) { //если выбран ИТОМ
	if (!vymave) por=pori440; else if (vymave==1) por=pori620; 
	else if (vymave==2) por=pori860; else if (vymave==3) por=pori1000;
	}
	else if (vyve==3) { //если выбран КВИ
	if (vymave==4) por=pork400; else if (vymave==5) por=pork500; 
	else if (vymave==6) por=pork600; else if (vymave==7) por=pork700; 
	else if (vymave==8) por=pork800; else if (vymave==9) por=pork900;
	else if (vymave==10) por=pork1000;
	}
	if (cemv>0) ete=new double[cemv]; if (!ete) { cout << snm << endl; k = getchar(); exit(1); }
	k=0; ete[k]=tnach; 
	for (k=1; k<cemv; k++) ete[k]=ete[k-1]+dete; //cout << "por = " << porver << endl;
	return por;
}
double **PoiskZavVelTem(int v, double **mu, int rmu, int n)
{
	int k=1, f=cemdu; double *poin=new double[k], delT=dete, Tnacha=tnac, cf=0.0, hf=1e0, nf=0.0;
	double *temvht=NULL, *temvct=NULL, *tepvt=NULL, *po=NULL, *temhq=NULL, *temcq=NULL, *temvs=NULL, *tepv=NULL; 
	int vyfrve=0, vysove=0, vymivmf=1, vyukve=1, c=0, nnvyfv=2, nnvysv=3, nnvmivmf=2;
	int nnvyuv=2, kvf=0, jvsv=0, qvmi=0, qvuv=0, d=0;
	double *nvyfv=new double[nnvyfv], *nvysv=new double[nnvysv], *nvmivmf=new double[nnvmivmf];
	double *nvyuv=new double[nnvyuv], **muv=new double*[f]; temhq=new double[n]; temcq=new double[n];
	temvs=new double[n]; tepv=new double[n];
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv) || (!muv) || (!temhq) || (!temcq) || (!temvs) || (!tepv)) 
	{ cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) { temhq[k]=0.0; temcq[k]=0.0; temvs[k]=0.0; tepv[k]=0.0; }
	k=0; nvyfv[k]=0; k++; nvyfv[k]=1; //фракции вермикулита
	k=0; nvysv[k]=0; k++; nvysv[k]=1; k++; nvysv[k]=2; //состояния вермикулита
	k=0; nvmivmf[k]=1; k++; nvmivmf[k]=2; //стационарные методы измерений - 2019 и 2020
	k=0; nvyuv[k]=1; k++; nvyuv[k]=2; //укладка вермикулита
	for (kvf=0; kvf<nnvyfv; kvf++) {
		vyfrve=nvyfv[kvf];
		for (jvsv=0; jvsv<nnvysv; jvsv++) {
			vysove=nvysv[jvsv];
			if ((!vyfrve) || (vyfrve==2)) {
				for (qvmi=0; qvmi<nnvmivmf; qvmi++) {
					vymivmf=nvmivmf[qvmi];
					muv=napMasEKTP(vyfrve, vysove, vymivmf, temvs, k, ys0, vyukve, n, k, muv, f, dmko);
					d=0; temvct=muv[d]; d++; temvht=muv[d]; d++; tepvt=muv[d];
					for (k=0; k<n; k++) { cf=(temvct[k]+temvht[k])/2e0; 
					temvs[k]=temvs[k]+cf; cf=tepvt[k]; temhq[k]=temhq[k]+cf*temvht[k]; 
					temcq[k]=temcq[k]+cf*temvct[k]; tepv[k]=tepv[k]+cf; }
					for (k=0; k<f; k++) { po=muv[k]; delete[]po;}
					c++; } }
			else if (vyfrve==1) {
				for (qvuv=0; qvuv<nnvyuv; qvuv++) {
					vyukve=nvyuv[qvuv];
					muv=napMasEKTP(vyfrve, vysove, vymivmf, temvs, k, ys0, vyukve, n, k, muv, f, dmko);
					d=0; temvct=muv[d]; d++; temvht=muv[d]; d++; tepvt=muv[d];
					for (k=0; k<n; k++) { cf=(temvct[k]+temvht[k])/2e0; cout << "t_s = " << cf << "\t";
					temvs[k]=temvs[k]+cf; cf=tepvt[k]; temhq[k]=temhq[k]+cf*temvht[k]; cout << "q_ob = " << cf << "\n";
					temcq[k]=temcq[k]+cf*temvct[k]; tepv[k]=tepv[k]+cf; } 
					for (k=0; k<f; k++) { po=muv[k]; delete[]po;}
					c++; } } } }
	cout << "c = " << c << endl;
	cf=0.0; for (k=0; k<c; k++) cf=cf+hf; 
	if (cf>0.0) { for (k=0; k<n; k++) {  temhq[k]=temhq[k]/tepv[k]; temcq[k]=temcq[k]/tepv[k];
	tepv[k]=tepv[k]/cf; temvs[k]=temvs[k]/cf; cout << "t_s = " << temvs[k] << "\tq = " << tepv[k] << endl; } } 
	cf=0.0; for (k=0; k<n; k++) cf=cf+hf;
	k=0; poin[k]=cf; mu[k]=tepv; k++; mu[k]=temvs; k++; mu[k]=poin; k++; mu[k]=temhq; k++; mu[k]=temcq;
	if (nvyfv) delete[]nvyfv; if (nvysv) delete[]nvysv; 
	if (nvmivmf) delete[]nvmivmf; if (nvyuv) delete[]nvyuv; if (muv) delete[]muv;
	if (k<=rmu) return mu;
}
double **opredtemphc(double *effktp, double *efftem, double *tgv, double *thv, double *qon, int dlma, int n, double h, double *qob, double *tgor, double *thol, double **mu) //n=3 - длина массива коэффициентов приближающего многочлена, tgv - температура горячей стенки, thv - температура холодной стенки, dlma - длина массива ЭКТП
{
	int k=0, j=0; double *kq=new double[n], *kho=new double[n], *kgo=new double[n];
	double *tesr=new double[n], ts=0.0, g=0.0, p=0.0, hf=1e0, r=0.0, s=0.0, tego=0.0, teho=0.0;
	if ((!kq) || (!kho) || (!kgo) || (!tesr)) { cout << snm << endl; j=getchar(); exit(1); }
	for (k=0; k<n; k++) tesr[k]=(tgv[k]+thv[k])/2e0;
	kq=koefPribSha(qon, tesr, n, kq, snm);
	kgo=koefPribSha(tgv, tesr, n, kgo, snm);
	kho=koefPribSha(thv, tesr, n, kho, snm);
	for (k=0; k<dlma; k++) {
		ts=efftem[k];
		g=0.0; p=0.0; for (j=0; j<n; j++) { g=g+pow(ts, p)*kq[j]; p=p+hf; } 
		if (g<0.0) g=0.0; if (qob) qob[k]=g;
		p=effktp[k];
		g=fabs(g*h/p); 
		tego=ts+g/2e0; teho=ts-g/2e0;
		if ((tego<0.0) || (teho<0.0)) {
		if (tego<0.0) {
			r=0.0; s=0.0; for (j=0; j<n; j++) { r=r+pow(ts, s)*kgo[j]; s=s+hf; } 
			if (r<0.0) r=0.0; tego=r; }
		if (teho<0.0) {
			r=0.0; s=0.0; for (j=0; j<n; j++) { r=r+pow(ts, s)*kho[j]; s=s+hf; } 
			if (r<0.0) r=0.0; teho=r; } }
	if (tgor) tgor[k]=tego; if (thol) thol[k]=teho; }
	if (tesr) delete[]tesr; if (kq) delete[]kq; 
	if (kho) delete[]kho; if (kgo) delete[]kgo;
	k=0; mu[k]=tgor; k++; mu[k]=thol; k++; mu[k]=qob;
}
double **opredTempHolGor(double *ektp, double *ete, int n, int l, double h0, int v, double **mu, int rmu, double *temvs, double *qob, int ni, double *ktp, double *temvh, double *temvc) //моделирование процесса теплообмена в образце
{
	int cemf=cemdu, k=0, j=1, w=0, rm=k; double del=1e0, g, p, *qon, *tena, nit=1e6, hf=1e0, *po=new double[j], nf=0.0;
	double ep=1e-2, d=1e-3, Thna=ete[w], Thnac=0.0, *koeq=new double[l], ts, laef, dt=hf, **muv=new double*[cemf];
	double tgor=0.0, thol=0.0, tsred=0.0, etem=0.0;
	if ((!koeq) || (!po) || (!muv)) { cout << snm << endl; k = getchar(); exit(1); }
	muv=PoiskZavVelTem(k, muv, cemf, n);
	k=0; qon=muv[k]; k++; tena=muv[k]; k++; po=muv[k]; 
	k=0; nf=po[k]; w=k; g=0.0; while (g<nf) { g=g+hf; w++; }
	for (k=0; k<w; k++) cout << "q ( " << k << " ) = " << qon[k] << "\tt = " << tena[k] << endl;
	for (k=0; k<l; k++) koeq[k]=0.0;
	koeq=koefPribSha(qon, tena, w, koeq, snm);
	if (n<=ni) {
		for (k=0; k<n; k++) {
			ts=ete[k]; g=0.0; p=g;
			for (j=0; j<l; j++) { 
				g=g+pow(ts, p)*koeq[j]; p=p+hf; } 
			qob[k]=g; } //плотность теплового потока, которую может создать лабораторная установка //for (k=0; k<n; k++) cout << "qob = " << qob[k] << "\tektp = " << ektp[k] << "\tte = " << ete[k] << endl; 
		if (koeq) delete[]koeq; for (k=0; k<cemf; k++) { po=muv[k]; delete[]po; } delete[]muv;
		g=0.0; for (k=0; k<n; k++) {
			if (qob[k]>0.0) {
				laef=ektp[k]; p=0.0; Thnac=Thna+g*dt; del=hf; etem=ete[k];
				while ((del>ep) && (p<nit)) {
					thol=Thnac+p*d; //Tg - массив температур горячих стенок, Th - массив температур холодных стенок
					tgor=thol+qob[k]*h0/laef;
					del=(2e0*etem-(thol+tgor));
					p=p+hf; } 
				g=g+hf; }
			else { thol=0.0; tgor=0.0; } tsred=(tgor+thol)/2e0; 
			temvc[k]=thol; temvh[k]=tgor; temvs[k]=tsred;
		}
	} 
	k=0; mu[k]=temvc; k++; mu[k]=temvh; k++; mu[k]=qob; k++; mu[k]=temvs; k++; mu[k]=ktp; 
}
void napMasEKTPVerNac()
{
	int k=1, nk=0, j=0, jk=0, f=cemdu; double hf=1e0, *etevv, *isv, *mkovv, *p, s;
	por=novNapMas(tnac, k, j, vmivmf, vyfv, vyuv, vysv, cem);
	double **mu=new double*[f]; if (!mu) { cout << snm << endl; k=getchar(); exit(1); } 
	mu=napMasEKTP(vyfv, vysv, vmivmf, ete, k, ys0, vyuv, cem, k, mu, f, dmko);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; ektp=mu[k]; k++; tsred=mu[k];
	for (k=0; k<cem; k++) cout << "t_hol = " << thol[k] << "\tt_gor = " << tgor[k] << "\tq_ob = " << qob[k] << "\tts = " << tsred[k] << "\tktp = " << ektp[k] << endl;
	k=0; mu=vydelPol(thol, tgor, qob, ektp, tsred, k, cem, mu, f);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; tsred=mu[k]; k++; ektp=mu[k]; cout << "cem = " << cem << endl;
	double nf=0.0, tn=0.0, g=0.0;
	for (k=0; k<cem; k++) tsred[k]=(thol[k]+tgor[k])/2e0;
	etevv=new double[cem]; isv=new double[cem]; mkovv=new double[cem]; 
	if ((!etevv) || (!isv) || (!mkovv)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<cem; k++) etevv[k]=ete[k]; 
	jk=cem; k=0; nf=p[k]; delete[]p; s=0.0; while (s<nf) { s=s+hf; k++; } nk=k; cem=nk; k=0; tn=tsred[k]; tna=tn; 
	if (ete) delete[]ete; ete=new double[cem]; if (!ete) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<cem; k++) ete[k]=0.0;
	k=0; ete[k]=tn; for (k=1; k<cem; k++) ete[k]=ete[k-1]+dete; 
	if (etevv) delete[]etevv; 
	for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektp[k] << endl; cout << endl;
}
double **napMasEKTP(int vyfrve, int vysove, int vymeizvemafr, double *te, int v, double h, int vyukve, int nt, int wa, double **mu, int rmu, int dmk)
{
	int k=0, nvm=dmk, j=0, vy=0, n84=0, n207=0, cedu=cemdu;
	double F=13.85*1e-4, nf84=0.0, nf207=0.0, *po=NULL, *tc84=NULL; 
	double *th84=NULL, hf=1e0, *tp84=NULL, *th207=NULL, *tc207=NULL, *tp207=NULL;
	double *koeft=NULL, *koefh=NULL, *koefc=NULL, *koefq=NULL, *koefs=NULL, **muv=NULL, s, t, r, *ma=NULL;
	double *ktp=NULL, *temvs=NULL, *temvh=NULL, *temvc=NULL, *tepv=NULL, *vm=NULL;
	//-------
	if ((!vyfrve) || (vyfrve==2)) { //фракция 2-0,7 мм или фракция 1,6-0,35 мм
		th207=arrTemHigh207(); 
		k=cemdum; muv=new double*[k]; if (!muv) { cout << snm << endl; k = getchar(); exit(1); }
		muv=arrTemCold207(muv); k=0; po=muv[k]; k++; tc207=muv[k]; 
		tp207=arrTepPot207();
		k=0; t=0.0; nf207 = po[k]; while (t<nf207) { t=t+hf; n207++; } 
		temvs=new double[n207]; temvh=new double[n207]; 
		temvc=new double[n207]; tepv=new double[n207]; 
		ktp=new double[n207]; 
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snm << endl; k=getchar(); exit(1); } //for (k=0; k<n207; k++) cout << "th = " << th207[k] << "\ttc = " << tc207[k] << "\tQ = " << tp207[k] << endl;
		for (k=0; k<n207; k++) { temvs[k]=(th207[k]+tc207[k])/2e0; 
		temvh[k]=th207[k]; temvc[k]=tc207[k]; tepv[k]=tp207[k]; }
		for (k=0; k<n207; k++) { ktp[k]=fabs(temvh[k]-temvc[k])/h; tepv[k]=tepv[k]/F; ktp[k]=tepv[k]/ktp[k]; }
		if (th207) delete[]th207; if (tp207) delete[]tp207; if (tc207) delete[]tc207; if (po) delete[]po; 
			if (!vysove) { //исходный
				if (muv) delete[]muv;
				muv=new double*[cedu]; if (!muv) { cout << snm << endl; k=getchar(); exit(1); }
				muv=oprkoefKTPiskhchao(vymeizvemafr, 0, temvs, tepv, h, temvh, temvc, ktp, nvm, muv, rmu, nt); 
				k=0; koeft=muv[k]; k++; koefh=muv[k]; k++; koefc=muv[k]; k++; koefs=muv[k]; k++; koefq=muv[k]; }
		else if (vysove == 1) { //после повторных измерений
			koefc = danPoTemTepl2072(temvs, temvc, nvm); koefh = danPoTemTepl2072(temvs, temvh, nvm); 
			koefq = danPoTemTepl2072(temvs, tepv, nvm); koeft = danPoTemTepl2072(temvs, ktp, nvm); 
			koefs = danPoTemTepl2072(temvs, temvs, nvm); }
		else if (vysove == 2) { //после обжига при 1000 °С
			koefc = danPoTemTepl2073(temvs, temvc, nvm); koefh = danPoTemTepl2073(temvs, temvh, nvm); 
			koefq = danPoTemTepl2073(temvs, tepv, nvm); koeft = danPoTemTepl2073(temvs, ktp, nvm); 
			koefs = danPoTemTepl2073(temvs, temvs, nvm); }
		else if (vysove == 3) { //после повторного обжига при 1000 °С
			koefc = danPoTemTepl2074(temvs, temvc, nvm); koefh = danPoTemTepl2074(temvs, temvh, nvm); 
			koefq = danPoTemTepl2074(temvs, tepv, nvm); koeft = danPoTemTepl2074(temvs, ktp, nvm); 
			koefs = danPoTemTepl2074(temvs, temvs, nvm); } }
	else if (vyfrve == 1) { //фракция 8-4 мм
		k=cemdum; muv=new double*[k]; muv=arrTemCold84(muv); k=0; po=muv[k]; k++; tc84=muv[k]; 
		th84=arrTemHigh84();
		tp84=arrTepPot84();
		k=0; t=0.0; nf84 = po[k]; while (t<nf84) { t=t+hf; n84++; } //for (k=0; k<n84; k++) cout << "th = " << th84[k] << "\ttc = " << tc84[k] << "\ttp = " << tp84[k] << endl;
		temvs=new double[n84]; temvh=new double[n84]; 
		temvc=new double[n84]; tepv=new double[n84]; ktp=new double[n84];
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snm << endl; k = getchar(); exit(1); }
		for (k=0; k<n84; k++) { temvs[k]=(th84[k]+tc84[k])/2e0; temvh[k]=th84[k]; 
		temvc[k]=tc84[k]; tepv[k]=tp84[k]; tepv[k]=tepv[k]/F; ktp[k]=fabs(th84[k]-tc84[k])/h; ktp[k]=tepv[k]/ktp[k]; }
		if (po) delete[]po; if (th84) delete[]th84; if (tp84) delete[]tp84; if (tc84) delete[]tc84; 
		if (vyukve==1) { //плоско-параллельная засыпка
			if (!vysove) { //исходный
				koefc = danPoTemTepl840(temvs, temvc, nvm); koefh = danPoTemTepl840(temvs, temvh, nvm); 
				koefq = danPoTemTepl840(temvs, tepv, nvm); koeft = danPoTemTepl840(temvs, ktp, nvm); 
				koefs = danPoTemTepl840(temvs, temvs, nvm); }
			else if (vysove == 1) { //после повторных измерений
				koefc = danPoTemTepl842(temvs, temvc, nvm); koefh = danPoTemTepl842(temvs, temvh, nvm); 
				koefq = danPoTemTepl842(temvs, tepv, nvm); koeft = danPoTemTepl842(temvs, ktp, nvm); 
				koefs = danPoTemTepl842(temvs, temvs, nvm); }
			else if (vysove == 2) { //после обжига
				koefc = danPoTemTepl845(temvs, temvc, nvm); koefh = danPoTemTepl845(temvs, temvh, nvm); 
				koefq = danPoTemTepl845(temvs, tepv, nvm); koeft = danPoTemTepl845(temvs, ktp, nvm); 
				koefs = danPoTemTepl845(temvs, temvs, nvm); } }
		else if (vyukve == 2) { //вертикальная засыпка
			if (!vysove) { //исходный
				koefc = danPoTemTepl841(temvs, temvc, nvm); koefh = danPoTemTepl841(temvs, temvh, nvm); 
				koefq = danPoTemTepl841(temvs, tepv, nvm); koeft = danPoTemTepl841(temvs, ktp, nvm); 
				koefs = danPoTemTepl841(temvs, temvs, nvm); }
			else if (vysove == 1) { //после повторных измерений
				koefc = danPoTemTepl844(temvs, temvc, nvm); koefh = danPoTemTepl844(temvs, temvh, nvm); 
				koefq = danPoTemTepl844(temvs, tepv, nvm); koeft = danPoTemTepl844(temvs, ktp, nvm); 
				koefs = danPoTemTepl844(temvs, temvs, nvm); }
			else if (vysove == 2) { //после обжига
				koefc = danPoTemTepl843(temvs, temvc, nvm); koefh = danPoTemTepl843(temvs, temvh, nvm); 
				koefq = danPoTemTepl843(temvs, tepv, nvm); koeft = danPoTemTepl843(temvs, ktp, nvm); 
				koefs = danPoTemTepl844(temvs, temvs, nvm); } } }
	if (temvs) delete[]temvs; if (temvh) delete[]temvh; if (temvc) delete[]temvc; 
	if (tepv) delete[]tepv; if (ktp) delete[]ktp; if (muv) delete[]muv;
	muv=new double*[cedu];
	k=0; muv[k]=koefc; k++; muv[k]=koefh; k++; muv[k]=koefq; k++; muv[k]=koeft; k++; muv[k]=koefs;
	for (k=0; k<cedu; k++) { ma=muv[k]; for (j=0; j<dmk; j++) cout << "koef ( " << j << " ) = " << ma[j] << "\t"; cout << endl; }
	for (vy=0; vy<rmu; vy++) { vm=muv[vy]; ma=new double[nt]; if (!ma) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<nt; k++) { t=te[k]; s=0.0; r=s; ma[k]=r; 
	for (j=0; j<nvm; j++) { s=s+vm[j]*pow(t, r); r=r+hf; } ma[k]=s; cout << "te_s = " << t << "\twa = " << s << endl; } 
	po=muv[vy]; if (po) delete[]po; mu[vy]=ma; 
	} 
	if (muv) delete[]muv;
	return mu;
}
double **oprkoefKTPiskhchao(int vmiv, int v, double *temvs, double *tepv, double h, double *temvh, double *temvc, double *ktp, int n, double **mu, int rmu, int cem) //n=4 - длина массива коэффициентов, vmiv - выбор метода измерений
{
	int nk=n, k=0, j=0, q=0, qn=6, nn=n; 
	double hf=1e0, nf207=0.0, *po=NULL, *koeft=NULL, *koefh=NULL, s=0.0;
	double *koefc=NULL, *koefs=NULL, *koefq=NULL, **muv=NULL, r=0.0, e=1e-1;
	double *tsv=NULL, *tgv=NULL, *thv=NULL, *qov=NULL, *ktpv=NULL, *mt=NULL, *ts=NULL; //for (k=0; k<n207; k++) cout << "th = " << th207[k] << "\ttc = " << tc207[k] << "\tts = " << ts207[k] << endl;
	if (!vmiv) {  //0 - установка Netzsch - нестационарный метод
		k=cemdum; muv=new double*[k]; if (!muv) { cout << snm << endl; k = getchar(); exit(1); }
		muv=arrTem_Netzsch(muv);
		k=0; po=muv[k]; k++; mt=muv[k];
		k=0; s=0.0; nf207=po[k]; while (s<nf207) { s=s+hf; k++; } nk=k; //nk=8
		ktpv=arrKTP_Netzsch(); 
		if (muv) delete[]muv; if (po) delete[]po;
		k=cemdu; muv=new double*[k]; if (!muv) { cout << snm << endl; k = getchar(); exit(1); }
		muv=opredTempHolGor(ktpv, mt, nk, nn, h, 0, muv, cemdu, temvs, tepv, cem, ktp, temvh, temvc);
		k=0; thv=muv[k]; k++; tgv=muv[k]; k++; qov=muv[k]; k++; tsv=muv[k]; k++; ktp=muv[k];
		muv=vydelPol(thv, tgv, qov, ktpv, tsv, 0, nk, muv, cemdu);
		k=0; thv=muv[k]; k++; tgv=muv[k]; k++; qov=muv[k]; k++; tsv=muv[k]; k++; ktp=muv[k]; k++; po=muv[k];
		q=0; nf207=po[q]; s=0.0; while (s<nf207) { s=s+hf; q++; } qn = q; //cout << "qn = " << qn << endl;
		koefh=new double[nn]; koefc=new double[nn]; koefq=new double[nn]; koeft=new double[nn]; koefs=new double[nn];
		if ((!koefh) || (!koefc) || (!koefq) || (!koeft) || (!koefs)) { cout << snm << endl; k = getchar(); exit(1); }
		for (k=0; k<nn; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; }
		koeft=koefPribSha(ktpv, tsv, qn, koeft, snm); 
		koefh=koefPribSha(tgv, tsv, qn, koefh, snm); 
		koefc=koefPribSha(thv, tsv, qn, koefc, snm);
		koefs=koefPribSha(tsv, tsv, qn, koefs, snm); 
		koefq=koefPribSha(qov, tsv, qn, koefq, snm);
		if (tsv) delete[]tsv; if (tgv) delete[]tgv; if (thv) delete[]thv; if (ktpv) delete[]ktpv; 
		if (qov) delete[]qov; if (po) delete[]po; if (mt) delete[]mt; if (ktp) delete[]ktp;
	}
	else if (vmiv==1) { //данные 2020 года - ГОСТ 12170 - стационарный метод
		ktp=arrKTP_2020(); 
		tgv=arrTem1_2020(); 
		k=cemdum; muv=new double*[k]; if (!muv) { cout << snm << endl; k=getchar(); exit(1); }
		muv=arrTem2_2020(muv);
		k=0; po=muv[k]; nf207=po[k]; k++; thv=muv[k]; 
		tsv=arrTem3_2020();
		k=0; s=e; while (s<nf207) { s=s+hf; k++; } nk=k; k=0; q=3; //cout << "nk 2020 = " << nk << endl;
		ktp=danIskh207(ktp, k, nk, q); 
		tgv=danIskh207(tgv, k, nk, q); 
		thv=danIskh207(thv, k, nk, q); 
		tsv=danIskh207(tsv, k, nk, q);
		ts=new double[q]; if (!ts) { cout << snm << endl; k = getchar(); exit(1); } //nk - размер массива ts
		if ((ts) && (tgv) && (thv)) for (k=0; k<nk; k++) ts[k]=(tgv[k]+thv[k])/2e0;
		koeft=new double[n]; koefh=new double[n];
		koefc=new double[n]; koefs=new double[n]; koefq=new double[n];
		if ((!koeft) || (!koefh) || (!koefc) || (!koefs) || (!koefq)) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<n; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; } nk=q;
		if (ktp) koeft=koefPribSha(ktp, ts, nk, koeft, snm); 
		if (tgv) koefh=koefPribSha(tgv, ts, nk, koefh, snm);
		if (thv) koefc=koefPribSha(thv, ts, nk, koefc, snm);
		if (tsv) koefs=koefPribSha(tsv, ts, nk, koefs, snm);
		qov=new double[n]; if (!qov) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<nk; k++) {
			s=0.0; r=0.0; for (j=0; j<n; j++) { s=s+koeft[j]*pow(ts[j], r); r=r+hf; }
			qov[k]=s*fabs(tgv[k]-thv[k])/h; }
		koefq=koefPribSha(qov, ts, nk, koefq, snm); if (muv) delete[]muv; 
	}
	else if (vmiv==2) { //данные 2019 года - ГОСТ 12170
		koefq = danPoTemTepl2071(temvs, tepv, n);
		koefh = danPoTemTepl2071(temvs, temvh, n);
		koefc = danPoTemTepl2071(temvs, temvc, n);
		koeft = danPoTemTepl2071(temvs, ktp, n);
		koefs = danPoTemTepl2071(temvs, temvs, n);
	}
	k=0; mu[k]=koeft; k++; mu[k]=koefh; k++; mu[k]=koefc; k++; mu[k]=koefs; k++; mu[k]=koefq;
	return mu;
}
double *danIskh207(double *ma, int v, int n, int np)
{
	if (ma) {
		int k, p, q; double s=0.0, r=0.0, hf=1e0, *m=new double[n];
		if (!m) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<np; k++) { s=s+ma[k]; r=r+hf; } s=s/r;
		k=0; m[k]=s; k++; p=3; m[k]=ma[p]; k++; p=4; q=5; m[k]=(ma[p]+ma[q])/2e0; //600, 800 и 1000 °C
		delete[]ma; return m;
	}
	else return ma;
}
//---------
double **arrTem_Netzsch(double **mu)
{
	int k = 0, le = 8; double *tem = new double[le], hf = 1e0, nf = 0.0, *po = new double[1];
	if ((!tem) || (!po)) { cout << snm << endl; k = getchar(); exit(1); }
	tem[k]=27.967;  k++; tem[k]=93.833;  k++; tem[k]=192.5; k++; 
	tem[k]=341.667; k++; tem[k]=491.467; k++; tem[k]=641.2; k++; 
	tem[k]=790.933; k++; tem[k]=991.133; k++;
	nf=0.0; for (k=0; k<le; k++) { tem[k]=tem[k]+te0; nf=nf+hf; } k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=tem;
	return mu;
}
double *arrKTP_Netzsch()
{
	int k=0, le=8; double *ktp=new double[le]; 
	if (!ktp) { cout << snm << endl; k = getchar(); exit(1); }
	ktp[k] = 7.0*1e-2;   k++; ktp[k] = 8.0*1e-2;   k++; ktp[k] = 103.0*1e-3; k++; 
	ktp[k] = 15.0*1e-2;  k++; ktp[k] = 207.0*1e-3; k++; ktp[k] = 283.0*1e-3; k++; 
	ktp[k] = 373.0*1e-3; k++; ktp[k] = 477.0*1e-3; 
	return ktp;
}
double *arrKTP_2020()
{
	int k=0, n=6; double *a = new double[n]; if (!a) { cout << snm << endl; k = getchar(); exit(1); }
	a[k] = 0.175566644058715; k++; a[k] = 0.176801537812368; k++; a[k] = 0.179324717653617; k++;
	a[k] = 0.211768953068592; k++; a[k] = 0.237194543621728; k++; a[k] = 0.237231989760775;
	return a;
}
double *arrTem1_2020()
{
	int k=0, n=6; double *a=new double[n]; if (!a) { cout << snm << endl; k = getchar(); exit(1); }
	a[k] = 585.0; k++; a[k] = 600.0;  k++; a[k] = 585.0; k++; 
	a[k] = 800.0; k++; a[k] = 1000.0; k++; a[k] = 1000.0;
	for (k = 0; k < n; k++) a[k] = a[k] + te0;
	return a;
}
double **arrTem2_2020(double **mu)
{
	int k=1, n=6; double *a=new double[n], nf=0.0, hf=1e0, *po=new double[k];
	if ((!a) || (!po)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k = 0; k < n; k++) nf = nf + hf; k=0; po[k]=nf;
	a[k] = 119.75; k++; a[k] = 138.0; k++; a[k] = 129.5; k++; 
	a[k] = 200.0;  k++; a[k] = 273.0; k++; a[k] = 261.0;
	for (k=0; k<n; k++) a[k]=a[k]+te0;
	k=0; mu[k]=po; k++; mu[k]=a;
	return mu;
}
double *arrTem3_2020()
{
	int k=0, n=6; double *a=new double[n]; if (!a) { cout << snm << endl; k = getchar(); exit(1); }
	a[k] = 377.0; k++; a[k] = 396.0; k++; a[k] = 383.5; k++; 
	a[k] = 548.0; k++; a[k] = 703.0; k++; a[k] = 697.25;
	for (k = 0; k < n; k++) a[k] = a[k] + te0;
	return a;
}
//---------
double *danPoTemTepl840(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, исходный
{
	double tho1=0.0, tho2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=14; q=16; tem1=temvs[p]+temvs[q];
	p=15; q=17; tem2=temvs[p]+temvs[q];
	p=14; q=16; tho1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1;
	p=15; q=17; tho2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; 
	if (fabs(tem2-tem1)>0.0) k1=(tho2-tho1)/(tem2-tem1); else k1=0.0; 
	k2=tho2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl841(double *temvs, double *temvh, int n) //Засыпка вертикальная, исходный
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan = new double[n]; 
	int k=0, p=0, q=0, r=0;
	if (!isdan) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=4; q=6; r=10; tem1=temvs[p]+temvs[q]+temvs[r];
	p=5; q=7; r=11; tem2=temvs[p]+temvs[q]+temvs[r];
	p=4; q=6; r=10; if (fabs(tem1)>0.0) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem1;
	p=5; q=7; r=11; if (fabs(tem2)>0.0) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem2;
	tem1=tem1/3e0; tem2=tem2/3e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl842(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, повторы
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0, r=0;
	if (!isdan) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=18; q=20; r=22; tem1=temvs[p]+temvs[q]+temvs[r];
	p=19; q=21; r=23; tem2=temvs[p]+temvs[q]+temvs[r];
	p=18; q=20; r=22; if (fabs(tem1)>0.0) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem1; else tc1=0.0;
	p=19; q=21; r=23; if (fabs(tem2)>0.0) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q]+temvh[r]*temvs[r])/tem2; else tc2=0.0;
	tem1=tem1/3e0; tem2=tem2/3e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl843(double *temvs, double *temvh, int n) //Засыпка вертикальная, после обжига при 1000 °С
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; int k=0, p=0, q=0;
	if (!isdan) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=0; q=2; tem1=temvs[p]+temvs[q];
	p=1; q=3; tem2=temvs[p]+temvs[q];
	p=0; q=2; if (fabs(tem1)>0.0) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1;
	p=1; q=3; if (fabs(tem2)>0.0) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl844(double *temvs, double *temvh, int n) //Засыпка вертикальная, повторы
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], k1=0.0, k2=0.0; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=8; q=12; tem1=temvs[p]+temvs[q];
	p=9; q=13; tem2=temvs[p]+temvs[q];
	p=8; q=12; if (fabs(tem1)>0.0) tc1=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem1;
	p=9; q=13; if (fabs(tem2)>0.0) tc2=(temvh[p]*temvs[p]+temvh[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl845(double *temvs, double *temvh, int n) //Засыпка плоско-параллельная, после обжига при 1000 °С
{
	double tc1=0.0, tc2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], hf=1e0; int k=0;
	if (!isdan) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	k=24; tem1=temvs[k];
	k=25; tem2=temvs[k];
	k=24; if (fabs(tem1)>0.0) tc1=(temvh[k]*temvs[k])/tem1;
	k=25; if (fabs(tem2)>0.0) tc2=(temvh[k]*temvs[k])/tem2;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tc2-tc1)/(tem2-tem1); else k1=0.0; k2=tc2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
//---------
double *danPoTemTepl2071(double *temvs, double *tepv, int n) //Засыпка исходная, фракция 2-0,7 мм
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], hf=1e0; int k=0;
	if (!isdan) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k] = 0.0;
	k=0; tem1=temvs[k];
	k=1; tem2=temvs[k];
	k=0; if (fabs(tem1)>0.0) tepv1=(tepv[k]*temvs[k])/tem1;
	k=1; if (fabs(tem2)>0.0) tepv2=(tepv[k]*temvs[k])/tem2;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tepv2-tepv1)/(tem2-tem1); else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1;
	return isdan;
}
double *danPoTemTepl2072(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм (повторные измерения)
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=2; q=4; tem1=temvs[p]+temvs[q];
	p=3; q=5; tem2=temvs[p]+temvs[q];
	p=2; q=4; if (fabs(tem1)>0.0) tepv1=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem1;
	p=3; q=5; if (fabs(tem2)>0.0) tepv2=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tepv2-tepv1)/(tem2-tem1); 
	else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k] = k1; return isdan;
}
double *danPoTemTepl2073(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм, после обжига при 1000 °С
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n]; 
	int k=0, p=0, q=0;
	if (!isdan) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=6; q=8; tem1=temvs[p]+temvs[q];
	p=7; q=9; tem2=temvs[p]+temvs[q];
	p=6; q=8; tepv1=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem1;
	p=7; q=9; tepv2=(tepv[p]*temvs[p]+tepv[q]*temvs[q])/tem2;
	tem1=tem1/2e0; tem2=tem2/2e0;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tepv2-tepv1)/(tem2-tem1); 
	else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
double *danPoTemTepl2074(double *temvs, double *tepv, int n) //Фракция 2-0,7 мм, после повторного обжига при 1000 °С
{
	double tepv1=0.0, tepv2=0.0, tem1=0.0, tem2=0.0, *isdan=new double[n], hf=1e0; 
	int k=0, p=0;
	if (!isdan) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) isdan[k]=0.0;
	p=2; tem1=temvs[p];
	p=3; tem2=temvs[p];
	p=2; tepv1=(tepv[p]*temvs[p])/tem1;
	p=3; tepv2=(tepv[p]*temvs[p])/tem2;
	tem1=tem1/hf; tem2=tem2/hf;
	double k1=0.0, k2=0.0; if (fabs(tem2-tem1)>0.0) k1=(tepv2-tepv1)/(tem2-tem1); 
	else k1=0.0; k2=tepv2-k1*tem2;
	k=0; isdan[k]=k2; k++; isdan[k]=k1; 
	return isdan;
}
//---------
double **arrTemCold84(double **mu)
{
	int k=0, n=26; double *temcol=new double[n], *po=new double[1], nf=0.0, hf=1e0;
	if ((!temcol) || (!po)) { cout << snm << endl; k = getchar(); exit(1); }
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
	for (k=0; k<n; k++) temcol[k]=temcol[k] + te0;
	for (k=0; k<n; k++) nf=nf+hf; k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=temcol;
	return mu;
}
double *arrTemHigh84()
{
	int k=0, n=26; double *temh=new double[n];
	if (!temh) { cout << snm << endl; k=getchar(); exit(1); }
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
	for (k=0; k<n; k++) temh[k]=temh[k]+te0;
	return temh;
}
double *arrTepPot84()
{
	int k=0, n=26; double *tep=new double[n];
	if (!tep) { cout << snm << endl; k=getchar(); exit(1); }
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
double **arrTemCold207(double **mu)
{
	int k=1, n=10; double *temcol=new double[n], *po=new double[k], nf=0.0, hf=1e0;
	if ((!temcol) || (!po)) { cout << snm << endl; k=getchar(); exit(1); } k=0;
	temcol[k] = 109.0;  k++; temcol[k] = 235.0; k++; 
	temcol[k] = 101.0;  k++; temcol[k] = 199.0; k++;
	temcol[k] = 108.75; k++; temcol[k] = 238.0; k++; 
	temcol[k] = 124.0;  k++; temcol[k] = 266.0; k++;
	temcol[k] = 111.0;  k++; temcol[k] = 262.0;
	for (k=0; k<n; k++) temcol[k]=temcol[k]+te0;
	for (k=0; k<n; k++) nf=nf+hf; k=0; po[k]=nf;
	k=0; mu[k]=po; k++; mu[k]=temcol;
	return mu;
}
double *arrTemHigh207()
{
	int k=0, n=10; double *temh=new double[n];
	if (!temh) { cout << snm << endl; k=getchar(); exit(1); }
	temh[k] = 585.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 603.0; k++; temh[k] = 1000.0; k++;
	temh[k] = 603.0; k++; temh[k] = 1000.0; k++; 
	temh[k] = 571.5; k++; temh[k] = 1000.75; k++;
	temh[k] = 583.0; k++; temh[k] = 1000.0;
	for (k=0; k<n; k++) temh[k]=temh[k]+te0;
	return temh;
}
double *arrTepPot207()
{
	int k=0, n=10; double *tepot=new double[n];
	if (!tepot) { cout << snm << endl; k = getchar(); exit(1); }
	tepot[k] = 3.9596;  k++; tepot[k] = 8.6377; k++; 
	tepot[k] = 2.3003;  k++; tepot[k] = 5.3674; k++;
	tepot[k] = 3.56149; k++; tepot[k] = 7.123;  k++; 
	tepot[k] = 2.12992; k++; tepot[k] = 7.6956; k++;
	tepot[k] = 2.3003;  k++; tepot[k] = 6.9009;
	return tepot;
}
double **vydelPol(double *temvcs, double *temvhs, double *qos, double *ektpvs, double *temvss, int v, int n, double **mu, int f)
{ int q=0, qn=0, k=0, m=0, w=f-1, j=1;
	q=0; for (k=0; k<n; k++) 
		if ((temvcs[k]>0.0) && (temvhs[k]>0.0) && (qos[k]>0.0) && (ektpvs[k]>0.0) && (temvss[k]>0.0) && (temvhs[k]<templa)) 
			q++; qn=q; cem=qn;
	double hf=1e0, nf=0.0, *p=new double[j], *vm=NULL; if (!p) { cout << snm << endl; k=getchar(); exit(1); } 
	for (k=0; k<qn; k++) nf=nf+hf; k=0; p[k]=nf;
	for (m=0; m<w; m++) {
	vm=new double[qn]; if (!vm) { cout << snm << endl; k=getchar(); exit(1); } 
	q=0; for (k=0; k<n; k++) {
		if ((temvcs[k]>0.0) && (temvhs[k]>0.0) && (qos[k]>0.0) && (ektpvs[k]>0.0) && (temvss[k]>0.0) && (temvhs[k]<templa)) {
			if (!m) vm[q]=temvcs[k]; else if (m==1) vm[q]=temvhs[k]; 
			else if (m==2) vm[q]=qos[k]; else if (m==3) vm[q]=temvss[k]; 
			else if (m==4) vm[q]=ektpvs[k]; q++; } } mu[m]=vm; } 
	if (temvcs) delete[]temvcs; if (temvhs) delete[]temvhs; if (qos) delete[]qos; if (ektpvs) delete[]ektpvs; if (temvss) delete[]temvss;
	k=0; temvcs=mu[k]; k++; temvhs=mu[k]; k++; qos=mu[k]; k++; temvss=mu[k]; k++; ektpvs=mu[k];
	return mu; 
}
double *koefPribSha(double *ktp, double *te, int le, double *ko, char *snome)
{
	int k, kem=3; double **A=new double*[kem], *AA, *b=new double[kem];
	double yx2=0.0, yx=0.0, p=0.0, hf=1e0, x4=0.0, x3=0.0, x2=0.0;
	double x=0.0, y=0.0, de=0.0, de1=0.0, de2=0.0, de3=0.0;
	if (A) { for (k=0; k<kem; k++) { AA=new double[kem]; 
	if (AA) A[k]=AA; else { cout << snome << endl; k=getchar(); exit(1); } } }
	else { cout << snome << endl; k=getchar(); exit(1); }
	if ((!ko) || (!b)) { cout << snome << endl; k=getchar(); exit(1); }
	for (k=0; k<le; k++) {
		yx2=yx2+ktp[k]*pow(te[k], 2e0); 
		yx=yx+ktp[k]*te[k]; 
		y=y+ktp[k];
		x4=x4+pow(te[k], 4e0); 
		x3=x3+pow(te[k], 3e0); 
		x2=x2+pow(te[k], 2e0); 
		x=x+te[k];
	} //применение метода наименьших квадратов
	k=0; b[k]=yx2; k++; b[k]=yx; k++; b[k]=y; 
	p = 0.0; for (k = 0; k < le; k++) p = p + hf;
	A[0][0] = x4; A[0][1] = x3; A[0][2] = x2; A[1][0] = x3; A[1][1] = x2; A[1][2] = x; A[2][0] = x2; A[2][1] = x; A[2][2] = p;
	de = A[0][0] * (A[2][2] * A[1][1] - A[2][1] * A[1][2]) - A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	de1 = b[0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) - A[0][1] * (b[1] * A[2][2] - b[2] * A[1][2]) + A[0][2] * (b[1] * A[2][1] - b[2] * A[1][1]);
	de2 = A[0][0] * (b[1] * A[2][2] - b[2] * A[1][2]) - b[0] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * b[2] - A[2][0] * b[1]);
	de3 = A[0][0] * (A[1][1] * b[2] - A[2][1] * b[1]) - A[0][1] * (A[1][0] * b[2] - A[2][0] * b[1]) + b[0] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	k=2; ko[k] = de1 / de; k--; ko[k] = de2 / de; k--; ko[k] = de3 / de;
	if (b) delete[]b; 
	for (k = 0; k < kem; k++) { AA = A[k]; if (AA) delete[]AA; } 
	return ko;
}*/