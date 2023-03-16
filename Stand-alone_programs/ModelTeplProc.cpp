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
#define vmivmf 0 //выбор метода измерений для фракции 2-0,7 мм: 0 - нестационарный, 1 - стационарный
#define vyfv 0 //выбор фракции: 0 - фракция 2-0,7 мм, 1 - фракция 8-4 мм, 2 - фракция 1,6-0,35 мм
#define vyuv 1 //выбор укладки: 1 - плоскопараллельная, 2 - вертикальная
#define vysv 0 //выбор состояния: 0 - исходное, 1 - после повторных измерений, 2 - после прокаливания при 1000 град С
using namespace std;
const int cemn=11, dso=12, dmko=4, isrp=0, vpkf=0, vpmf=0, cemdu=6, cemdum=2; //dmko - длина массива коэффициентов
const int vybitom=0, vybkvi=4;
const double ys0=3e1*1e-3, dete=1e2, te0=273.15, tn=22.0+te0, tnac=2e2+te0, templa=134.0*1e1+te0; 
const double pori440=(8e1+82e0)*1e-2/2e0, pori620=(75e0+78e0)*1e-2/2e0, pori860=(65e0+68e0)*1e-2/2, pori1000=(62e0+65e0)*1e-2/2e0;
struct derevo 
{
	int otre; //1 - отражение или 0 - пропускание
	int ste; //номер стенки
	int vis; //видимость: 1 - виден, 0 - нет
	int lev; //номер уровня
	struct derevo *back; //указатель назад
	struct derevo *next; //указатель вперед
};
class KoePog 
{ 
public: double alp, tem; 
KoePog *nex; 
};
//----------------------------------------
double *qob = NULL, *ete = NULL, tna=0.0, *tgor = NULL, *thol = NULL;
double por = 0.0, *tsred=NULL, *ektpv = NULL, *kektp = NULL, porsha=0.0;
char *snm = NULL;
int cem=cemn, vybsha=0, vystsha=0;
//----------------------------------------
double **oprkoefKTPiskhchao(int, int, double *, double, int, double **, int, int, int);
double **opredtemphc(double *, double *, double *, double *, double *, int, int, double, double *, double *, double *, double **, char *);
double **vydelPol(int, int, double **, double **, int, int);
double **opredTempHolGor(double *, double *, int, int, double, int, double **, int, int, double *, int);
double **napMasEKTPVer(int, int, int, double *, int, double, int, int, int, double **, int, int);
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
double *danIskh207(double *, double*, int, int, int);
double **PoiskZavVelTem(int, double **, int, int, double *, double);
double **napMasEKTPVerNac();
double *koefPribSha(double *, double *, int, double *, char *);
double novNapMas(double, int, int, int, int, int, int, int);
void vyvodfile(double *, int, int, double, char *);
double *poisMasKoefItom(int, double *, int);
int napMasEKTPitomNac();
int napMasEKTPkviNac();
double *poisMasKoefkvi(int, double *, int);
double *poisMasKoefsha(double *, int);
int napMasEKTPshaNac();
//----------------
double novNapMas(double tnach, int vyve, int vymave, int vmimfv, int vyfrve, int vyukve, int vysove, int cemv)
{	
	int k=0;
	snm=new char[dso]; if (!snm) { cout << "No memory"; k=getchar(); exit(1); } 
	k=0; snm[k]='N'; k++; snm[k]='o'; k++; snm[k]='_'; k++; snm[k]='m'; k++; snm[k]='e'; k++; 
	snm[k]='m';  k++; snm[k]='o'; k++; snm[k]='r'; k++; snm[k]='y'; k++; snm[k]='!'; k++; snm[k]='\0';
	if (!vyve) //если выбран шамот
	{
	double ko=1e-2, porsha0=21.8*ko, porsha1=11.014*ko, porsha2=25.2*ko, porsha3=26.5*ko, porsha4=11.5*ko, porsha5=16.5*ko; //0 - Роучка, 1 - ШПД-41, 2 - ШБ-1 2-1, 3 - ШВ-1 1-1, 4 - ШПД, 5 - ШКУ-32 3-1
	if (!vymave) por=porsha0; else if (vymave==1) por=porsha1; else if (vymave==2) por=porsha2;
	else if (vymave==3) por=porsha3; else if (vymave==4) por=porsha4; else if (vymave==5) por=porsha5;
	}
	else if (vyve==1) 
	{ //если выбран вермикулит
	double ko=1e-2, por207=66.35*ko, poris84=55.75*ko, porin84=81.53*ko;
	double poro84=86.61*ko, poro16035=84.36*ko, por16035=83.97*ko;
	if ((!vysv) || (vysv==1)) 
	{  //исходный или после повторных измерений
	if (!vyfv) 
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
	if (!vymave) por=pori440; 
	else if (vymave==1) por=pori620; 
	else if (vymave==2) por=pori860; 
	else if (vymave==3) por=pori1000;
	}
	else if (vyve==3) 
	{ //если выбран КВИ
	double ko=1e-2, pork400=51.74*ko, pork500=52.07*ko, pork600=51.51*ko, pork700=39.75*ko;
	double pork800=40.85*ko, pork900=39.37*ko, pork1000=36.07*ko;
	if (vymave==4) por=pork400; 
	else if (vymave==5) por=pork500; 
	else if (vymave==6) por=pork600; 
	else if (vymave==7) por=pork700; 
	else if (vymave==8) por=pork800; 
	else if (vymave==9) por=pork900;
	else if (vymave==10) por=pork1000;
	}
	if (cemv>0) ete=new double[cemv]; if (!ete) { cout << snm << endl; k = getchar(); exit(1); }
	k=0; ete[k]=tnach; 
	for (k=1; k<cemv; k++) ete[k]=ete[k-1]+dete; //cout << "por = " << porver << endl;
	return por;
}
double **PoiskZavVelTem(int v, double **mu, int rmu, int n, double *efte, double h)
{
	int k=1, f=rmu, j=0, q=j, w=j, b=10, l=0; double *poin=new double[k], cf=0.0, hf=1e0, nf=cf, tf=cf, e=1e-1, cfp=cf, cft=cf;
	double *temvht=NULL, *temvct=NULL, *tepvt=NULL, *po=NULL, *temhq=NULL, *ts=NULL, ***muu=NULL;
	double *temcq=NULL, *temvs=NULL, *tepv=NULL, *cemt=NULL, *ktpq=NULL, *ktp=NULL, **muv=NULL;
	int vyfrve=0, vysove=0, vymivmf=1, vyukve=1, c=0, nnvyfv=2, nnvysv=3, nnvmivmf=2;
	int nnvyuv=2, kvf=0, jvsv=0, qvmi=0, qvuv=0, d=0;
	double *nvyfv=new double[nnvyfv], *nvysv=new double[nnvysv], *nvmivmf=new double[nnvmivmf];
	double *nvyuv=new double[nnvyuv]; temhq=new double[n]; temcq=new double[n]; temvs=new double[n]; 
	tepv=new double[n]; cemt=new double[n]; ktpq=new double[n]; ktp=new double[n]; muu=new double**[b]; 
	if ((!nvyfv) || (!nvysv) || (!nvmivmf) || (!nvyuv) || (!temhq) || (!temcq) || (!temvs) || (!tepv) || (!cemt) || (!ktpq) || (!ktp) || (!muu)) 
	{ cout << snm << endl; k=getchar(); exit(1); } 
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
					muv=new double*[f]; if (!muv) { cout << snm << endl; k=getchar(); exit(1); }
					muv=napMasEKTPVer(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmko); muu[b]=muv; b++; if (vysove>0) break; } }
			else if (vyfrve==1) {
				for (qvuv=0; qvuv<nnvyuv; qvuv++) {
					vyukve=nvyuv[qvuv];
					muv=new double*[f]; if (!muv) { cout << snm << endl; k=getchar(); exit(1); }
					muv=napMasEKTPVer(vyfrve, vysove, vymivmf, efte, k, h, vyukve, n, k, muv, f, dmko); muu[b]=muv; b++; } } } } 
	for (j=0; j<b; j++) { muv=muu[j];
	d=0; temvct=muv[d]; d++; temvht=muv[d]; d++; tepvt=muv[d]; d++; ktp=muv[d]; d++; ts=muv[d]; d++; po=muv[d]; k=0; nf=po[k]; 
					cf=e; while (cf<nf) { cf=cf+hf; k++; } q=k; 
					d=1; k=0; cfp=tepvt[k]; for (k=d; k<q; k++) { cft=tepvt[k]; if ((cft<=cfp) && (d>0)) { d=-1; break; } cfp=cft; } //{ for (k=0; k<q; k++) cout << "c = " << c << "\tqo = " << tepvt[k] << "\tts = " << ts[k] << "\t"; cout << endl; } //vyvodfile(ts, q, k, nf, "C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\tmp.txt");
						if (d>0) { for (w=0; w<q; w++) { 
							tf=ts[w]; cf=tepvt[w]*tf;
							for (k=0; k<n; k++) {
							if (fabs(tf-efte[k])<=hf) {
								temvs[k]=temvs[k]+tf;
								temhq[k]=temhq[k]+tf*temvht[w]; 
								temcq[k]=temcq[k]+tf*temvct[w]; 
								ktpq[k]=ktpq[k]+tf*ktp[w];
								tepv[k]=tepv[k]+cf; 
								cemt[k]=cemt[k]+hf; //cout << "tf = " << tf << "\tptp = " << tepvt[k] << "\t";
								break; } } } } c++; }
	b=c; for (k=0; k<b; k++) { muv=muu[k]; for (j=0; j<f; j++) { po=muv[j]; delete[]po; } if (muv) delete[]muv; } //for (k=0; k<n; k++) cout << "cemt = " << cemt[k] << "\t";
		for (k=0; k<n; k++) {
			cf=temvs[k]; if (cf>e) {
				temhq[k]=temhq[k]/cf; temcq[k]=temcq[k]/cf; ktpq[k]=ktpq[k]/cf; tepv[k]=tepv[k]/cf; }
			cf=cemt[k]; if (cf>e) temvs[k]=temvs[k]/cf; else temvs[k]=0.0; }
	cf=0.0; for (k=0; k<n; k++) cf=cf+hf;
	k=0; poin[k]=cf;
	mu[k]=temcq; k++; mu[k]=temhq; k++; mu[k]=tepv; k++; 
	mu[k]=ktpq; k++; mu[k]=temvs; k++; mu[k]=poin; //for (k=0; k<n; k++) cout << "cemt = " << cemt[k] << "\tth = " << temhq[k] << "\ttc = " << temcq[k] << "\tq = " << tepv[k] << "\tktp = " << ktpq[k] << "\tts = " << temvs[k] << endl;
	if (nvyfv) delete[]nvyfv; if (nvysv) delete[]nvysv;
	if (nvmivmf) delete[]nvmivmf; if (nvyuv) delete[]nvyuv;
	if (muu) delete[]muu; if (cemt) delete[]cemt;
	return mu; 
}
//----------------
double **opredtemphc(double *effktp, double *efftem, double *tgv, double *thv, double *qon, int dlma, int n, double h, double *qob, double *tgor, double *thol, double **mu, char *snome) //n=3 - длина массива коэффициентов приближающего многочлена, tgv - температура горячей стенки, thv - температура холодной стенки, dlma - длина массива ЭКТП
{
	int k=0, j=0; double *kq=new double[n], *kho=new double[n], *kgo=new double[n];
	double *tesr=new double[n], ts=0.0, g=0.0, p=0.0, hf=1e0, r=0.0, s=0.0, tego=0.0, teho=0.0;
	if ((!kq) || (!kho) || (!kgo) || (!tesr)) { cout << snome << endl; j=getchar(); exit(1); }
	for (k=0; k<n; k++) tesr[k]=(tgv[k]+thv[k])/2e0;
	kq=koefPribSha(qon, tesr, n, kq, snm);
	kgo=koefPribSha(tgv, tesr, n, kgo, snm);
	kho=koefPribSha(thv, tesr, n, kho, snm);
	for (k=0; k<dlma; k++) 
	{
		ts=efftem[k];
		g=0.0; p=0.0; for (j=0; j<n; j++) { g=g+pow(ts, p)*kq[j]; p=p+hf; } 
		if (g<0.0) g=0.0; if (qob) qob[k]=g;
		p=effktp[k];
		g=fabs(g*h/p); 
		tego=ts+g/2e0; teho=ts-g/2e0;
		if ((tego<0.0) || (teho<0.0)) 
		{
		if (tego<0.0) 
		{
			r=0.0; s=0.0; for (j=0; j<n; j++) { r=r+pow(ts, s)*kgo[j]; s=s+hf; } 
			if (r<0.0) r=0.0; tego=r; 
		}
		if (teho<0.0) 
		{
			r=0.0; s=0.0; for (j=0; j<n; j++) { r=r+pow(ts, s)*kho[j]; s=s+hf; } 
			if (r<0.0) r=0.0; teho=r; 
		} 
		}
	if (tgor) tgor[k]=tego; if (thol) thol[k]=teho; }
	if (tesr) delete[]tesr; if (kq) delete[]kq; 
	if (kho) delete[]kho; if (kgo) delete[]kgo;
	k=0; mu[k]=thol; k++; mu[k]=tgor; k++; mu[k]=qob; k++; mu[k]=effktp; k++; mu[k]=efftem;
	return mu;
}
//----------------
double **opredTempHolGor(double *ektp, double *ete, int n, int l, double h0, int v, double **mu, int rmu, int ni, double *efte, int dmkoef) //моделирование процесса теплообмена в образце
{ //n - длина массива ektp, l - длина массивов temvs, qob, ktp, temvh, temvc, ni - длина efte, qob - плотность теплового потока, которую может создать лабораторная установка 
	int cemf=rmu, k=0, j=1, w=0, rm=k, nk=n; 
	double g=0.0, p=0.0, *qon=NULL, *tena=NULL, nit=1e10, hf=1e0, *po=new double[j], del=hf;
	double ep=1e-3, d=1e-4, Thna=tn, Thnac=0.0, *koeq=new double[dmkoef], ts=0.0, laef=0.0;
	double tgor=0.0, thol=0.0, tsred=0.0, etem=0.0, dt=hf, **muv=new double*[cemf], nf=0.0, s=0.0;
	if ((!koeq) || (!po) || (!muv)) { cout << snm << endl; k=getchar(); exit(1); }
	muv=PoiskZavVelTem(k, muv, cemf, ni, efte, h0);
	k=2; qon=muv[k]; k=4; tena=muv[k]; k++; po=muv[k];
	k=0; nf=po[k]; w=k; g=ep; while (g<nf) { g=g+hf; w++; } //cout << "w = " << w << endl; for (k=0; k<w; k++) cout << "q ( " << k << " ) = " << qon[k] << "\tt = " << tena[k] << endl;
	for (k=0; k<dmkoef; k++) koeq[k]=0.0; koeq=koefPribSha(qon, tena, w, koeq, snm);
	double *qob=new double[n], *temvs=new double[n], *temvc=new double[n], *temvh=new double[n], *ktp=new double[n];
	if ((!qob) || (!temvs) || (!temvc) || (!temvh) || (!ktp)) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<n; k++) { qob[k]=0.0; temvs[k]=0.0; temvc[k]=0.0; temvh[k]=0.0; ktp[k]=0; }
		for (k=0; k<ni; k++) {
			ts=ete[k]; g=0.0; p=g;
			for (j=0; j<dmkoef; j++) { g=g+pow(ts, p)*koeq[j]; p=p+hf; } 
			qob[k]=g; } //cout << "n = " << n << endl; for (k=0; k<n; k++) cout << "qob = " << qob[k] << "\tektp = " << ektp[k] << "\tte = " << ete[k] << endl; 
		if (koeq) delete[]koeq; for (k=0; k<cemf; k++) { po=muv[k]; delete[]po; } delete[]muv; 
		g=0.0; for (k=0; k<n; k++) {
			s=qob[k];
			if (s>ep) {
				laef=ektp[k]; p=0.0; Thnac=Thna+g*dt; del=hf; etem=ete[k]; ktp[k]=laef;
				while ((del>ep) && (p<nit)) {
					thol=Thnac+p*d; //Tg - массив температур горячих стенок, Th - массив температур холодных стенок
					tgor=thol+qob[k]*h0/laef;
					del=fabs(2e0*etem-(thol+tgor));
					p=p+hf; } 
				g=g+hf; }
			else { thol=0.0; tgor=0.0; qob[k]=0.0; ktp[k]=0.0; } tsred=(tgor+thol)/2e0; 
				temvc[k]=thol; temvh[k]=tgor; temvs[k]=tsred; }
	g=0.0; for (k=0; k<n; k++) g=g+hf; //cout << "n = " << g << endl; for (k=0; k<n; k++) cout << "qob = " << qob[k] << "\tektp = " << ektp[k] << "\tte = " << ete[k] << "\tth = " << temvh[k] << "\ttc = " << temvc[k] << endl;
	k=1; po=new double[k]; muv=new double*[cemf]; if ((!po) || (!muv)) { cout << snm << endl; k=getchar(); exit(1); } k=0; po[k]=g;
	muv[k]=temvc; k++; muv[k]=temvh; k++; muv[k]=qob; k++; muv[k]=ktp; k++; muv[k]=temvs; k++; muv[k]=po; for (k=0; k<cemf; k++) mu[k]=muv[k]; //for (k=0; k<n; k++) cout << "tc = " << temvc[k] << "\ttg = " << temvh[k] << "\tts = " << temvs[k] << endl; 
	k=0; j=1; mu=vydelPol(k, n, muv, mu, rmu, j); //cout << muv << "\t" << mu << "\t"; //temvc=mu[k]; k++; temvh=mu[k]; k++; qob=mu[k]; k++; ktp=mu[k]; k++; temvs=mu[k]; k++; po=mu[k]; k=0; g=po[k]; nf=ep; while (nf<g) { nf=nf+hf; k++; } nk=k; for (k=0; k<nk; k++) cout << "tc = " << temvc[k] << "\ttg = " << temvh[k] << "\tts = " << temvs[k] << endl; for (k=0; k<cemf; k++) { po=muv[k]; if (po) delete[]po; }
	if (temvc) delete[]temvc; if (temvh) delete[]temvh; //if (qob) delete[]qob; 
	if (ktp) delete[]ktp; if (temvs) delete[]temvs; if (po) delete[]po; if (muv) delete[]muv; 
	return mu;
}
//----------------
int napMasEKTPitomNac()
{
	int vpi=vybitom, k=0, j=0, f=cemdu, cemi=cem, dmkoi=dmko, vyve=2;
	double *kektpi = new double[dmkoi]; if (!kektpi) { cout << snm << endl; k=getchar(); exit(1); }
	kektpi=poisMasKoefItom(vpi, kektpi, dmkoi);
	por=novNapMas(tnac, vyve, vpi, vmivmf, vyfv, vyuv, vysv, cemi);
	double hf = 1e0, g=0.0, s=0.0, t=0.0, *po=NULL, nf = 0.0, tn=0.0;
	double *ektpi=new double[cemi], **mu=new double*[f], **muv=new double*[f];
	double *qob=NULL, *thol=NULL, *tgor=NULL, *tsred=NULL;
	if ((!ektpi) || (!mu) || (!muv)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<f; k++) { po=new double[cemi]; if (!po) { cout << snm << endl; k=getchar(); exit(1); } muv[k]=po; } 
	for (k = 0; k < cemi; k++) { t=ete[k]; s=0.0; g=0.0; 
	for (j = 0; j<dmkoi; j++) { s = s + kektpi[j] * pow(t, g); g = g + hf; } ektpi[k] = s; } //for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpi[k] << "\t"; cout << endl;
	muv=opredTempHolGor(ektpi, ete, cemi, k, ys0, k, muv, f, cemi, ete, dmkoi); //cem - длина массива efte
	k=0; j=1; mu=vydelPol(k, cemi, muv, mu, f, j);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; 
	ektpi=mu[k]; k++; tsred=mu[k]; ete=tsred; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cemi=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cemi; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpi[k] << "\tth = " << thol[k] << "\ttg = " << tgor[k] << "\tqo = " << qob[k] << "\t"; cout << endl;
	for (k=0; k<f; k++) { po=muv[k]; delete[]po; } for (k=0; k<f; k++) { po=mu[k]; delete[]po; }
	if (muv) delete[]muv; if (mu) delete[]mu; if (kektpi) delete[]kektpi;
	return 0;
}
int napMasEKTPkviNac()
{
	int vpk=vybkvi, k=0, j=0, f=cemdu, cemk=cem, dmkok=dmko, vyve=3;
	double *kektpk=new double[dmkok]; if (!kektpk) { cout << snm << endl; k=getchar(); exit(1); }
	kektpk=poisMasKoefkvi(vpk, kektpk, dmkok);
	por=novNapMas(tnac, vyve, vpk, vmivmf, vyfv, vyuv, vysv, cemk);
	double hf=1e0, g=0.0, s=0.0, t=0.0, *po=NULL, nf = 0.0, tn=0.0;
	double *ektpk=new double[cemk], **mu=new double*[f], **muv=new double*[f];
	double *qob=NULL, *thol=NULL, *tgor=NULL, *tsred=NULL;
	if ((!ektpk) || (!mu) || (!muv)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<f; k++) { po=new double[cemk]; if (!po) { cout << snm << endl; k=getchar(); exit(1); } muv[k]=po; } 
	for (k=0; k<cemk; k++) { t=ete[k]; s=0.0; g=0.0; 
	for (j=0; j<dmkok; j++) { s=s+kektpk[j]*pow(t, g); g=g+hf; } ektpk[k]=s; } //for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpi[k] << "\t"; cout << endl;
	muv=opredTempHolGor(ektpk, ete, cemk, k, ys0, k, muv, f, cemk, ete, dmkok); //cem - длина массива efte
	k=0; j=1; mu=vydelPol(k, cemk, muv, mu, f, j);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; 
	ektpk=mu[k]; k++; tsred=mu[k]; ete=tsred; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cemk=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cemk; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpk[k] << "\tth = " << thol[k] << "\ttg = " << tgor[k] << "\tqo = " << qob[k] << "\t"; cout << endl;
	for (k=0; k<f; k++) { po=muv[k]; delete[]po; } for (k=0; k<f; k++) { po=mu[k]; delete[]po; }
	if (muv) delete[]muv; if (mu) delete[]mu; if (kektpk) delete[]kektpk;
	return 0;
}
int napMasEKTPshaNac()
{
	int vps=vybsha, k=0, j=0, f=cemdu, cems=cem, dmkos=dmko, vyve=0;
	double *kektps=new double[dmkos]; if (!kektps) { cout << snm << endl; k=getchar(); exit(1); }
	kektps=poisMasKoefsha(kektps, dmkos);
	por=novNapMas(tnac, vyve, vps, vmivmf, vyfv, vyuv, vysv, cems);
	double hf=1e0, g=0.0, s=0.0, t=0.0, *po=NULL, nf = 0.0, tn=0.0;
	double *ektps=new double[cems], **mu=new double*[f], **muv=new double*[f];
	double *qob=NULL, *thol=NULL, *tgor=NULL, *tsred=NULL;
	if ((!ektps) || (!mu) || (!muv)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<f; k++) { po=new double[cems]; if (!po) { cout << snm << endl; k=getchar(); exit(1); } muv[k]=po; } 
	for (k=0; k<cems; k++) { t=ete[k]-te0; s=0.0; g=0.0; 
	for (j=0; j<dmkos; j++) { s=s+kektps[j]*pow(t, g); g=g+hf; } ektps[k]=s; } //for (k=0; k<cem; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpi[k] << "\t"; cout << endl;
	muv=opredTempHolGor(ektps, ete, cems, k, ys0, k, muv, f, cems, ete, dmkos); //cem - длина массива efte
	k=0; j=1; mu=vydelPol(k, cems, muv, mu, f, j);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; 
	ektps=mu[k]; k++; tsred=mu[k]; ete=tsred; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cems=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cems; k++) cout << "tem = " << ete[k] << "\tktp = " << ektps[k] << "\tth = " << thol[k] << "\ttg = " << tgor[k] << "\tqo = " << qob[k] << "\t"; cout << endl;
	for (k=0; k<f; k++) { po=muv[k]; delete[]po; } for (k=0; k<f; k++) { po=mu[k]; delete[]po; }
	if (muv) delete[]muv; if (mu) delete[]mu; if (kektps) delete[]kektps;
	return 0;
}
double **napMasEKTPVerNac()
{
	int k=1, nk=0, j=0, jk=0, f=cemdu, qn=0, q=0, w=f-1, cemv=cem; 
	double hf=1e0, *po=NULL, s=0.0;
	double nf=0.0, t=0.0, g=0.0, e=1e-1;
	por=novNapMas(tnac, k, j, vmivmf, vyfv, vyuv, vysv, cemv);
	double **mu=new double*[f], **muv=new double*[f]; if ((!mu) || (!muv)) { cout << snm << endl; k=getchar(); exit(1); } 
	muv=napMasEKTPVer(vyfv, vysv, vmivmf, ete, k, ys0, vyuv, cemv, k, muv, f, dmko); //for (k=0; k<cem; k++) cout << "te = " << ete[k] << "\t"; cout << endl; //k=w; po=muv[k]; k=0; nf=po[k]; t=e; k=0; while (t<nf) { t=t+hf; k++; } cem=k; cout << "cem = " << cem << "\tnf = " << nf << "\n"; for (k=0; k<w; k++) { po=muv[k]; for (j=0; j<cem; j++) cout << " ( " << j << " ) = " << po[j] << "\t"; cout << endl; }
	k=0; j=1; mu=vydelPol(k, cemv, muv, mu, f, j);
	k=0; thol=mu[k]; k++; tgor=mu[k]; k++; qob=mu[k]; k++; ektpv=mu[k]; k++; tsred=mu[k]; ete=tsred; k++; po=mu[k]; 
	k=0; nf=po[k]; t=0.0; k=0; while (t<nf) { t=t+hf; k++; } cemv=k; //cout << "cem = " << cem << "\tnf = " << nf << "\t";	
	for (k=0; k<cemv; k++) cout << "tem = " << ete[k] << "\tktp = " << ektpv[k] << "\tth = " << thol[k] << "\ttg = " << tgor[k] << "\tqo = " << qob[k] << "\t"; cout << endl;
	//for (k=0; k<f; k++) { po=muv[k]; delete[]po; } if (!((!vysv) && (!vmivmf) && (!vyfv))) for (k=0; k<f; k++) { po=mu[k]; delete[]po; }
	if (muv) delete[]muv; //if (mu) delete[]mu;
	return mu; //return 0;
}
//----------------
double **napMasEKTPVer(int vyfrve, int vysove, int vymeizvemafr, double *te, int v, double h, int vyukve, int nt, int wa, double **mu, int rmu, int dmk)
{
	int k=0, nvm=dmk, j=0, vy=0, n84=0, n207=0, cedu=rmu, q=0, cedumi=cemdum, nn=0;
	double F=13.85*1e-4, nf84=0.0, nf207=0.0, *po=NULL, *tc84=NULL, s=0.0, t=0.0, e=1e-1; 
	double *th84=NULL, hf=1e0, *tp84=NULL, *th207=NULL, *tc207=NULL, *tp207=NULL;
	double *koeft=NULL, *koefh=NULL, *koefc=NULL, *koefq=NULL, *koefs=NULL, **muv=NULL, r=0.0; //cout << "vf = " << vyfrve << "\tvs = " << vysove << "\tvmi = " << vymeizvemafr << "\tvu = " << vyukve << endl; //for (k=0; k<nt; k++) cout << "te ( " << k << " ) = " << te[k] << "\t";
	double *ktp=NULL, *temvs=NULL, *temvh=NULL, *temvc=NULL, *tepv=NULL, *vm=NULL, *ma=NULL; 
	//-------
	if ((!vyfrve) || (vyfrve==2)) { //фракция 2-0,7 мм или фракция 1,6-0,35 мм
		th207=arrTemHigh207(); 
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snm << endl; k = getchar(); exit(1); }
		muv=arrTemCold207(muv); k=0; po=muv[k]; k++; tc207=muv[k]; 
		tp207=arrTepPot207();
		k=0; t=0.0; nf207=po[k]; while (t<nf207) { t=t+hf; k++; } n207=k; nn=k; //cout << "n207 = " << n207 << "\t";
		temvs=new double[n207]; temvh=new double[n207]; temvc=new double[n207]; 
		tepv=new double[n207]; ktp=new double[n207]; 
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snm << endl; k=getchar(); exit(1); } //for (k=0; k<n207; k++) cout << "th = " << th207[k] << "\ttc = " << tc207[k] << "\tQ = " << tp207[k] << endl;
		for (k=0; k<n207; k++) {
		temvh[k]=th207[k]; temvc[k]=tc207[k]; temvs[k]=(temvh[k]+temvc[k])/2e0; 
		tepv[k]=tp207[k]/F; ktp[k]=fabs(temvh[k]-temvc[k])/h; ktp[k]=tepv[k]/ktp[k]; }
		if (tc207) delete[]tc207; if (th207) delete[]th207; if (tp207) delete[]tp207; if (po) delete[]po; //for (k=0; k<n207; k++) cout << "ktp = " << ktp[k] << "\tqob = " << tepv[k] << endl;
		if (!vysove) {
		if (!vymeizvemafr) { //установка Netzsch
		if (muv) delete[]muv; muv=new double*[cedu]; if (!muv) { cout << snm << endl; k=getchar(); exit(1); }
		muv=oprkoefKTPiskhchao(vymeizvemafr, k, te, h, n207, muv, cedu, nt, dmk); 
		k=0; koefc=muv[k]; k++; koefh=muv[k]; k++; koefq=muv[k]; k++;
		koeft=muv[k]; k++; koefs=muv[k]; k++; po=muv[k]; }
		else if (vymeizvemafr==1) { //данные 2020 года - ГОСТ 12170 - стационарный метод
		if (ktp) delete[]ktp; ktp=arrKTP_2020(); if (temvh) delete[]temvh; temvh=arrTem1_2020(); 
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snm << endl; k=getchar(); exit(1); }
		muv=arrTem2_2020(muv);
		k=0; po=muv[k]; nf207=po[k]; k++; if (temvc) delete[]temvc; temvc=muv[k]; 
		if (temvs) delete[]temvs; temvs=arrTem3_2020(); 
		k=0; s=e; while (s<nf207) { s=s+hf; k++; } n207=k; k=0; q=3; //cout << "nk 2020 = " << n207 << endl; //for (k=0; k<n207; k++) cout << "ts ( " << k << " ) = " << temvs[k] << "\tktp = " << ktp[k] << "\tth = " << temvh[k] << "\ttc = " << temvc[k] << endl;
		double *ts=new double[n207], *tst=new double[n207]; if ((!ts) || (!tst)) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<n207; k++) { tst[k]=(temvh[k]+temvc[k])/2e0; ts[k]=tst[k]; temvs[k]=ts[k]; }
		ts=danIskh207(ts, tst, k, q, q); ktp=danIskh207(ktp, tst, k, q, q); temvh=danIskh207(temvh, tst, k, q, q); 
		temvc=danIskh207(temvc, tst, k, q, q); temvs=danIskh207(temvs, tst, k, q, q); if (tst) delete[]tst; //for (k=0; k<q; k++) cout << "ts ( " << k << " ) = " << ts[k] << "\tktp = " << ktp[k] << "\tth = " << temvh[k] << "\ttc = " << temvc[k] << endl;
		koeft=new double[dmk]; koefh=new double[dmk];
		koefc=new double[dmk]; koefs=new double[dmk]; koefq=new double[dmk];
		if ((!koeft) || (!koefh) || (!koefc) || (!koefs) || (!koefq)) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; } 
		if (ktp) koeft=koefPribSha(ktp, ts, q, koeft, snm); 
		if (temvh) koefh=koefPribSha(temvh, ts, q, koefh, snm);
		if (temvc) koefc=koefPribSha(temvc, ts, q, koefc, snm);
		if (temvs) koefs=koefPribSha(temvs, ts, q, koefs, snm);
		tepv=new double[q]; if (!tepv) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<q; k++) {
			s=0.0; r=0.0; for (j=0; j<dmk; j++) { s=s+koeft[j]*pow(ts[k], r); r=r+hf; }
			tepv[k]=s*fabs(temvh[k]-temvc[k])/h; }
			koefq=koefPribSha(tepv, ts, q, koefq, snm); }
		else if (vymeizvemafr==2) { //данные 2019 года - ГОСТ 12170
		koefq=danPoTemTepl2071(temvs, tepv, nvm); koefh=danPoTemTepl2071(temvs, temvh, nvm);
		koefc=danPoTemTepl2071(temvs, temvc, nvm); koeft=danPoTemTepl2071(temvs, ktp, nvm);
		koefs=danPoTemTepl2071(temvs, temvs, nvm); } }
				if (vysove==1) { //после повторных измерений
			koefc=danPoTemTepl2072(temvs, temvc, nvm); koefh=danPoTemTepl2072(temvs, temvh, nvm); 
			koefq=danPoTemTepl2072(temvs, tepv, nvm); koeft=danPoTemTepl2072(temvs, ktp, nvm); 
			koefs=danPoTemTepl2072(temvs, temvs, nvm); }
				else if (vysove == 2) { //после обжига при 1000 °С
			koefc=danPoTemTepl2073(temvs, temvc, nvm); koefh=danPoTemTepl2073(temvs, temvh, nvm); 
			koefq=danPoTemTepl2073(temvs, tepv, nvm); koeft=danPoTemTepl2073(temvs, ktp, nvm); 
			koefs=danPoTemTepl2073(temvs, temvs, nvm); }
				else if (vysove == 3) 
				{ //после повторного обжига при 1000 °С
			koefc=danPoTemTepl2074(temvs, temvc, nvm); koefh=danPoTemTepl2074(temvs, temvh, nvm); 
			koefq=danPoTemTepl2074(temvs, tepv, nvm); koeft=danPoTemTepl2074(temvs, ktp, nvm); 
			koefs=danPoTemTepl2074(temvs, temvs, nvm); } }
	else if (vyfrve==1) { //фракция 8-4 мм
		k=cedumi; muv=new double*[k]; muv=arrTemCold84(muv); k=0; po=muv[k]; k++; tc84=muv[k]; 
		th84=arrTemHigh84(); tp84=arrTepPot84();
		k=0; t=0.0; nf84=po[k]; while (t<nf84) { t=t+hf; k++; } n84=k; nn=k; 
		temvs=new double[n84]; temvh=new double[n84]; temvc=new double[n84]; tepv=new double[n84]; ktp=new double[n84];
		if ((!temvs) || (!temvh) || (!temvc) || (!tepv) || (!ktp)) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<n84; k++) { 
		temvs[k]=(th84[k]+tc84[k])/2e0; 
		temvh[k]=th84[k]; temvc[k]=tc84[k]; tepv[k]=tp84[k]; 
		tepv[k]=tepv[k]/F; ktp[k]=fabs(th84[k]-tc84[k])/h; 
		ktp[k]=tepv[k]/ktp[k]; } //for (k=0; k<n84; k++) cout << "th = " << temvh[k] << "\ttc = " << temvc[k] << "\ttp = " << tepv[k] << endl;
		if (po) delete[]po; if (th84) delete[]th84; if (tp84) delete[]tp84; if (tc84) delete[]tc84; 
		if (vyukve==1) { //плоско-параллельная засыпка
			if (!vysove) { //исходный
				koefc=danPoTemTepl840(temvs, temvc, nvm); koefh=danPoTemTepl840(temvs, temvh, nvm); 
				koefq=danPoTemTepl840(temvs, tepv, nvm); koeft=danPoTemTepl840(temvs, ktp, nvm); 
				koefs=danPoTemTepl840(temvs, temvs, nvm); }
			else if (vysove==1) { //после повторных измерений
				koefc=danPoTemTepl842(temvs, temvc, nvm); koefh=danPoTemTepl842(temvs, temvh, nvm); 
				koefq=danPoTemTepl842(temvs, tepv, nvm); koeft=danPoTemTepl842(temvs, ktp, nvm); 
				koefs=danPoTemTepl842(temvs, temvs, nvm); }
			else if (vysove==2) { //после обжига
				koefc=danPoTemTepl845(temvs, temvc, nvm); koefh=danPoTemTepl845(temvs, temvh, nvm); 
				koefq=danPoTemTepl845(temvs, tepv, nvm); koeft=danPoTemTepl845(temvs, ktp, nvm); 
				koefs=danPoTemTepl845(temvs, temvs, nvm); } }
		else if (vyukve==2) { //вертикальная засыпка
			if (!vysove) { //исходный
				koefc=danPoTemTepl841(temvs, temvc, nvm); koefh=danPoTemTepl841(temvs, temvh, nvm); 
				koefq=danPoTemTepl841(temvs, tepv, nvm); koeft=danPoTemTepl841(temvs, ktp, nvm); 
				koefs=danPoTemTepl841(temvs, temvs, nvm); }
			else if (vysove == 1) { //после повторных измерений
				koefc = danPoTemTepl844(temvs, temvc, nvm); koefh = danPoTemTepl844(temvs, temvh, nvm); 
				koefq = danPoTemTepl844(temvs, tepv, nvm); koeft = danPoTemTepl844(temvs, ktp, nvm); 
				koefs = danPoTemTepl844(temvs, temvs, nvm); }
			else if (vysove == 2) { //после обжига
				koefc = danPoTemTepl843(temvs, temvc, nvm); koefh = danPoTemTepl843(temvs, temvh, nvm); 
				koefq = danPoTemTepl843(temvs, tepv, nvm); koeft = danPoTemTepl843(temvs, ktp, nvm); 
				koefs = danPoTemTepl843(temvs, temvs, nvm); } } }
	if (temvs) delete[]temvs; if (temvh) delete[]temvh; 
	if (temvc) delete[]temvc; if (tepv) delete[]tepv; 
	if (ktp) delete[]ktp; if (muv) delete[]muv; 
	t=0.0; for (k=0; k<cedu; k++) t=t+hf; 
	k=1; muv=new double*[cedu]; po=new double[k]; if ((!muv) || (!po)) { cout << snm << endl; k=getchar(); exit(1); }
	k=0; po[k]=t; muv[k]=koefc; k++; muv[k]=koefh; k++; 
	muv[k]=koefq; k++; muv[k]=koeft; k++; muv[k]=koefs; k++; muv[k]=po; //vy=cedu-1; for (k=0; k<vy; k++) { ma=muv[k]; for (j=0; j<dmk; j++) cout << "k = " << k << "\tkoef ( " << j << " ) = " << ma[j] << "\t"; cout << endl; }
	q=cedu-1; for (vy=0; vy<q; vy++) { 
	vm=muv[vy]; ma=new double[nt]; if (!ma) { cout << snm << endl; k=getchar(); exit(1); }
	for (k=0; k<nt; k++) { t=te[k]; s=0.0; r=s;
	for (j=0; j<nvm; j++) { s=s+vm[j]*pow(t, r); r=r+hf; } 
	ma[k]=s; //if (!vy) cout << "te_s = " << t << "\tt_c = " << s << "\t"; if (vy==1) cout << "te_s = " << t << "\tt_h = " << s << "\t"; if (vy==2) cout << "te_s = " << t << "\tq = " << s << "\t"; if (vy==3) cout << "te_s = " << t << "\tktp = " << s << "\t";if (vy==4) cout << "te_s = " << t << "\tts = " << s << "\t";
	} //cout << endl;
	mu[vy]=ma; if (vm) delete[]vm; } po=muv[q]; if (po) delete[]po; if (muv) delete[]muv;
	k=1; po=new double[k]; if (!po) { cout << snm << endl; k=getchar(); exit(1); }
	s=0.0; for (k=0; k<nt; k++) s=s+hf; k=0; po[k]=s; mu[q]=po; //cout << "nt = " << nt << "\t";
	muv=new double*[cedu]; if (!muv) { cout << snm << endl; k=getchar(); exit(1); }
	k=0; j=-1; muv=vydelPol(k, nt, mu, muv, rmu, j); for (k=0; k<cedu; k++) { po=muv[k]; mu[k]=po; } //k=5; po=mu[k]; k=0; t=po[k]; //cout << "s = " << s << "\tt = " << t << "\t"; 
	if (muv) delete[]muv; 
	return mu;
}
//----------------
double **oprkoefKTPiskhchao(int vmiv, int v, double *efte, double h, int n207, double **mu, int rmu, int cem, int dmk) //vmiv - выбор метода измерений
{
	int nk=0, k=0, j=0, q=0, qn=0, nn=n207, cedumi=cemdum, w=rmu-1, f=rmu; 
	double hf=1e0, nf207=0.0, *po=NULL, *koeft=NULL, *koefh=NULL, s=0.0, t=0.0;
	double *koefc=NULL, *koefs=NULL, *koefq=NULL, **muv=NULL, **muvv=NULL, r=0.0, e=1e-1;
	double *tsv=NULL, *tgv=NULL, *thv=NULL, *qov=NULL, *ktpv=NULL, *mt=NULL, *ts=NULL; 
	if (!vmiv) {  //0 - установка Netzsch - нестационарный метод
		k=cedumi; muv=new double*[k]; if (!muv) { cout << snm << endl; k = getchar(); exit(1); }
		muv=arrTem_Netzsch(muv);
		k=0; po=muv[k]; k++; mt=muv[k];
		k=0; s=e; nf207=po[k]; while (s<nf207) { s=s+hf; k++; } nk=k; //cout << "nk = " << nk << "\t"; //nk=8 - длина массива ktpv
		ktpv=arrKTP_Netzsch(); 
		if (muv) delete[]muv; if (po) delete[]po;
		k=rmu; muv=new double*[k]; k=0; if (!muv) { cout << snm << endl; k=getchar(); exit(1); }
		muv=opredTempHolGor(ktpv, mt, nk, nn, h, k, muv, rmu, cem, efte, dmko); //cem - длина массива efte //for (k=0; k<f; k++) { po=muv[k]; delete[]po; }
		if (ktpv) delete[]ktpv; if (mt) delete[]mt;
		k=0; thv=muv[k]; k++; tgv=muv[k]; k++; qov=muv[k]; k++; ktpv=muv[k]; k++; 
		tsv=muv[k]; k++; po=muv[k]; if (muv) delete[]muv;
		q=0; nf207=po[q]; s=0.0; while (s<nf207) { s=s+hf; q++; } qn=q; 
		koefh=new double[dmk]; koefc=new double[dmk]; koefq=new double[dmk]; koeft=new double[dmk]; koefs=new double[dmk];
		if ((!koefh) || (!koefc) || (!koefq) || (!koeft) || (!koefs)) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<dmk; k++) { koefh[k]=0.0; koefc[k]=0.0; koefq[k]=0.0; koeft[k]=0.0; koefs[k]=0.0; }
		koefc=koefPribSha(thv, tsv, qn, koefc, snm); koefh=koefPribSha(tgv, tsv, qn, koefh, snm); 
		koefq=koefPribSha(qov, tsv, qn, koefq, snm); koeft=koefPribSha(ktpv, tsv, qn, koeft, snm); 
		koefs=koefPribSha(tsv, tsv, qn, koefs, snm); 
		if (thv) delete[]thv; if (tgv) delete[]tgv; if (qov) delete[]qov; 
		if (ktpv) delete[]ktpv; if (tsv) delete[]tsv; if (po) delete[]po; 
	}
	k=1; po=new double[k]; if (!po) { cout << snm << endl; k=getchar(); exit(1); }
	s=0.0; for (k=0; k<rmu; k++) s=s+hf; k=0; po[k]=s;
	k=0; mu[k]=koefc; k++; mu[k]=koefh; k++; mu[k]=koefq; k++; mu[k]=koeft; k++; mu[k]=koefs; k++; mu[k]=po;
	return mu; 
}
//----------------
double *danIskh207(double *ma, double *x, int v, int n, int np)
{
	if (ma) 
	{
		int k=0, p=0, q=0; 
		double s=0.0, t=0.0, hf=1e0, *m=new double[n];
		if (!m) { cout << snm << endl; k=getchar(); exit(1); }
		for (k=0; k<np; k++) { s=s+ma[k]*x[k]; t=t+x[k]; } s=s/t;
		k=0; m[k]=s; k++;
		p=3; m[k]=ma[p]; k++;
		p=4; q=5; m[k]=(ma[p]*x[p]+ma[q]*x[q])/(x[p]+x[q]); //600, 800 и 1000 °C
		if (ma) delete[]ma; return m;
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
double **vydelPol(int v, int n, double **mu, double **muv, int f, int fl)
{ 
	int q=0, qn=n, k=0, m=0, w=f-1, j=1, nk=n, *mui=NULL, x=0, qk=0;
	double *po=NULL, nf=0.0, t=0.0, hf=1e0, *vm=NULL, e=1e-4, r=0.0, s=0.0;
	k=w; po=mu[k]; k=0; nf=po[k]; 
	t=e; while (t<nf) { k++; t=t+hf; } nk=k; mui=new int[nk]; if (!mui) { cout << snm << endl; k=getchar(); exit(1); }
	q=0; for (k=0; k<nk; k++) {
		x=1; for (m=0; m<w; m++) { 
			po=mu[m]; if ((po[k]<e) && (x>0)) { x=-1; break; } } 
		if (x>0) { mui[k]=k; q++; } else mui[k]=-1; } qn=q;
		if (fl>0) { m=1; po=mu[m]; for (k=0; k<nk; k++) if (po[k]>templa) { mui[k]=-1; qn--; } } //if (fl>0) { cout << "qn = " << qn << "\tnf = " << nf << "\tw = " << w << "\tn = " << n << "\t"; } 
	for (m=0; m<w; m++) {
	vm=new double[qn]; if (!vm) { cout << snm << endl; k=getchar(); exit(1); } //for (k=0; k<qn; k++) vm[k]=0.0;
	q=0; po=mu[m]; for (k=0; k<nk; k++) {
		x=mui[k]; if (x>=0) { vm[q]=po[x]; q++; } } muv[m]=vm; } qn=q; qk=q; if (mui) delete[]mui; 
	nf=0.0; for (k=0; k<qn; k++) nf=nf+hf; 
	k=1; po=new double[k]; if (!po) { cout << snm << endl; k=getchar(); exit(1); } k=0; po[k]=nf; muv[w]=po; //if (fl>0) { q=w; po=muv[q]; q=0; t=po[q]; cout << "r = " << t << endl; s=e; q=0; while (s<t) { s=s+hf; q++; } qk=q; q=0; po=muv[q]; for (q=0; q<qk; q++) cout << "temc = " << po[q] << "\t"; cout << endl; q=1; po=muv[q]; for (q=0; q<qk; q++) cout << "temh = " << po[q] << "\t"; cout << endl; q=2; po=muv[q]; for (q=0; q<qk; q++) cout << "qo = " << po[q] << "\t"; cout << endl; q=3; po=muv[q]; for (q=0; q<qk; q++) cout << "ktp = " << po[q] << "\t"; cout << endl; q=4; po=muv[q]; for (q=0; q<qk; q++) cout << "tems = " << po[q] << "\t"; cout << endl; q=w; po=mu[q]; q=0; t=po[q]; cout << "r = " << t << endl; s=e; q=0; while (s<t) { s=s+hf; q++; } qk=q; q=0; po=mu[q]; for (q=0; q<qk; q++) cout << "temc = " << po[q] << "\t"; cout << endl; q=1; po=mu[q]; for (q=0; q<qk; q++) cout << "temh = " << po[q] << "\t"; cout << endl; q=2; po=mu[q]; for (q=0; q<qk; q++) cout << "qo = " << po[q] << "\t"; cout << endl; q=3; po=mu[q]; for (q=0; q<qk; q++) cout << "ktp = " << po[q] << "\t"; cout << endl; q=4; po=mu[q]; for (q=0; q<qk; q++) cout << "tems = " << po[q] << "\t"; cout << endl;  cout << mu << "\t" << muv << "\t"; }
	for (k=0; k<f; k++) { vm=mu[k]; //if ((fl>0) && (vm) && (k<w)) { for (q=0; q<nk; q++) cout << "k = "<< k << "\tq = " << q << "\tp = " << vm[q] << "\t"; cout << endl; } 
	if ((fl<0) && (vm)) delete[]vm; } //if (fl>0) { cout << "qn = " << qn << "\tn = " << nf << "\t"; } 
	return muv;
}
double *koefPribSha(double *ktp, double *te, int le, double *ko, char *snome)
{
	int k, kem=3; double **A=new double*[kem], *AA, *b=new double[kem];
	double yx2=0.0, yx=0.0, p=0.0, hf=1e0, x4=0.0, x3=0.0, x2=0.0;
	double x=0.0, y=0.0, de=0.0, de1=0.0, de2=0.0, de3=0.0;
	if (A) 
	{ 
	for (k=0; k<kem; k++) 
	{ 
	AA=new double[kem]; 
	if (AA) A[k]=AA; else { cout << snome << endl; k=getchar(); exit(1); } 
	} 
	}
	else { cout << snome << endl; k=getchar(); exit(1); }
	if ((!ko) || (!b)) { cout << snome << endl; k=getchar(); exit(1); }
	for (k=0; k<le; k++) 
	{
		yx2=yx2+ktp[k]*pow(te[k], 2e0); 
		yx=yx+ktp[k]*te[k]; 
		y=y+ktp[k];
		x4=x4+pow(te[k], 4e0); 
		x3=x3+pow(te[k], 3e0); 
		x2=x2+pow(te[k], 2e0); 
		x=x+te[k];
		p=p+hf;
	} //применение метода наименьших квадратов
	k=0; b[k]=yx2; k++; b[k]=yx; k++; b[k]=y; 
	A[0][0] = x4; A[0][1] = x3; A[0][2] = x2; A[1][0] = x3; A[1][1] = x2; A[1][2] = x; A[2][0] = x2; A[2][1] = x; A[2][2] = p;
	de = A[0][0] * (A[2][2] * A[1][1] - A[2][1] * A[1][2]) - A[0][1] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	de1 = b[0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) - A[0][1] * (b[1] * A[2][2] - b[2] * A[1][2]) + A[0][2] * (b[1] * A[2][1] - b[2] * A[1][1]);
	de2 = A[0][0] * (b[1] * A[2][2] - b[2] * A[1][2]) - b[0] * (A[1][0] * A[2][2] - A[2][0] * A[1][2]) + A[0][2] * (A[1][0] * b[2] - A[2][0] * b[1]);
	de3 = A[0][0] * (A[1][1] * b[2] - A[2][1] * b[1]) - A[0][1] * (A[1][0] * b[2] - A[2][0] * b[1]) + b[0] * (A[1][0] * A[2][1] - A[2][0] * A[1][1]);
	k=2; ko[k] = de1 / de; k--; ko[k] = de2 / de; k--; ko[k] = de3 / de;
	if (b) delete[]b; 
	for (k = 0; k < kem; k++) { AA = A[k]; if (AA) delete[]AA; } 
	return ko;
}
void vyvodfile(double *ao, int sao, int zf, double hvv, char *nafa)
{
	int k; FILE *fo = fopen(nafa, "at");
	if (!fo) { cout << "No memory" << endl; k = getchar(); exit(1); }
	fprintf(fo, "%0.6lf\n", hvv);
	for (k = 0; k < sao; k++) fprintf(fo, "%0.6lf\n", ao[k]); fprintf(fo, "\n"); fclose(fo);
}
double *poisMasKoefItom(int no, double *kti, int n)
{
	int f=n, k=0; double tei0=te0, t1=2e2, t2=38e1, kk=1e-2, te200=t1+tei0, te380=t2+tei0, *ktpit=NULL;
	if (!no) {
		double *kitom440 = new double[f]; ktpit = new double[f];
		if ((!kitom440) || (!ktpit)) { cout << snm << endl; getchar(); exit(1); }
		for (k=0; k<f; k++) { ktpit[k]=0.0; kitom440[k]=0.0; }
		k=0; ktpit[k] = 9e0*kk; k++; ktpit[k] = 12e0*kk;
		kitom440[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); 
		kitom440[0] = ktpit[0] - kitom440[1] * te200;
		for (k=0; k<f; k++) kti[k] = kitom440[k]; if (kitom440) delete[]kitom440;
	}
	else if (no == 1) {
		double *kitom620 = new double[f]; ktpit = new double[f];
		if ((!kitom620) || (!ktpit)) { cout << snm << endl; getchar(); exit(1); } 
		for (k=0; k<f; k++) { ktpit[k]=0.0; kitom620[k]=0.0; }
		//ktpit[0] = 12.0*1e-2; ktpit[1] = 139.0*1e-3; //из Диссертации
		k=0; ktpit[k] = 18.0*kk; k++; ktpit[k] = 19.0*kk; //Данные 2017 года
		kitom620[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); 
		kitom620[0] = ktpit[1] - kitom620[1] * te380;
		for (k=0; k<f; k++) kti[k] = kitom620[k]; if (kitom620) delete[]kitom620;
	}
	else if (no == 2) {
		double *kitom860 = new double[f]; ktpit = new double[f];
		if ((!kitom860) || (!ktpit)) { cout << snm << endl; getchar(); exit(1); } 
		for (k=0; k<f; k++) { ktpit[k]=0.0; kitom860[k]=0.0; }
		//ktpit[0] = 18.3*kk; ktpit[1] = 19.4*kk; //из Диссертации
		ktpit[0] = 26.0*kk; ktpit[1] = 37.0*kk; //Данные 2017 года
		kitom860[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); kitom860[0] = ktpit[1] - kitom860[1] * te380;
		for (k=0; k<f; k++) kti[k] = kitom860[k]; if (kitom860) delete[]kitom860;
	}
	else if (no == 3) {
		double *kitom1000 = new double[f]; ktpit = new double[f];
		if ((!kitom1000) || (!ktpit)) { cout << snm << endl; getchar(); exit(1); } 
		for (k=0; k<f; k++) { ktpit[k]=0.0; kitom1000[k]=0.0; }
		//ktpit[0] = 23.0*kk; ktpit[1] = 25.0*kk; //из Диссертации
		ktpit[0] = 42.0*kk; ktpit[1] = 52.0*kk; //Данные 2017 года
		kitom1000[1] = (ktpit[1] - ktpit[0]) / (te380 - te200); kitom1000[0] = ktpit[1] - kitom1000[1] * te380;
		for (k=0; k<f; k++) kti[k] = kitom1000[k]; if (kitom1000) delete[]kitom1000;
	}
	else { cout << "Net takoy marki ITOM!" << endl; k = getchar(); exit(1); }
	if (ktpit) delete[]ktpit;
	return kti;
}
double *poisMasKoefkvi(int vyb, double *kktp, int n)
{
	int f=n, k=0, j=0; double t1=25e0, t2=5e2, dt=t2-t1, kn=0.0;
	for (k=0; k<f; k++) kktp[k]=0.0;
	if (vyb==3) { kktp[1]=0.00015; kktp[0] = 0.068; } //350
	else if (vyb == 4) { kktp[1] = 0.000125; kktp[0] = 0.082; //1
	kn=(0.156-0.087)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.156-kn; //2
	kn=(0.155-0.087)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.155-kn; //3 
	kktp[1] = 0.00016; kktp[0] = 0.14; //4  //Spirina //выбор
	//kn=(0.146-0.087)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.146-kn; //к вопросу о стандартизации КВИ
	} //400
	else if (vyb == 5) { kktp[1] = 0.0001; kktp[0] = 0.103; //из Ахтямова, 1991
	kn=(0.165-0.105)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.165-kn; //выбор
	//kn=(0.178-0.105)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.178-kn; //Двухслойные
	//kn=(0.152-0.105)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.152-kn; //к вопросу о стандартизации КВИ
	} //500
	else if (vyb == 6) { kktp[1] = 0.00015; kktp[0] = 0.116; //выбор
	kn=(0.201-0.12)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.201-kn;
	kn=(0.195-0.12)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.195-kn;
	kktp[1] = 0.00015; kktp[0] = 0.17; //Spirina
	kn=(0.196-0.12)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.196-kn; //Двухслойные //выбор
	//kn=(0.178-0.12)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.178-kn; //к вопросу о стандартизации КВИ
	} //600
	else if (vyb == 7) { kktp[1] = 0.00017; kktp[0] = 0.146; 
	kn=(0.216-0.15)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.216-kn; 
	kn=(0.235-0.15)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.235-kn; //к вопросу о стандартизации КВИ
	kn=(0.251-0.16)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.251-kn; //Двухслойные //выбор
	} //700
	else if (vyb == 8) { kktp[1] = 0.00018;  kktp[0] = 0.156; 
	kn=(0.226-0.16)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.226-kn; //выбор
	//kn=(0.25-0.16)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.25-kn;
	//kktp[1] = 0.00014; kktp[0] = 0.21; //Spirina
	//kn=(0.23-0.16)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.23-kn; //к вопросу о стандартизации КВИ
	} //800
	else if (vyb == 9) { kktp[1] = 0.00019;  kktp[0] = 0.185; 
	kn=(0.23-0.195)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.23-kn; //выбор
	//kn=(0.29-0.195)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.29-kn;
	} //900
	else if (vyb == 10) { kktp[1] = 0.00025;  kktp[0] = 0.246; 
	kn=(0.287-0.25)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.287-kn; //выбор
	//kn=(0.36-0.25)/dt; kktp[1] = kn; kn=kn*t2; kktp[0] = 0.36-kn;
	//kn=(0.35-0.25)/dt; kktp[1]=kn; kn=kn*t2; kktp[0] = 0.35-kn; //к вопросу о стандартизации КВИ
	} //1000
	return kktp;
}
double *poisMasKoefsha(double *kektp, int n)
{
	double sfeosha=0.0, smgosha=0.0, salosha=0.0, ssiosha=0.0;
	double ko=1e-2, p20=2e1*ko, p24=24.0*ko, p30=3e1*ko, p33=33.0*ko, p16=16.0*ko, p10=1e1*ko; int k=0, vrsh=vybsha;
	if (!vrsh) { sfeosha = 1.24*ko; smgosha = 0.29*ko, salosha = 41.9*ko, ssiosha = 54.0*ko; porsha=21.8*ko; } //Роучка
	else if (vrsh==1) { sfeosha = 1.21*ko; smgosha = 0.29*ko; salosha = 42.6*ko; ssiosha = 1e0-salosha-sfeosha-smgosha; porsha=11.0144*ko; } //ШПД-41
	else if (vrsh==2) { sfeosha = 1.64*ko; smgosha = 0.36*ko; salosha = 35.9*ko; ssiosha = 59.1e-2;  porsha=25.2*ko; } //ШБ-1 2-1
	else if (vrsh==3) { sfeosha = 1.66*ko; smgosha = 0.4*ko; salosha = 37.3*ko; ssiosha = 57.4*ko;  porsha=26.5*ko; } //ШВ-1 1-1
	else if (vrsh==4) { sfeosha = 1.24*ko; smgosha = 0.29*ko; salosha = 41.9*ko; ssiosha = 54*ko;  porsha=11.5*ko; } //ШПД
	else if (vrsh==5) { sfeosha = 1.54*ko; smgosha = 0.3*ko; salosha = 38.6*ko; ssiosha = 56.5*ko;  porsha=16.5*ko; } //ШКУ-32 3-1
	for (k=0; k<n; k++) kektp[k]=0.0;
	if ((salosha >= 28e-2) && (salosha <= 38e-2)) {
		if ((porsha >= p20) && (porsha < p24)) { vybsha = 0; vystsha = 0; 
		kektp[3] = -0.435e-9; kektp[2] = 0.685e-6; kektp[1] = 0.134e-3; kektp[0] = 0.725; }
		else if ((porsha >= p24) && (porsha < p30)) { vybsha = 0; vystsha = 1; 
		kektp[3] = -0.867e-9; kektp[2] = 1.77e-6; kektp[1] = -0.523e-3; kektp[0] = 0.806; } //задание коэффициентов - шамот средней пористости
		else if ((porsha >= p16) && (porsha < p20)) { vybsha = 1; 
		kektp[3] = -0.397e-9; kektp[2] = 0.71e-6; kektp[1] = 0.011e-3; kektp[0] = 0.851; } //уплотненный шамот
		else if ((porsha >= p30) && (porsha <= p33)) { vybsha = 2; 
		kektp[3] = -0.377e-9; kektp[2] = 0.918e-6; kektp[1] = -0.338e-3; kektp[0] = 0.77; } //низкоплотный шамот
		else if ((porsha >= p10) && (porsha<p16)) { vybsha = 3; 
		kektp[3] = 0.0; kektp[2] = -0.607e-6; kektp[1] = 1.14e-3; kektp[0] = 0.641; }
	} //повышенной плотности шамот
	if ((salosha>38e-2) && (salosha <= 45e-2)) {
		if ((porsha >= p20) && (porsha < p24)) { vybsha = 0; vystsha = 0; 
		kektp[3] = -0.124e-9; kektp[2] = 0.215e-6; kektp[1] = 0.125e-3; kektp[0] = 1.01; }
		else if ((porsha >= p24) && (porsha < p30)) { vybsha = 0; vystsha = 1; 
		kektp[3] = -0.333e-9; kektp[2] = 0.805e-6; kektp[1] = -0.289e-3; kektp[0] = 0.903; } //задание коэффициентов - шамот средней пористости
		else if ((porsha >= p16) && (porsha < p20)) { vybsha = 1; 
		kektp[3] = 0.0; kektp[2] = -0.154e-6; kektp[1] = 0.369e-3; kektp[0] = 1.03; } //уплотненный шамот
		else if ((porsha >= p30) && (porsha < p33)) { vybsha = 2; 
		kektp[3] = -0.377e-9; kektp[2] = 0.918e-6; kektp[1] = -0.338e-3; kektp[0] = 0.77; }  //низкоплотный шамот
		else if ((porsha >= p10) && (porsha < p16)) { vybsha = 3; 
		kektp[3] = 0.0; kektp[2] = -0.141e-6; kektp[1] = 0.437e-3; kektp[0] = 1.32; }
	} //повышенной плотности шамот
	return kektp;
}