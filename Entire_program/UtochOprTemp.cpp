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
double **opredTempLPSten(double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *, double);
double **oprRasTemNach(int, double, double, double, double, double, int, double, double, double, char *);
double opredKTPTKToch(double *, double *, double, int);
double *opredTempTvKarFragm(int, double *, double *, int, double, double, double, double, double, char *, double, double *);
double opredTempStenFragmUchIzl(double, int, double *, double *, double *, double, int, int, double, double, double, double, double *, double, double, double, double *, double *, double);
double opredTempStenFragm(double, int, double *, double *, double *, double *, int, int, double, double, double, double);
double *opredTempStenFragmMasTem(int, double *, double *, double *, double *, int, int, double, double, double, double, double, char *, int, int);
//-------
double **oprRasTemNach(int kost, double htch, double hvozd, double hko, double a, double b, int w, double thol, double tgor, double tol, char *snm)
{
	int k=0, j=0, c=2; 
	double ht=hvozd+htch, tsre=0.0, e=1e-10, kx=a, bx=b, hkx=hko;
	double **mu=new double*[c], *xi=new double[kost], *teks=new double[kost];
	if ((!mu) || (!xi) || (!teks)) { cout << snm << endl; j=getchar(); exit(1); }
	if (!w) { kx=(thol-tgor)/tol; bx=tgor; }
	k=0; xi[k]=hkx+htch/2e0; for (k=1; k<kost; k++) xi[k]=xi[k-1]+ht; //массив середин каждой из стенок по толщине
	for (k=0; k<kost; k++) teks[k]=kx*xi[k]+bx;
	k=0; j=kost-1; tsre=(teks[k]+teks[j])/2e0; //qobsha = opredKTPTKToch(qob, etesha, tsre, ce); //for (k=0; k<kost; k++) cout << "k = " << k << "\txk = " << xi[k] << "\ttex = " << teks[k] << endl;
	mu[k]=teks; k++; mu[k]=xi;
	return mu;
}
double **opredTempLPSten(double Tnach, int n, double *ktpvo, double *vte, double *ete, double *ktptk, int ce, int dmv, double htk, 
	double hvo, double ptpob, int kost, char *snm, double zn)
{
	int q=2, m=0, k=0; 
	double tp=0.0, tt=tp, lam=tt, lamvo=tt, lamt=tt, ht=tt, **mu=new double*[q];
	double *Ts=new double[kost], *Tss=new double[kost], eps=1e-6; 
	if ((!Ts) || (!Tss) || (!mu)) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<kost; k++) { Ts[k]=0.0; Tss[k]=0.0; }
	k=0; tp=Tnach; tt=tp; Ts[k]=tt; q=1; m=0;
	for (k=1; k<n; k++)  { 
		lamt=opredKTPTKToch(ktptk, ete, tt, ce); if (lamt<0.0) lamt=0.0; //ищем КТП твердого скелета при данной температуре
		lamvo=opredKTPTKToch(ktpvo, vte, tt, dmv); if (lamvo<0.0) lamvo=0.0; //ищем КТП воздуха при данной температуре
		if (k%2) { ht=htk; lam=lamt; } else { ht=hvo; lam=lamvo; }
		if (fabs(lam)>eps) tt=tp+zn*ptpob*ht/lam; tp=tt; //рассчитываем температуру на границах стенок 
		if (k%2) { Tss[m]=tt; m++; } else { Ts[q]=tt; q++; } }
	k=0; mu[k]=Ts; k++; mu[k]=Tss; 
	return mu;
}
double *opredTempStenFragmMasTem(int ksu, double *ktpvo, double *vte, double *ete, double *ktptk, int ce, int dmv, double htk, 
	double hvo, double ptpob, double Tnach, double zn, char *snm, int kost, int vyb)
{
	int j=0, k=j+1;
	double *po=NULL; 
	double **mu=opredTempLPSten(Tnach, ksu, ktpvo, vte, ete, ktptk, ce, dmv, htk, hvo, ptpob, kost, snm, zn);
	double *Ts=mu[j], *Tss=mu[k];
	if (!vyb) { po=Ts; if (Tss) delete[]Tss; } else if (vyb==1) { po=Tss; if (Ts) delete[]Ts; } 
	else if (vyb==2) { if (Ts) delete[]Ts; if (Tss) delete[]Tss; po=new double[ksu]; if (!po) { cout << snm << endl; k=getchar(); exit(1); } 
	j=0; for (k=0; k<ksu; k++) if (!(k%2)) po[k]=Ts[j]; else { po[k]=Tss[j]; j++; } }
	if (mu) delete[]mu;
	return po;
}
double opredTempStenFragm(double Tna, int n, double *ktpvo, double *vte, double *ete, double *ktptk, int ce, int dmv, double htk, 
	double hvo, double ptpob, double zn)
{
	int k=0;
	double Tk=Tna, lamt=Tk, lamvo=Tk, lam=Tk, ht=Tk, Tn=Tk, Ts=Tk, Tkon=Ts, x=0.0;
	double qm=ptpob, hf=1e0, dt=hf, sig=5.67e-8, miep=1e-6;
	for (k=0; k<n; k++)  {
		Ts=(Tn+Tk)/2e0;
		lamt=opredKTPTKToch(ktptk, ete, Ts, ce); //ищем КТП твердого скелета при данной температуре
		lamvo=opredKTPTKToch(ktpvo, vte, Ts, dmv); //ищем КТП воздуха при данной температуре
		if (!(k%2)) { ht=htk; lam=lamt; } //твердый каркас
		else { ht=hvo; lam=lamvo; } //воздушная прослойка
		Tn=Tk; 
		Tk=Tk+zn*qm*ht/lam; } //рассчитываем температуру на границах стенок
	return fabs(Tna-Tk);
}
double opredTempStenFragmUchIzl(double Tna, int n, double *ktpvom, double *vtem, double *ete, double ktptk, int ce, int dmv, double htk, double hvo, double ptpob, double zn, double *stch, double ektp, double tolmss, double ksf, double *ktptkm, double *ektpm, double nenu)
{
	int k=0; 
	double Tk=Tna, lamt=Tk, lamvo=Tk, lam=Tk, ht=Tk, Tn=Tk, Ts=Tk, Tkon=Ts, x=0.0;
	double ql=0.0, qm=ql, eps=ql, epr=eps, hf=1e0, efktp=hf, dt=hf, sig=5.67e-8, miep=1e-6; //ektp=opredKTPTKToch(ektpm, ete, Tna, ce); 
	Tkon=fabs(Tna-fabs(ptpob*tolmss/ektp));
	Ts=(Tna+Tkon)/2e0;
	eps=opredKTPTKToch(stch, ete, Ts, ce); //ищем степень черноты при данной температуре
	if (fabs(eps)>miep) epr=hf/(2e0/eps-hf); else epr=0.0;
	ql=epr*sig*fabs(pow(Tkon, 4e0)-pow(Tna, 4e0))*nenu/ksf;
	qm=ptpob-ql; //cout << "\tqm = " << qm << "\tqo = " << ptpob << "\tql = " << ql << "\teps = " << eps << "\tepr = " << epr << "\ttk = " << tkon << "\tektp = " << efktp << "\ttn = " << Tna << "\ttol = " << tolmss << endl;
	for (k=0; k<(n-1); k++)  {
		Ts=(Tn+Tk)/2e0;
		lamt=opredKTPTKToch(ktptkm, ete, Ts, ce); //ищем КТП твердого скелета при данной температуре 
		lamvo=opredKTPTKToch(ktpvom, vtem, Ts, dmv); //ищем КТП воздуха при данной температуре
		if (!(k%2)) { ht=htk; lam=lamt; } //твердый каркас
		else { ht=hvo; lam=lamvo; } //воздушная прослойка
		Tn=Tk; 
		Tk=Tk+zn*qm*ht/lam; } //рассчитываем температуру на границах стенок
	return fabs(Tna-Tk);
}
double *opredTempTvKarFragm(int kost, double *ete, double *kttkve, int cemve, double htch, double hvoz, double qobtk, double temc, 
	double zn, char *snm, double Tnach, double *Tpctve)
{
	int k=0, n=kost; //cout << "k = " << k << "\t";
	double lamt=0.0, tx=Tnach; //cout << "n = " << n << "\t";
	double eps=1e-8, lam=eps, ht=htch+hvoz/2e0, ho=htch+hvoz; //cout << "ho = " << ho << "\tht = " << ht << "\tTn = " << Tnach << "\tqo_tk = " << qobtk << "\t";
	for (k=0; k<n; k++)  { 
		lam=opredKTPTKToch(kttkve, ete, tx, cemve); if (lam<0.0) lam=0.0; //ищем КТП твердого скелета при данной температуре
		if (fabs(lam)>eps) tx=Tnach+zn*qobtk*ht/lam; //cout << "tx ( " << k << " ) = " << tx << "\tlam = " << lam << "\t";
		Tpctve[k]=tx; ht=ht+ho; } //рассчитываем температуру на границах стенок твердогокаркаса //cout << endl; 
	return Tpctve;
}