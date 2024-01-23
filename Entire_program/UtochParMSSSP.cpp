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
const double pi = acos(-1e0);
double oprSrRazPor(int, int, double, int, int, int, int, int, int, char *);
double oprMaxRazPor(int, int, double, int, int, int, int, int, int, char *);
double **opredTolVozdPros(int, int, int, int, int, int, int, int, int, int, int, int, char *, int, double, double *, double *, 
	double *, double *, double *, double *, double *, double *, double *, int, int, double, double);
double opredKTPTKToch(double *, double *, double, int);
double opredTempStenFragmUchIzl(double, int, double *, double *, double *, double, int, int, double, double, double, double, double *, double, 
	double, double, double *, double *, double);
double *oprEffDoliTepPerenObFig(double, double, double, double, int, double, double, double *, double *, int, char *, double *, double, int, 
	double, double *, double, double *, double *, int, double, double, int, int, int, int);
double oprEffDoliTepPerenPodQobshObFig(double, double, double, double, int, double, double, double *, double *, int, char *, double *, 
	double, int, double, double *, double *, double *, double, int, double, double);
double oprEffDoliTepPerenPodTempObFig(double, double, double, double, int, double, double, double *, double *, int, char *, double *, 
	double, int, double, double *, double *, double *, double, int, double, double);
double oprEffDoliTepPerenPodEffKTPObFig(double, double, double, double, int, double, double, double *, double *, int, char *, double *, 
	double, int, double, double *, double *, double *, double);
double opredTolTverChasObFigObPor(double, double, double, double, double, int, double, double);
double oprOtnSechPeremMSS(double, double, double, double, double, double, double);
//-----
double **opredTolVozdPros(int vrsh, int vystsha, int kost, int vybves, int vmv, int vfv, int vsv, int vuv, int isrp, int vpmf, 
	int vmivmf, int vpkf, char *snm, int cem, double por, double *ktpvo, double *tevo, double *ete, double *ektp, double *tgor, 
	double *thol, double *qob, double *ktptkv, double *stch, int dmkvoz, int vogf, double p1, double p2)
{
	FILE *fo=NULL;
	int k=0, n=6, j=k+1, q=10000, qk=100, w=0, cel=0, u=0, cema=2, f1=0, f2=0;
	int f3=0, f4=0, p=0, y=0, yk=5, mo=-1, z=1, kon=10, pr=1, u1=1, u2=1, u3=1, u4=1;
	double e=1e-9, ep=1e-8, d=0.0, srp=d, ce=d, hf=1e0, nf=d, t=d, qobsh=d, ee=1e-10, gk=1e3;
	double mrp=0.0, tho=d, tgo=d, efktp=d, kttk=d, *rs=new double[cem], *po=NULL, *mx=new double[cem]; 
	double *htc=new double[cem], **mu=new double*[cema], *v=new double[cem], xa=0.0, xb=xa, ya=xa, yb=xa, yc=xa;
	double ra=0.0, g=2e0, x=hf/g, su=1e-3, hmi=1e-3, ko=1e-6;
	double faa=d, fbb=d, fcc=d, vz=d, *tvp=new double[cem];
	double sa=0.0, saa=sa, sb=0.0, sbb=sb, sc=0.0, scc=sc, tol=sc, plc=sc, plcc=sc, dtmsc=sc, ecktpc=sc; 
	double dttkc=0.0, *pro=new double[cem], htcv=0.0, htaa=0.0, htbb=0.0, htcc=0.0;	
	if ((!rs) || (!htc) || (!v) || (!pro) || (!mx)) { cout << snm << endl; k=getchar(); exit(1); } 
	for (k=0; k<cem; k++) { v[k]=0; rs[k]=0.0; htc[k]=0.0; mx[k]=0.0; pro[k]=0.0; }
	d=0.0; for (k=0; k<kost; k++) d=d+hf; 
	if (vybves<=1) { g=2e0; x=hf/g; }
	if (vybves>1) { g=4e0; x=hf/g; }
		srp=oprSrRazPor(vybves, vmv, por, vfv, vsv, isrp, vpkf, vrsh, vystsha, snm); 
		mrp=oprMaxRazPor(vybves, vmv, por, vfv, vsv, isrp, vpkf, vrsh, vystsha, snm); 
		if (srp>hmi) srp=srp*ko; if (mrp>hmi) mrp=mrp*ko; cout << "srp = " << srp << "\tmrp = " << mrp << "\t";
	p2=p1; xa=p1*srp; xb=p2*srp;
	for (k=0; k<cem; k++) { 
		t=ete[k]; cout << "t = " << t << "\t";
		qobsh=opredKTPTKToch(qob, ete, t, cem); 
		tho=opredKTPTKToch(thol, ete, t, cem); 
		tgo=opredKTPTKToch(tgor, ete, t, cem);
		efktp=opredKTPTKToch(ektp, ete, t, cem); 
		kttk=opredKTPTKToch(ktptkv, ete, t, cem);
	ra=hf; ya=x*srp; yb=g*srp; w=1; f1=1; f2=0; f3=0; f4=0;
	while ((ra>=e) && (w<=q)) {
		yc=(ya+yb)/2e0;
	po=oprEffDoliTepPerenObFig(ya, d, por, efktp, kost, qobsh, e, ktpvo, tevo, cem, snm, ete, kttk, dmkvoz, t, stch, mrp, 
		ktptkv, ektp, vogf, xa, xb, f1, f2, f3, f4); 
	u=0; sa=po[u]; u++; htaa=po[u]; u++; faa=po[u]; u++; v[k]=po[u];
	if (po) delete[]po; 
	po=oprEffDoliTepPerenObFig(yb, d, por, efktp, kost, qobsh, e, ktpvo, tevo, cem, snm, ete, kttk, dmkvoz, t, stch, mrp, 
		ktptkv, ektp, vogf, xa, xb, f1, f2, f3, f4); 
	u=0; sb=po[u]; u++; htbb=po[u]; u++; fbb=po[u]; u++; v[k]=po[u]; 
	if (po) delete[]po; 
	po=oprEffDoliTepPerenObFig(yc, d, por, efktp, kost, qobsh, e, ktpvo, tevo, cem, snm, ete, kttk, dmkvoz, t, stch, mrp, 
		ktptkv, ektp, vogf, xa, xb, f1, f2, f3, f4); 
	u=0; sc=po[u]; u++; htcc=po[u]; u++; fcc=po[u]; u++; v[k]=po[u];
	if (po) delete[]po; 
	//-----
	cout << "w = " << w << "\tfaa = " << faa << "\tfbb = " << fbb << "\tfcc = " << fcc << "\tya = " << ya << 
		"\tyb = " << yb << "\tyc = " << yc << "\tv = " << v[k] << "\tg = " << g << endl; 
	ra=fabs(ya-yb);
	if ((!(w%qk)) && (w<=(kon*qk))) { g=g*2e0; if (g>gk) break; 
	x=hf/g; ra=hf; ya=x*srp; yb=g*srp; w++; continue; }
	w++; 
	if ((fabs(faa)<ee) && (fabs(fbb)<ee) && (fabs(fcc)<ee)) 
	{ y++; if (y>=yk) { g=g*2e0; y=0; } if (g>gk) break; 
	x=hf/g; ra=hf; ya=x*srp; yb=g*srp; } } 
	htc[k]=htcc; //толщина перегородки
	mx[k]=sc*xa/(hf-2e0*sc); //толщина перемычки
	rs[k]=sc; //относительный размер перемычки
	tvp[k]=yc; //толщина воздушной прослойки
	if (pr) {
		if (fabs(v[k]-hf)<e) { 
			tol=xa*(d-hf)+htc[k]*d; 
			plc=(2e0*sc+xa)*(2e0*sc+xb);
			plcc=(plc-xa*xb)/plc; 
	dttkc=(qob[k]*plcc)*(plc*plcc);
	ecktpc=qob[k]*(hf-plcc)*(xa*xb); 
	fcc=dttkc+ecktpc-qob[k]*plc; 
	pro[k]=fcc; } } }
	for (k=0; k<cem; k++) cout << "k = " << k << "\thtc = " << htc[k] << "\trs = " << rs[k] << "\tv = " << v[k] << 
		"\tProverka = " << pro[k] << endl; 
	if (z) { 
		fo=fopen("Tolschina_peregorodki.txt", "at"); 
	if (!fo) { cout << "File is not open" << endl; k=getchar(); exit(1); }
	for (k=0; k<cem; k++) fprintf(fo, "T = %0.2lf\thtc = %0.20lf\tvyfv = %d\tvysv = %d\tvyuv = %d\tvmimf = %d\tvybmar = %d\tvybves = %d\tsrp = %0.10lf\tmrp = %0.10lf\n", 
		ete[k], htc[k], vfv, vsv, vuv, vmivmf, vmv, vybves, srp, mrp); fprintf(fo,"\n");
	fclose(fo); 
	fo=fopen("Otnositelnyi_razmer_peremychki.txt", "at"); 
	if (!fo) { cout << "File is not open" << endl; k=getchar(); exit(1); }
	for (k=0; k<cem; k++) fprintf(fo, "T = %0.2lf\tOPP = %0.20lf\tvyfv = %d\tvysv = %d\tvyuv = %d\tvmimf = %d\tvybmar = %d\tvybves = %d\n", 
		ete[k], rs[k], vfv, vsv, vuv, vmivmf, vmv, vybves); fprintf(fo,"\n");
		fclose(fo); } 
	if (v) delete[]v; if (pro) delete[]pro; if (mx) delete[]mx; if (tvp) delete[]tvp;
	k=0; mu[k]=rs; k++; mu[k]=htc;
	return mu;
}
double *oprEffDoliTepPerenObFig(double hvozd, double d, double por, double ektp, int kost, double qobsh, double epsilon, double *ktpvo, 
	double *vte, int cem, char *snm, double *ete, double kttk, int dmkvo, double tsr, double *stch, double marp, double *kttkm, 
	double *ektpm, int vogf, double pa, double pb, int f1, int f2, int f3, int f4)
{
	int jt=0, qt=0, kt=10000, kont=100, ce=4, k=0, f=0; 
	double hf=1e0, lavo=opredKTPTKToch(ktpvo, vte, tsr, dmkvo);
	double e=1e-9, htcha=0.0, htchb=0.0, htchc=0.0, emi=e/1e1, rat=hf; 
	double sa=0.0, sb=sa, sc=sa, fat=sa, fbt=sa, fct=sa, fdt=sa, pt=sa, ho=sa;
	double sc1=sa, sc2=sa, sc3=sa, sc4=sa, *vp=new double[ce], saa=sa, r=0.0;
	if (!vp) { cout << snm << endl; jt=getchar(); exit(1); } for (k=0; k<ce; k++) vp[k]=0.0;
	htcha=e; htchb=marp*(hf-por)/por;
	while ((qt<=kt) && (rat>e)) {
	htchc=(htcha+htchb)/2e0; 
	if (f1) sa=oprEffDoliTepPerenPodQobshObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htcha, vogf, pa, pb); //определение относительной толщины перемычки (ОТП) через равенство тепловых потоков
	if (f2) sa=oprEffDoliTepPerenPodEffKTPObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htcha); //определение ОТП через равенство коэффициентов теплопроводности
	if (f3) sa=oprEffDoliTepPerenPodTempObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htcha, vogf, pa, pb); //определение ОТП через равенство средних температур на конце МССП
	if (f4) sa=oprOtnSechPeremMSS(d, lavo, hvozd, qobsh, htcha, kttk, ektp); //определение ОТП в обычной МССП через чистую теплопроводность
	if (sa<0.0) sa=0.0; if (sa>hf) sa=hf;
	fat=opredTolTverChasObFigObPor(sa, d, hvozd, por, marp, vogf, pa, pb); //подтягивание толщины перемычки к пористости
	fat=htcha-fat;
	if (f1) sb=oprEffDoliTepPerenPodQobshObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htchb, vogf, pa, pb);
	if (f2) sb=oprEffDoliTepPerenPodEffKTPObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htchb);
	if (f3) sb=oprEffDoliTepPerenPodTempObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htchb, vogf, pa, pb);
	if (f4) sb=oprOtnSechPeremMSS(d, lavo, hvozd, qobsh, htchb, kttk, ektp);
	if (sb<0.0) sb=0.0;
	fbt=opredTolTverChasObFigObPor(sb, d, hvozd, por, marp, vogf, pa, pb);
	fbt=htchb-fbt; 
	if (f1) { sc=oprEffDoliTepPerenPodQobshObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htchc, vogf, pa, pb); if (sc<0.0) sc=0.0;
	sc1=sc; }
	if (f2) { sc=oprEffDoliTepPerenPodEffKTPObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htchc); if (sc<0.0) sc=0.0;
	sc2=sc; }
	if (f3) { sc=oprEffDoliTepPerenPodTempObFig(hvozd, d, por, ektp, kost, qobsh, epsilon, ktpvo, vte, cem, snm, ete, kttk, dmkvo, 
		tsr, stch, kttkm, ektpm, htchc, vogf, pa, pb); if (sc<0.0) sc=0.0;
	sc3=sc; }
	if (f4) { sc=oprOtnSechPeremMSS(d, lavo, hvozd, qobsh, htchc, kttk, ektp); if (sc<0.0) sc=0.0;
	sc4=sc; }
	fct=opredTolTverChasObFigObPor(sc, d, hvozd, por, marp, vogf, pa, pb);
	fct=htchc-fct;
	if ((sa<emi) || (sb<emi) || (sc<emi)) { fat=fbt; fct=fbt; qt++; continue; }
	if ((fct*fbt>0.0) && (fat*fct<0.0)) htchb=htchc;
	if ((fct*fat>0.0) && (fbt*fct<0.0)) htcha=htchc; 
	if (((fat<0.0) && (fbt<0.0) && (fct>0.0)) || ((fat>0.0) && (fbt>0.0) && (fct<0.0))) htchb=htchc;
	if (fbt*fat>0.0) { 
		if (jt>0) { f1=1; f2=0; f3=0; f4=0; }
		if (jt>kont) { f1=0; f2=1; f3=0; f4=0; } 
		if (jt>(2*kont)) { f1=0; f2=0; f3=1; f4=0; } 
		if (jt>(3*kont)) { f1=0; f2=0; f3=0; f4=1; } 
		if (jt>(4*kont)) { f1=0; f2=f1; f3=f1; f4=f1; 
		sc1=0.0; sc2=0.0; sc3=0.0;
		break; } jt++; } 
	rat=fabs(htcha-htchb); 
	qt++; }
	if ((sc1>e) && (f1)) { sc=sc1; f=f1; }
	if ((sc2>e) && (f2)) { sc=sc2; f=f2; }
	if ((sc3>e) && (f3)) { sc=sc3; f=f3; }
	if ((sc4>e) && (f4)) { sc=sc4; f=f4; }
	r=0.0; k=0; while (k<f) { k++; r=r+hf; }
	fct=opredTolTverChasObFigObPor(sc, d, hvozd, por, marp, vogf, pa, pb); //толщина перегородки //1
	fct=hvozd*(d-hf)+fct*d; //толщина стенки //sc1
	sa=pa*sc/(hf-2e0*sc); saa=pb*sc/(hf-2e0*sc);
	fat=(pa+2e0*sa)*(pb+2e0*saa)*fat; //полный объем МССП //sc1
	fat=pa*pb*hvozd*(d-hf)/fat-por; //проверка пористости //sc1
	k=0; vp[k]=sc; k++; vp[k]=htchc; k++; vp[k]=fat; k++; vp[k]=r;
	return vp;
}
double oprOtnSechPeremMSS(double d, double lavo, double hvozd, double qobsh, double htcha, double kttk, double ektp)
{
	double dt=d, dt2=d, a1=d, a2=d, tpv=lavo/hvozd, sc=d, ho=hvozd+htcha;
	dt=qobsh*ho/ektp; dt2=tpv*dt/(kttk/htcha+tpv); 
	a1=kttk*dt/ho; a2=kttk*dt2/htcha; sc=(qobsh-a2)/(a1-a2); 
	return sc;
}
double oprEffDoliTepPerenPodQobshObFig(double hvozd, double d, double por, double ektp, int kost, double qobsh, double epsilon, 
	double *ktpvo, double *vte, int cem, char *snm, double *ete, double kttk, int dmkvo, double tsr, double *stch, double *kttkm, 
	double *ektpm, double htch, int vogf, double xa, double xb)
{ //подстройка к общему тепловому потоку
	int k=0, j=0, f=1000, ksu=2*kost;
	double r=0.0, hf=1e0, bm=2e1, e=epsilon/bm, pla=r, plb=r, plc=r, plaa=r, plbb=r, plcc=r;
	double ecktp=r, dtms=r, p=r, saa=e, sba=r, sca=r, t=r, fa=r, fb=r, fc=r, mo=-hf;
	double dtmsc=r, dtmsa=r, dtmsb=r, dttka=r, dttkb=r, dttkc=r, ecktpa=r, ecktpb=r, ecktpc=r;
	double sab=r, sbb=r, scb=r, ra=r, rb=r;
	p=hvozd*(d-hf)+htch*d; 
	saa=e*xa; sba=(bm*hf-e)*xa; //sab=e*xb; sbb=(bm*hf-e)*xb; 
	k=0; ra=hf; rb=hf;
	while ((ra>epsilon) && (k<=f) && (rb>epsilon)) {
		sca=(sba+saa)/2e0; scb=(sab+sbb)/2e0;
			pla=(2e0*saa+xa)*(2e0*sab+xb); //площадь сечения МССП поперек направления ПТП
			plb=(2e0*sba+xa)*(2e0*sbb+xb);
			plc=(2e0*sca+xa)*(2e0*scb+xb);
		plaa=(pla-xa*xb)/pla; //относительная площадь ТК
		plbb=(plb-xa*xb)/plb;
		plcc=(plc-xa*xb)/plc;
	dttka=(qobsh*plaa)*(plaa*pla); //тепловой поток через перемычку твердого каркаса
	dttkb=(qobsh*plbb)*(plbb*plb);
	dttkc=(qobsh*plcc)*(plcc*plc);
		ecktpa=(qobsh*(hf-plaa))*xa*xb; //тепловые потоки через МССП
		ecktpb=(qobsh*(hf-plbb))*xa*xb;
		ecktpc=(qobsh*(hf-plcc))*xa*xb;
				fa=dttka+ecktpa-qobsh*pla; //уравновешиваем тепловые потоки
				fb=dttkb+ecktpb-qobsh*plb;
				fc=dttkc+ecktpc-qobsh*plc;
			if (vogf==1) { r=pi/4e0; fa=fa*r; fb=fb*r; fc=fc*r; }
			if ((fc*fb>0.0) && (fa*fc<0.0)) { sba=sca; sbb=scb; }
			if ((fc*fa>0.0) && (fb*fc<0.0)) { saa=sca; sab=scb; }
			ra=fabs(saa-sba); 
			rb=fabs(sab-sbb); 
			k++; } 
	sca=sca/(xa+2e0*sca);
	scb=scb/(xb+2e0*scb); //sca=(sca+scb)/2e0;
	if (k<f) return sca; else return 0.0;
}
double oprEffDoliTepPerenPodTempObFig(double hvozd, double d, double por, double ektp, int kost, double qobsh, double epsilon, 
	double *ktpvo, double *vte, int cem, char *snm, double *ete, double kttk, int dmkvo, double tsr, double *stch, double *kttkm, 
	double *ektpm, double htch, int vogf, double xa, double xb)
{ //подстройка к конечной температуре стенки
	int k=0, j=0, f=1000, ksu=2*kost; 
	double r=0.0, hf=1e0, tef=r, kctp=r, tctp=r, tms=r, bm=2e1, e=epsilon/bm;
	double ecktp=r, dt=r, p=r, saa=r, sba=r, sca=r, t=r, fa=r, fb=r, fc=r, mo=-hf, zn=mo;
	double sab=r, sbb=r, scb=r, ra=r, rb=r, pla=r, plb=r, plc=r, plaa=r, plbb=r, plcc=r;
	p=hvozd*(d-hf)+htch*d; 
	tef=tsr+zn*fabs(qobsh*p/ektp); //разность температур
	r=(tsr+tef)/2e0;
	kctp=opredKTPTKToch(kttkm, ete, r, cem);
	saa=e*xa; sba=(bm*hf-e)*xa; 
	sab=e*xb; sbb=(bm*hf-e)*xb; 
	k=0; ra=hf; rb=hf; epsilon=1e-6;
	while ((ra>epsilon) && (k<=f) && (rb>epsilon)) {
			sca=(saa+sba)/2e0;
			scb=(sab+sbb)/2e0;
				pla=(2e0*saa+xa)*(2e0*sab+xb);
				plb=(2e0*sba+xa)*(2e0*sbb+xb);
				plc=(2e0*sca+xa)*(2e0*scb+xb);
		plaa=(pla-xa*xb)/pla; //Относительная площадь перемычки для переноса тепла
		plbb=(plb-xa*xb)/plb;
		plcc=(plc-xa*xb)/plc;
			tctp=fabs(qobsh*plaa*p/kctp);
			tctp=tsr+zn*tctp; //температура на конце перемычки через всю многослойную стенку
			tms=opredTempStenFragmUchIzl(tsr, ksu, ktpvo, vte, ete, kttk, cem, dmkvo, htch, hvozd, qobsh*(hf-plaa), mo, stch, ektp, p, d, 
				kttkm, ektpm, hf);
			tms=tsr+zn*fabs(tms); //температура на конце при переносе тепла через многослойную стенку
			fa=(pla*plaa)*tctp+tms*(xa*xb);
			fa=fa/pla-tef; //средние температуры по поверхности должны сравняться
		tctp=fabs(qobsh*plbb*p/kctp);
		tctp=tsr+zn*tctp;
			tms=opredTempStenFragmUchIzl(tsr, ksu, ktpvo, vte, ete, kttk, cem, dmkvo, htch, hvozd, qobsh*(hf-plbb), mo, stch, ektp, p, d, 
				kttkm, ektpm, hf);//максимальная разность температур - на краях стенок //в середине слоя 
			tms=tsr+zn*fabs(tms);
			fb=(plb*plbb)*tctp+tms*(xa*xb);
			fb=fb/plb-tef;
		tctp=fabs(qobsh*plcc*p/kctp);
		tctp=tsr+zn*tctp;
			tms=opredTempStenFragmUchIzl(tsr, ksu, ktpvo, vte, ete, kttk, cem, dmkvo, htch, hvozd, qobsh*(hf-plcc), mo, stch, ektp, p, d, 
				kttkm, ektpm, hf);//максимальная разность температур - на краях стенок //в середине слоя 
			tms=tsr+zn*fabs(tms);
			fc=(plc*plcc)*tctp+tms*(xa*xb);
			fc=fc/plc-tef;
			if (vogf==1) { r=pi/4e0; fa=r*fa; fb=r*fb; fc=r*fc; }
			if ((fc*fb>0.0) && (fa*fc<0.0)) { sba=sca; sbb=scb; }
			if ((fc*fa>0.0) && (fb*fc<0.0)) { saa=sca; sab=scb; }
			ra=fabs(saa-sba); 
			ra=fabs(sab-sbb); 
			k++; } 
	sca=sca/(xa+2e0*sca);
	scb=scb/(xb+2e0*scb); //sca=(sca+scb)/2e0;
	if (k<f) return sca; else return 0.0;
}
double oprEffDoliTepPerenPodEffKTPObFig(double hvozd, double d, double por, double ektp, int kost, double qobsh, double epsilon, 
	double *ktpvo, double *vte, int cem, char *snm, double *ete, double kttk, int dmkvo, double tsr, double *stch, double *kttkm, 
	double *ektpm, double htch)
{ //подстройка к эффективному КТП многослойной стенки
	int k=0, j=0, f=1000, ksu=2*kost; 
	double r=0.0, hf=1e0, bm=2e1, e=epsilon/bm;
	double p=r, sa=r, sb=r, sc=r, t=r, fa=r, fb=r, fc=r, mo=-hf;
	double dta=r, dtb=r, dtc=r, ecktpa=r, ecktpb=r, ecktpc=r;
	p=hvozd*(d-hf)+htch*d; 
	sa=e; sb=hf/2e0-e; 
	k=0; r=hf;
	while ((r>epsilon) && (k<=f)) {
			sc=(sa+sb)/2e0;
			dta=opredTempStenFragmUchIzl(tsr, ksu, ktpvo, vte, ete, kttk, cem, dmkvo, htch, hvozd, qobsh*(hf-2e0*sa), mo, stch, 
		ektp, p, d, kttkm, ektpm, hf)/p; //температурный градиент
	ecktpa=qobsh*(hf-2e0*sa)/dta;
	fa=kttk*(2e0*sa)+ecktpa*(hf-2e0*sa)-ektp; //КТП многослойной стенки с перемычкой должны сравняться с ЭКТП
			dtb=opredTempStenFragmUchIzl(tsr, ksu, ktpvo, vte, ete, kttk, cem, dmkvo, htch, hvozd, qobsh*(hf-2e0*sb), mo, stch, 
		ektp, p, d, kttkm, ektpm, hf)/p;
	ecktpb=qobsh*(hf-2e0*sb)/dtb; 
	fb=kttk*(2e0*sb)+ecktpb*(hf-2e0*sb)-ektp; 
	dtc=opredTempStenFragmUchIzl(tsr, ksu, ktpvo, vte, ete, kttk, cem, dmkvo, htch, hvozd, qobsh*(hf-2e0*sc), mo, stch, 
		ektp, p, d, kttkm, ektpm, hf)/p; 
	ecktpc=qobsh*(hf-2e0*sc)/dtc; 
	fc=kttk*(2e0*sc)+ecktpc*(hf-2e0*sc)-ektp;
		if ((fc*fb>0.0) && (fa*fc<0.0)) sb=sc;
		if ((fc*fa>0.0) && (fb*fc<0.0)) sa=sc;
			r=fabs(sa-sb); k++;	} 
	if (k<f) return sc; else return 0.0;
}
double opredTolTverChasObFigObPor(double sc, double d, double hvo, double por, double marp, int vogf, double xa, double xb)
{
	double hf=1e0, e=1e-9, hta=e, htb=marp, htc=hf, fa=0.0, fb=fa, fc=fa, ra=hf, hmi=1e-3;
	double pta=hf, ptb=hf, ptc=hf, sa=xa*sc/(hf-2e0*sc), sb=xb*sc/(hf-2e0*sc), r=hf, ko=1e-6;
	int j=0, k=1000, q=0, kon=k/10;
	if (hvo>hmi) hvo=hvo*ko; if (marp>hmi) marp=marp*ko;
	while ((j<=k) && (ra>e)) {
	htc=(hta+htb)/2e0;
	pta=hvo*(d-hf)+hta*d; //толщина стенки
	pta=(xa+2e0*sa)*(xb+2e0*sb)*pta; //полный объем
	fa=xa*xb*hvo*(d-hf)/pta; //объем пор
	fa=fa-por;
	ptb=hvo*(d-hf)+htb*d;
	ptb=(xa+2e0*sa)*(xb+2e0*sb)*ptb; //полный объем
	fb=xa*xb*hvo*(d-hf)/ptb; //объем пор
	fb=fb-por;
	ptc=hvo*(d-hf)+htc*d;
	ptc=(xa+2e0*sa)*(xb+2e0*sb)*ptc; //полный объем
	fc=xa*xb*hvo*(d-hf)/ptc; //объем пор
	fc=fc-por;
	if (vogf==1) { r=pi/4e0; fa=fa*r; fb=fb*r; fc=fc*r; }
	if ((fc*fb>0.0) && (fa*fc<0.0)) htb=htc;
	if ((fc*fa>0.0) && (fb*fc<0.0)) hta=htc; 
	if ((fc*fa>0.0) && (fb*fc>0.0)) { q++; if (q>kon) break; }
	ra=fabs(htb-hta); j++; }
	if (j<k) return htc; else return 0.0;
}