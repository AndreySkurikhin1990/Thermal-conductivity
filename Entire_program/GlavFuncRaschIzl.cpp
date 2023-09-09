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
double opredKTPTKToch(double *, double *, double, int);
double LuchKTPChudnovsky(double *, double, int, double);
double **opredTempLPSten(double *, double *, double *, double *, double, int, double *, double *, double *, double *, int, int, double, double, double, int, char *);
double **RaschVneshIzluch(double *, double *, double *, double *, double *, int, double, double, char *);
double *reshnewtrafs(double *, double *, double *, double *, int, double *, double *, int, double *, double *, int, double *, double *, int, double, double *, double *, int, double, double, double, double *, int, double, double, double, char *);
double **RaschSobLuchPlotTepPot(int, double *, double *, double *, double *, double *, double *, double, double, double *, int, char *);
double **RasLuchPloTepPot(int, double *, double *, double *, double *, double *, double *, int, double *, double *, double *, char *, double);
double RasIzlSerStenNac(double *, double *, double *, double *, double *, double *, double *, double, double, double, double, int, double, double, double, double, int, double *, double *, int, int, double *, double *, double *, char *, int, double *, double *, int);
double *FuncRaschIzl(double, double, double, double, int, double *, double, double *, double, double, double, double, int, double *, double *, int, int, char *, double, int, double *, double, double *, double *, int, double *, double *, int, double *, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, int, double, double, int, double *, double *, double, double, double, char *, int, double);
void vyvodfile(double *, int, int, double, char *);
double *KoefPoglRosselNac(double *, int, int, double, double, double, double *, double *, int, double, double, double *, double *, int, int, int);
//-----
/*double *FuncRaschIzl(double ta, double tb, double d, double tc, int m, double *prot, double de, double *ao, double rez, double rap, double b, double c, int vybo, double *qobve, double *eteve, int cemve, int vyve, char *snms, double dpct, int vyte, double *ektpve, double y0ve, double *ktrve, double *kttkve, int kost, double *ktpvove, double *vteve, int dmkvove, double *ecktpve, double htch, double hvoz, double *Refe, double *Trae, double *Abse, double *Reft, double *Trat, double *Abst, double *txve, double *qxve, double *dtxve, int sctve, double tocrasve, double poristost, int dmao, double *stchsrve, double *Tpctve, double b0, double c0, double mk, char *navyfa, int fla, int ks, double hk)
{
	int j=0, zf=1, cels=6, k=0, nvm=11, ksu=2*kost, nfm=28;
	double qrr=0.0, qr=0.0, gtm = 0.0, dtma = 0.0, t = 0.0, hlr0 = 0.0;
	double hrl0 = 0.0, rhlr0 = 0.0, slr = 0.0, srl = 0.0, ra10, **mu, ra9, ra, ra7;
	double *Ts = NULL, *Tss = NULL, *sislr = NULL, *sisrl = NULL, *Tna = NULL, *reiz = NULL, *hrl1 = NULL, *hlr1 = NULL, *Tsr = NULL, qis;
	double *fm = new double[nfm], ktptk, ecktp, ktpr, dp = (1e0 - dpct), qs, *tevy = NULL, r, tesr, dote, p, *vm = new double[nvm];
	double qobvee, gg, ep = tocrasve, rr, ra0, tcc=(ta+tb)/2e0, *Evsv=NULL;
	if (!vm) { cout << snms << endl; k = getchar(); exit(1); }
	gtm = LuchKTPChudnovsky(Abse, tc, kost, hvoz); //по Чудновскому
	qobvee = opredKTPTKToch(qobve, eteve, tc, cemve); if (qobvee<=0.0) qobvee=opredKTPTKToch(qobve, eteve, eteve[vyte], cemve);
	qs = dp*qobvee;
	p = opredKTPTKToch(ektpve, eteve, tc, cemve);
	dtma = fabs(tb - ta);
	gtm = dtma / y0ve;
	ktpr = opredKTPTKToch(ktrve, eteve, tc, cemve);
	qr = ktpr*gtm; qrr = qr;
	ktptk = opredKTPTKToch(kttkve, eteve, tc, cemve);
	ecktp = opredKTPTKToch(ecktpve, eteve, tc, cemve);
	if (!m) { t = (dpct*ecktp*b + (1e0 - dp*b)*ktptk); if (fabs(t) > ep) hlr0 = qs*ktpr / t; else hlr0 = 0.0; }
	else hlr0 = rez; hrl0 = 0.0; rhlr0 = hlr0 - hrl0;
	if (qs<rhlr0) { rhlr0 = qs; hlr0 = rhlr0 + hrl0; }
	mu = opredTempLPSten(Ts, Tss, Tsr, Tna, tc, ksu, ktpvove, vteve, eteve, kttkve, cemve, dmkvove, htch, hvoz, qs - rhlr0, kost, snms);
	k = 0; Ts = mu[k]; k++; Tss = mu[k]; k++; Tsr = mu[k]; k++; Tna = mu[k]; if (mu) delete[]mu;
	mu = RaschVneshIzluch(Trae, Refe, prot, hlr1, hrl1, ks, hlr0, hrl0, snms);
	k = 0; hlr1 = mu[k]; k++; hrl1 = mu[k]; k++; reiz = mu[k]; slr = reiz[0]; srl = reiz[1]; ra0 = slr - srl; delete[]mu;
	tevy = reshnewtrafs(prot, Tna, hlr1, hrl1, zf, Trae, Abse, kost, kttkve, eteve, cemve, ktpvove, vteve, dmkvove, qs, txve, qxve, sctve, htch, hvoz, tocrasve, ektpve, cels + 2 * kost, tb, ta, y0ve, snms);
	k=2; tesr = tevy[k]; txve[sctve] = tesr; qxve[sctve] = tevy[1]; dtxve[sctve] = tevy[0]; gtm = tevy[4];
	if ((m>0) && (fabs(tevy[4]) > 0.0)) r = tevy[1] * b*dp / tevy[4] + opredKTPTKToch(kttkve, eteve, tesr, cemve)*(1e0 - b*dp); else r = 0.0; dote = tevy[5]; delete[]tevy;
	if (r<0.0) r = 0.0; if (fabs(p) > ep) gg = fabs((r - p) / p)*1e2; else gg = mk; if (mk < gg) { b0 = b; c0 = c; mk = gg; }
	for (j = 0; j < dmao; j++) ao[j] = 0.0;
	mu = RaschSobLuchPlotTepPot(kost, prot, Ts, Tss, Trae, Refe, Abse, slr, srl, ao, 0, snms); 
	k = 0; ao = mu[k]; k++; sislr = mu[k]; k++; sisrl = mu[k]; delete[]mu;
	mu = RasLuchPloTepPot(kost, hrl1, hlr1, ao, Trae, Refe, Abse, 2, sislr, sisrl, Tsr, snms); 
	k = 0; ao = mu[k]; k++; delete[]reiz; reiz = mu[k]; k++; Evsv=mu[k]; delete[]mu; 
	k = 0; slr = ao[k]; k++; srl = ao[k]; ra = slr - srl;
	qr = RasIzlSerStenNac(Refe, Trae, Abse, prot, Reft, Trat, Abst, tc, hlr0, hrl0, b, vyve, poristost, dpct, htch, hvoz, kost, eteve, kttkve, cemve, vyte, qobve, ktpvove, vteve, snms, dmao, stchsrve, Tpctve, dmkvove); qis = qr; //0 - шамот, 1 - вермикулит, 2 - ИТОМ
	qr = (qr*(d - 1e0) + ra)*dp*b / d;
	k = 7; ra7 = (ra + ao[k] * (d - 1e0)) / d; ra7 = ra7*dp*c*b; //ra7=qrr+ra7; 
	k = 9; ra9 = (ra + ao[k] * (d - 1e0)) / d; ra9 = ra9*dp*c*b; //ra9=ra9+qrr; //c - поправка на разрешающие угловые коэффициенты из-за небесконечности стенок, b - их аспектное соотношение (поправка на их некубичность)
	k = 10; ra10 = (ra + ao[k] * (d - 1e0)) / d; ra10 = ra10*dp*c*b; //ra10=qrr+ra10; 
	rez = ra; t = opredKTPTKToch(kttkve, eteve, tc, cemve); 
	if (fabs(p) > ep) gg = fabs((r - p) / p)*1e2; else gg = mk; if (mk>gg) { b0 = b; c0 = c; mk = gg; }
	if (fabs(qr) > ep) rr = ra9/qr; else rr = 0.0;
	k = 0; fm[k] = rr; k++; fm[k] = dote;   k++; fm[k] = qr;      k++; fm[k] = qrr;      k++; fm[k] = r;         k++; fm[k] = gg;    k++;
	fm[k] = b;         k++; fm[k] = c;      k++; fm[k] = qis;     k++; fm[k] = ao[5];    k++; fm[k] = ao[7];     k++; fm[k] = ao[9]; k++;
	fm[k] = ao[10];    k++; fm[k] = qobvee; k++; fm[k] = rap;     k++; fm[k] = t;        k++; fm[k] = ra9;       k++; fm[k] = ra10;  k++;
	fm[k] = ta;        k++; fm[k] = tb;     k++; fm[k] = ra0;	  k++; fm[k] = ra9 - qr; k++; fm[k] = ra10 - qr; k++; fm[k] = ra7;   k++; 
	fm[k] = ra7 - qr;  k++; fm[k] = ra;     k++; fm[k] = hk;      k++; fm[k] = tcc;	     k++;
	if (fla) { vyvodfile(fm, nfm, 2, 0.0, navyfa); //cout << "tc = " << tc << "\ttNR = " << dote << "\tqr = " << qr << "\tqrr = " << qrr; cout << "\tktp = " << r << "\topktp = " << gg << "\tb = " << b << "\tc = " << c << "\tqss = " << qis; cout << "\tao5 = " << ao[5] << "\tao7 = " << ao[7] << "\tao9 = " << ao[9] << "\tao10 = " << ao[10] << "\tqob = " << qobvee; cout << "\trap = " << rap << "\tktp_tk = " << t << "\tra9 = " << ra9 << "\tra10 = " << ra10; cout << "\tta = " << ta << "\ttb = " << tb << "\tslr-srl = " << slr - srl << "\tra9-qr = " << (ra9 - qr) << "\tra10-qr = " << ra10 - qr << endl; 
	} delete[]fm; delete[]sislr; delete[]sisrl; delete[]hlr1; delete[]hrl1; 
	delete[]reiz; delete[]Tsr; delete[]Ts; delete[]Tss; delete[]Tna; delete []Evsv;
	k = 0; vm[k] = ra; k++; vm[k] = qr; k++; vm[k] = qrr; k++; vm[k] = ra9 - qr; k++; vm[k] = dote; k++; 
	vm[k] = b0; k++; vm[k] = c0; k++; vm[k] = mk; k++; vm[k] = rr; k++; vm[k] = gg; k++; vm[k] = r; return vm;
}
double opredLuchSost(double *prot, double *Tsr, int le, int zf, int vyte, int N, int kost, double bn, double bk, double cn, double ck, char *snm, double qob, int vt, int vtn, int sctx, int nxt, int nnxt, double y0, double *dkosct, double *dkoscm, int dkoscl, double dkosp, double dkoal, double *tkusc, double *kusc, int dmkoosc, double *thol, double *tgor, double *ete, int cem, double smgo, double ssio, double salo)
{
	int k=0, j=0, q=0, m=0, l=3, mma=3000, mm=0, mt=0, vy=9, idr=1;
	double **mauk=NULL, *pt=NULL, s=0.0, p=0.0, rez=0.0, rap=0.0, *temras=NULL;
	double *ao=new double[2*N], *po=NULL, d=0.0, hb=5e-1, hb0=hb, hf=1e0, *ktpr=NULL;
	d=0.0; for (j=0; j<kost; j++) d=d+hf;
	double de=hf, ta=0.0, tb=ta, tc=ta, fa=ta, fb=tb, fc=tc;
	double ba=bn, bb=bk, bc=bn, xk=0.0, dxk=1e-3, co=0.0, bo=co;
	double ca=cn, cb=ck, cc=cb, hc=1e-1, hc0=hc, tt=0.0, raz=1e2;
	double mk=raz, mkm=mk, c00=0.0, b0=bc, c0=cc, kdc=5e0, kdb=5e0, tmin=0.0, tmax=0.0;
	double *R=NULL, *T=NULL, *A=NULL, *Ra=NULL, *Ta=NULL, *Aa=NULL, *Rt=NULL, *Tt=NULL, *At=NULL;
	if (!ao) { cout << snm << endl; k = getchar(); exit(1); }
	for (k=0; k<(2*N); k++) ao[k]=0.0; //cout << "qo = " << qob << "\n";
	if ((vt==vtn) && (sctx<=2) && (!temras)) {
		j=nxt-nnxt; temras=new double[j]; if (!temras) { cout << snm << endl; k = getchar(); exit(1); }
		ta=(thol[vt]-tgor[vt])/y0; tb=tgor[vt];
		xk=0.0; for (k=0; k<j; k++) { temras[k]=tb+ta*xk; xk=xk+dxk; }	}
	if ((vt==vtn) && (sctx<=2) && (!ktpr)) {
		ktpr=KoefPoglRosselNac(ete, idr, cem, smgo, ssio, salo, dkosct, dkoscm, dkoscl, dkosp, dkoal, tkusc, kusc, dmkoosc, idr, 1);
		for (k = 0; k < cem; k++) cout << "ktr = " << ktpr[k] << "\ttem = " << ete[k] << endl; } 
	l=vt; tmin=thol[l]; tmax=tgor[l];
	ta=tmin; tb=tmax; tc=temras[sctx];
	mauk=RaschRTA(kost, h, 0.0, 0.0, 1, vt, hvo, 0, tc, 0, 0, 1);
	k=0; Ra=mauk[k]; k++; Ta=mauk[k]; k++; Aa=mauk[k]; k++; R=mauk[k]; k++; T=mauk[k]; k++;
	A=mauk[k]; k++; Rt=mauk[k]; k++; Tt=mauk[k]; k++; At=mauk[k]; delete[]mauk;
	ca=cn; cb=ck; m=0; q=0;
	while ((m<mma) && (raz>1e-3) && (!q)) {
		cc=(ca+cb)/2e0;
		pt=FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cc, vy, qob, ete, cem, 0, snm, dpct, vt, ektp, y0, ktpr, kttk, kost, ktpvo, vte, dmkvo, ecktp, h, hvo, R, T, A, Rt, Tt, At, tx, qx, dtx, sctx, tocras, por, 2*N, stchsr, Tpct, b0, c0, mk, svfdv, 0, hk);
		k=3; fc=pt[k]; delete[]pt;
		if (fabs(fc)<hf) { c0=cc; break; }
		pt=FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, ca, vy, qob, ete, cem, 0, snm, dpct, vt, ektp, y0, ktpr, kttk, kost, ktpvo, vte, dmkvo, ecktp, h, hvo, R, T, A, Rt, Tt, At, tx, qx, dtx, sctx, tocras, por, 2*N, stchsr, Tpct, b0, c0, mk, svfdv, 0, hk);
		k=3; fa=pt[k]; delete[]pt;
		if (fabs(fa)<hf) { c0=ca; break; }
		pt=FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, cb, vy, qob, ete, cem, 0, snm, dpct, vt, ektp, y0, ktpr, kttk, kost, ktpvo, vte, dmkvo, ecktp, h, hvo, R, T, A, Rt, Tt, At, tx, qx, dtx, sctx, tocras, por, 2*N, stchsr, Tpct, b0, c0, mk, svfdv, 0, hk);
		k=3; fb=pt[k]; delete[]pt;
		if (fabs(fb)<hf) { c0 = cb; break; }
		if ((fa*fc)<0.0) { cb = cc; q = 0; }
		if ((fb*fc)<0.0) { ca = cc; q = 0; }
		if (((fb*fc)<0.0) && ((fa*fc)<0.0)) { cb = cc; q = 0; }
		if (((fb*fc)>0.0) && ((fa*fc)>0.0)) q = 1;
		m++; raz=fabs(cb-ca);
	}
	c0=cc; c00=c0; co=c00; cout << "c0 = " << c0 << "\ttb = " << tb << "\ttc = " << tc; 
	m=0; ba=bn; bb=bk; l=vt; ta=thol[l]; tb=tgor[l]; tc=temras[sctx];
	while ((m<mma) && ((bb-ba)>1e-1)) {
		bc=(ba+bb)/2e0; b0=bc; c0=c00; cc=c0; 
		pt=FuncRaschIzl(ta, tb, d, tc, m, prot, de, ao, rez, rap, bc, c0, vy, qob, ete, cem, 0, snm, dpct, vt, ektp, y0, ktpr, kttk, kost, ktpvo, vte, dmkvo, ecktp, h, hvo, R, T, A, Rt, Tt, At, tx, qx, dtx, sctx, tocras, por, 2*N, stchsr, Tpct, b0, c0, mk, svfdv, 0);
		k=9; s=pt[k]; k=7; mk=pt[k]; k=5; b0=pt[k]; k=4; tt=pt[k]; k=10; r=pt[k]; delete[]pt;
		if ((mk<mkm) && (mk>0.0)) { mkm=mk; bo=bc; b0=bc; lx[sctxv] = r; oox[sctxv] = s; tx[sctx] = tt; qx[sctxv] = qob; } //cout << "\tbc = " << bc << "\tmk = " << mk << "\tmkm = " << mkm << endl;
		m++; if (bc<hf) hb=hb0/kdb; else hb=hb0; bb=bb-hb;
	} cout << "\tbop = " << bo << "\tmop = " << mkm << "\tcop = " << co << endl;
	if (fabs(tnr-tt)>1e-1)
		tnr=tt; else {
			cn=co-hc0; if (cn<=0.0) cn = 1e-2;
			ck=co+2e0*hc0;
			bn=bo-hb0; if (bn<=0.0) bn=1e-2;
			bk=bo+2e0*hb0;
			tnr=tt;
		}
	delete[]ao; return rez;
}*/