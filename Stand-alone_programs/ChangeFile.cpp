#include <fstream>
#include <iostream>
#include <conio.h>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <cstdlib>
using namespace std;
const int N = 1;
class words {
	public:
	char *str;
	words *pp; };
struct chisla {
	double zn;
	struct chisla *sl; };
void Razdeleniye();
void TestFunct();
void main()
{ 	int i=0;
	char *s=new char[N];
	//Razdeleniye(); //TestFunct(); i=getchar(); exit(1);
	ifstream fin;
	fin.open("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\tmp2.txt");
	if (!(fin.is_open())) { cout << "File is not found!" << endl; i=getchar(); exit(1); }
	ofstream fout;
	fout.open("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\2022\\tmp2_.txt");
	words *prev=NULL, *nex=NULL, *root=NULL;
	prev=new words;
	root=prev;
	while (!fin.eof()) {
	nex = new words;
	nex->str = new char[N];
	fin.read(nex->str,N);
	nex->str[N] = '\0';
	s = nex->str;
	for (i=0; i<N; i++)
	if (nex->str[i]=='.') 
	s[i]=',';
	nex->str = s;
	fout << nex->str;
	prev->pp = nex; }
	fin.close();
	fout.close(); }
void Razdeleniye()
{
int k=28, q=100, j, f=k, kk, no=7;
char *s=new char[q];
struct chisla *te, *pr, **mp=new struct chisla*[k], **mt, **mr=new struct chisla*[k];
if ((!s) || (!mr) || (!mp)) { cout << "No memory!" << endl; k=getchar(); exit(1); } 
j=0; while (j<k) { te=new struct chisla; 
if (!te) { cout << "No memory!" << endl; f=getchar(); exit(1); } 
mr[j]=te; mp[j]=te; j++; }
double t;
for (j=0; j<q; j++) s[j]='\0';
ifstream fin; fin.open("D:\\_Аспирантура\\tmp\\Doli_prop_Shamot-A.txt"); 
if (!fin.is_open()) { cout << "File is not open!" << endl; k=getchar(); exit(1); }
k=0; while (!fin.eof()) 
{ j=0; mt=new struct chisla*[f];
if (!mt) { cout << "No memory!" << endl; k=getchar(); exit(1); } 
while (j<f) {
fin.getline(s,q,'\n');
pr=mp[j];
t=atof(s);
te=new struct chisla;
if (!te) { cout << "No memory!" << endl; f=getchar(); exit(1); } 
mt[j]=te;
te->zn=t;
te->sl=NULL;
pr->sl=te;
j++; }
for (j=0; j<f; j++) mp[j]=mt[j];
k++; }
kk=k;
ofstream fout; fout.open("D:\\_Аспирантура\\tmp\\Doli_Prop_Shamot-1.txt", ios_base::out | ios_base::trunc); 
if (!fout.is_open()) { cout << "File is not open!" << endl; k = getchar(); exit(1); }
te=mr[no]; k=0; 
while ((k<=kk) && (te)) 
{ t=te->zn;
fout << t << endl;
te=te->sl;
k++; } 
fout.close();
j=0;
while (j<f) {
k=0; pr=mr[j]; te=pr;
while ((k<=kk) && (te)) {
te=pr->sl;
delete []pr;
k++;
pr=te;
} j++; }
fin.close(); }
void TestFunct()
{
	int n=10, n2=2*n, i, j, m=2;
char *s1=new char[n2], *s2=new char[n];
for (i=0; i<n2; i++) s1[i]='a'+i; //for (i=0; i<n; i++) s2[i]='\0';
for (i=m, j=0; (i<n2) && (j<n); i++, j++)
	{ cout << i << "\t" << j << "\t" << s1[i] << endl;
s2[j]=s1[i]; }
s2[j-1]='\0';
cout << s2 << endl;
delete []s1; delete []s2;
}