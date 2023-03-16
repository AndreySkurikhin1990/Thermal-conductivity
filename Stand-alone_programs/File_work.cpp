// File_work.cpp : Defines the entry point for the console application.
//
// ChngFile.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <conio.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string.h>
using namespace std;
#define N 16
int main()
{
	ifstream fi;
	long r=1, i=0, q=0, g=0, d=0, j=0, k=0;
	fi.open("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\Vermiculite_data_1.txt");
	FILE *fp, *fm;
	fm=fopen("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\verm00.txt","w");
	fp=fopen("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\verm01.txt","w");
	if ((!fp) || (!fi.is_open()) || (!fm)) { printf("Can not open file!"); getch(); exit(1); }
	char *buf, c, *st;
	buf=(char*)malloc(sizeof(char)*N);
	if (!buf) {cout << "No memory"; getch(); exit(1);}
	for (i=0; i<N; i++) buf[i]='0';
	while (!fi.eof())
	{
		fi.read(buf,N);
	for (i=0; i<N; i++) if (buf[i]=='\n') g++;
	r++;
	}
	double t=0, *a;
	a=(double*)malloc(sizeof(double)*2*g);
	st=(char*)malloc(sizeof(char)*N*r);	
	if ((!st) || (!a)) {cout << "No memory"; getch(); exit(1);}
	for (i=0; i<2*g; i++) a[i]=0.0;	
	for (i=0; i<N; i++) buf[i]='\0';
	for (i=0; i<r*N; i++) st[i]='\0';
	st[r*N-1]='\0';
	fi.close();
	fi.open("C:\\Users\\Андрей\\Documents\\_Аспирантура\\tmp\\itom_data_1.txt");
	if (!fi.is_open()) { printf("Can not open file!"); getch(); exit(1); }
			while (!fi.eof())
	{
		fi.read(buf,N);
		for (i=0; i<N; i++) if (buf[i]==',') buf[i]='.';
		for (i=0; i<N; i++) { st[j]=buf[i]; j++; }
		for (i=0; i<N; i++) buf[i]='\0';
	}
				fi.close();
	i=0;
			while (i<r*N)
	{
	for (j=0; j<N; j++)
		{ 
			c=st[i];
			if ((c=='\n') || (c=='\t') || (c==' '))
		{
	 buf[q]='\0';
					t=atof(buf); 
						if (d<2*g) a[d]=t; 
									d++; q=0;
					for (k=0; k<N; k++) buf[k]='\0'; 
					if ((q>N) || (i>r*N)) break;
					break;
			}
			buf[q]=c; q++; i++;
	if ((q>N) || (i>r*N)) break;
	}
	if ((st[i]=='\0') || (i>r*N)) break;
					i++;
		} 
	for (i=0; i<2*g; i++) if (!(i%2)) fprintf(fm,"%0.7lf\t",a[i]); else fprintf(fp,"%0.7lf\t",a[i]);
	free(buf);
	free(st);
	free(a);
	fclose(fm);
	fclose(fp);
	fi.close();
//getch();
return 0;
}