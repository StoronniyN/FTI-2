#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <limits> 
#include <iomanip> 
#include <math.h>
#include <fstream>
#include <stdio.h>
//#include <cstdlib>
using namespace std;
int mt, mtt, n=1000;
double u[1000], ut[1000], c, tim, sum, x;
const double t=0.00005, h= 0.01, nu=1;
int main() {
	ofstream test;
	test.open("test.dat");
	for (int i=0; i<n; i++){
		x=h*(i-1-n/2);
		u[i]=exp(-pow(x,2)/pow(0.5,2));
	}
	c=t*nu/pow(h,2);
	tim=0;
	mt=10;
	mtt=2000;
	for (int j=0; j<mt; j++){
		for (int jj=0; jj<mtt; jj++){
			ut[0]=u[0]+c*(u[1]-2*u[0]+u[1]);
			int i;
			for (i=1; i<n-1; i++){
				ut[i]=u[i]+c*(u[i-1]-2*u[i]+u[i+1]);
			}
			i=n;
			ut[n]=u[n]+c*(u[n-1]-2*u[n]+u[n-1]);
			tim+=t;
			for (int i=0; i<n; i++){
				u[i]=ut[i];
			}
		}
		sum=0;
		for (int i=0; i<n-1; i++){
			sum+=(h/2)*(u[i]+u[i+1]);
		}
		//test << tim << " " << sum << endl;
		cout << "meow  " << tim << " " << sum << endl;
		for (int i=0; i<n; i++){
			test << h*(i-1-n/2) << " " << u[i] << endl;
			//cout << h*(i-1-n/2) << " " << u[i] << endl;
		}
		test << endl;
		test << endl;
		cout << endl;
	}
		
	return 0;
}
