#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <limits> 
#include <iomanip> 
#include <math.h>
#include <fstream>
using namespace std;
const int N=5, a=5, number_of_chains=pow(2,N);
int chain[N][number_of_chains],dir,x[N][number_of_chains], y[N][number_of_chains], x_prev, y_prev;
int main(){
	ofstream coord;
	coord.open("coordinates");
		for (int j=0; j<number_of_chains; j++){
		x_prev=a*j;
		y_prev=a*j;
		x[0][j]=x_prev;
		y[0][j]=y_prev;
		cout << x[0][j] << " " << y[0][j] << endl;
		coord << x[0][j] << " " << y[0][j] << endl; 
		srand (time(0));
		for (int i=1; i<N+1; i++){
			dir=rand()%2;
			dir = (dir==0) ? -1: 1;
			if (i%2==0){
				x[i][j]=x_prev+dir;
				y[i][j]=y_prev;
			}
			else {
				y[i][j]=y_prev+dir;
				x[i][j]=x_prev;
			}
			for (int q=0; q<N+1; q++){
				if (q==i){
					continue;
				}
				if (x[i][j]==x[q][j] and y[i][j]==y[q][j]){
					dir=-dir;
				}
				if (i%2==0){
					x[i][j]=x_prev+dir;
					y[i][j]=y_prev;
				}
				else {
					y[i][j]=y_prev+dir;
					x[i][j]=x_prev;
				}
			}

			//cout << "x_prev = "<< x_prev << ", y_prev = " << y_prev << endl;
			//cout << "dir = " << dir << endl;
			cout << x[i][j] << " " << y[i][j] << endl;
			coord << x[i][j] << " " << y[i][j] << endl;
			//cout << endl;
			x_prev=x[i][j];
			y_prev=y[i][j];
		}
	
	coord << endl;
	}
	
	
	
	return 0;
	}
