#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <limits> 
#include <iomanip> 
#include <math.h>
#include <fstream>
using namespace std;
const int L = 32;
const double dp = 0.025;
double p=0;
int N, p_filled;
bool table[L][L];
double random(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}
int main ()
{
	if (L==4){
		p=0.2;
		}
	else{
		p=0.45;
		}
	N=(1-p-0.2)/0.025-1;
	ofstream coordinates, pc;
	coordinates.open("coordinates");
	cout << "MODELING OF PERCOLLATION IN " << L << "x" << L << " system" << endl;
	//
	//RANDON NUMBER GENERATOR
    srand((unsigned int)time(0));
	for (int i=0; i<L; i++){
		for (int j=0; j<L; j++){
			table[i][j]=false;
			//start << i+0.5 << " " << j+0.5 << endl;
		}

	}
	for (int iter=0; iter<N; iter++){
		if (iter==0){
			coordinates << endl;
			coordinates << endl;	
		}
		coordinates << p << endl;
		cout << p << endl;
		for (int i=0; i<L; i++){
			for (int j=0; j<L; j++){
				double number;
				number = random(0,1);
				if (number<=p){
					table[i][j]=true;
					p_filled++;
					coordinates << i+0.5 << " " << j+0.5 << endl;
					}
				cout << table[i][j] << " ";
				}
			cout << endl;
		}
		coordinates << endl;
		coordinates << endl;
		//cout << endl;
		//cout << endl;
		cout << p << endl;
		p+=0.025;
		cout << "p = " << p << endl;
		cout << endl;
		cout << p_filled << p_filled;
		
	
	}	

}

