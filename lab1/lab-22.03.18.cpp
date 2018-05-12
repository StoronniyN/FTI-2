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

const int number_of_particles = 36, u=sqrt(number_of_particles); 
const int step = 1000; //skip 1000 coordinates to write in an output file
const int NTIME=1000000 /* number of steps */, control_temp=1 /* control of temperature enable/disable 1/0 */;
int o = 0; //control temperature steps
double length[number_of_particles][number_of_particles];
const double sigma=2.74, epsilon=0.0031, lattice_constant=40, dt=0.0001, primary_velocity = 1;
double temperature_average, temperature_desired, temperature_prev, temperature, temperature_sum, kin_energy[NTIME], pot_energy[NTIME],  podgonian = 0;
int stop = 0;
double R, Z, dx, dy;
const int control_point_1=10/dt, control_point_2=50/dt, control_point_3=100/dt; //steps at which temperature control starts
double random(double min, double max)
{
    return (double)(rand())/RAND_MAX*(max - min) + min;
}
//STRUCTURE
struct primary_pairs
{
	public:
		double x;
		double y;
};
 	primary_pairs coord[number_of_particles][NTIME], acceleration[number_of_particles][NTIME], velocity[number_of_particles][NTIME];
	primary_pairs delta[number_of_particles][number_of_particles], velocity_average, sum_of_velocities;
int main() 
{
ofstream pulse, write, kinwrite, potwrite, fullwrite, tempwrite, avtempwrite, partwrite, rwrite, zwrite;
write.open("data.txt");
kinwrite.open("kinenergy");
potwrite.open("potenergy");
fullwrite.open("fullenergy");
tempwrite.open("temperature");
avtempwrite.open("avtemperature");
partwrite.open("particle");
rwrite.open("R");
zwrite.open("Z");
pulse.open("pulse");
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			length[i][j]=0;
		}
	}	
	//RANDOM NUMBER GENERATOR
    srand((unsigned int)time(0));

	for (int i=0; i<number_of_particles; i++)
	{
		for (int N=0; N<NTIME; N++)
		{
			coord[i][N].x = 0;
			coord[i][N].y = 0;
			acceleration[i][N].x = 0;
			acceleration[i][N].y = 0;
			velocity[i][N].x = 0;
			velocity[i][N].y = 0;
			kin_energy[N] = 0;
			pot_energy[N] = 0;
		}
		for (int j=0; j<number_of_particles; j++)
		{
			delta[i][j].x=0;
			delta[i][j].y=0;
		}
	}
    for (int i=0; i<number_of_particles; i++)
    {
		velocity[i][0].x =random(-primary_velocity, primary_velocity);
		//////cout << "v_x_" << i << "= " << velocity[i].x << endl;
		velocity[i][0].y = random(-primary_velocity, primary_velocity);
		//////cout << "v_y_" <<i << "= " << velocity[i].y << endl;
	}	
//CONDITION SUM[v_i] - number_of_particles*<V> = 0, <V>=SUM[v_i]

	velocity_average.x = 0;
	velocity_average.y = 0;
	sum_of_velocities.x = 0;
	sum_of_velocities.y = 0;
	for (int i=0; i<number_of_particles; i++)
	{ 
		//calculate a sum of average speed
		sum_of_velocities.x += velocity[i][0].x;
		sum_of_velocities.y += velocity[i][0].y;
	}
	//calculate the average speed
	velocity_average.x = sum_of_velocities.x/number_of_particles;
	velocity_average.y = sum_of_velocities.y/number_of_particles;
	for (int i=0; i<number_of_particles; i++)
	{
		velocity[i][0].x = velocity[i][0].x - velocity_average.x;
		velocity[i][0].y = velocity[i][0].y - velocity_average.y;
	}

double k=0;
for (int i=0; i<number_of_particles; i++)
{
	k+=velocity[i][0].x+velocity[i][0].y;
}
k=0;
pulse << 0 << " " << k <<endl;
	
// KIN ENERGIES 0
	for (int i=0; i<number_of_particles; i++)
	{
		kin_energy[0]+=(pow(velocity[i][0].x,2)+pow(velocity[i][0].y,2))/2;
	}	
//COORDINATES 0
	for (int a=0; a<u;  a++)
	{
		for (int i=0; i<u; i++)
		{
			coord[a*u+i][0].x=-lattice_constant/2+0.1*lattice_constant+i*(0.8*lattice_constant)/(sqrt(number_of_particles)-1);//random(-lattice_constant/2, lattice_constant/2);
			coord[a*u+i][0].y=-lattice_constant/2+0.1*lattice_constant+a*(0.8*lattice_constant)/(sqrt(number_of_particles)-1);//random(-lattice_constant/2, lattice_constant/2);
			//cout <<  coord[i][0].x << " " << coord[i][0].y << endl;
			write <<  coord[i][0].x << " " << coord[i][0].y << endl;
		}
	}
	//cout << "" << endl;
	//cout << "" << endl;
	write << "" << endl;
	write << "" << endl;
	for (int i=0; i<number_of_particles; i++) //CHOOSE A PARTICLE TO FIND A DISTANCE
	{
			for ( int j=0; j<number_of_particles; j++)
			{
				
				if (i==j)
				{
					continue;
				}
				delta[i][j].x=(coord[i][0].x-coord[j][0].x);
				delta[i][j].y=(coord[i][0].y-coord[j][0].y);
				length[i][j]=sqrt(abs(pow(delta[i][j].x,2)+pow(delta[i][j].y,2)));
				if (abs(delta[i][j].x)>lattice_constant/2 )
					{
						delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-lattice_constant);
					}
				if (abs(delta[i][j].y)>lattice_constant/2 )
					{
						delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-lattice_constant);
					}
			}
	}
	for (int i=0; i<number_of_particles; i++) //CHOOSE A PARTICLE TO FIND A DISTANCE
	{
		for ( int j=i+1; j<number_of_particles; j++)
		{
			pot_energy[0]+=4*epsilon*( pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
		}
	}
	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=0; j<number_of_particles; j++)
		{
			if (i==j)
			{
					continue;
			}
			
			acceleration[i][0].x+=-epsilon*(24*delta[i][j].x/pow(length[i][j],2))*( ( (pow(sigma/length[i][j],6) )) - ( (2*pow(sigma/length[i][j],12))) );
			acceleration[i][0].y+=-epsilon*(24*delta[i][j].y/pow(length[i][j],2))*( ( (pow(sigma/length[i][j],6) )) - ( (2*pow(sigma/length[i][j],12))) );
		}
			
	}
	for (int i=0; i<number_of_particles; i++)
	{
		coord[i][1].x=coord[i][0].x+velocity[i][0].x*dt+pow(dt,2)/2*acceleration[i][0].x;
		coord[i][1].y=coord[i][0].y+velocity[i][0].y*dt+pow(dt,2)/2*acceleration[i][0].y;
		if (coord[i][1].x>lattice_constant/2)
		{
			coord[i][1].x=coord[i][1].x-lattice_constant;
		}
		else if (coord[i][1].x<-lattice_constant/2)
		{
			coord[i][1].x=coord[i][1].x+lattice_constant;
		}
		if (coord[i][1].y>lattice_constant/2)
		{
			coord[i][1].y=coord[i][1].y-lattice_constant;
		}
		else if (coord[i][1].y<-lattice_constant/2)
		{
			coord[i][1].y=coord[i][1].y+lattice_constant;
		}		
		
		//////////cout << "i: " << coord[i][N+1].x << " , " << coord[i][N+1].y << endl;
		//cout << coord[i][1].x << " " << coord[i][1].y << endl;
		//cout << coord[i][N].x << endl;
		//cout << coord[i][N].y << endl;
	}
		//cout << "" << endl;
		//cout << "" << endl;
	kinwrite << 0 << " " << kin_energy[0] << endl;
	potwrite <<  0 << " " <<  pot_energy[0] << endl; 
	fullwrite << 0  << " " << kin_energy[0]+pot_energy[0] << endl; 
//CYCLE

for (int N=1; N<NTIME; N++)
{
	if (stop==1)
	{
		cout << "ЗАВЕРШЕНИЕ" << endl;
		cout << N-1 << endl;
		break;
	}
	for (int i=0; i<number_of_particles; i++) //CHOOSE A PARTICLE TO FIND A DISTANCE
	{
			for ( int j=0; j<number_of_particles; j++)
			{
				
				if (i==j)
				{
					continue;
				}
				delta[i][j].x=(coord[i][N].x-coord[j][N].x);
				delta[i][j].y=(coord[i][N].y-coord[j][N].y);
				length[i][j]=sqrt(pow(delta[i][j].x,2)+pow(delta[i][j].y,2));
				
				//cout << "dx["<< i << "," << j<< "] = " << delta[i][j].x << endl;
				//cout << "dy["<< i << "," << j<< "] = " << delta[i][j].y << endl;
				//cout << "length["<< i << "," << j<< "] = " << length[i][j] << endl;
				if (abs(delta[i][j].x)>lattice_constant/2 )
					{
						delta[i][j].x=((delta[i][j].x)/(abs(delta[i][j].x)))*(abs(delta[i][j].x)-lattice_constant);
					}
				if (abs(delta[i][j].y)>lattice_constant/2 )
					{
						delta[i][j].y=((delta[i][j].y)/(abs(delta[i][j].y)))*(abs(delta[i][j].y)-lattice_constant);
					}
				acceleration[i][N].x+=-(24*delta[i][j].x/pow(length[i][j],2))*epsilon*( ( (pow(sigma/length[i][j],6) )) - ( (2*pow(sigma/length[i][j],12))) );
				acceleration[i][N].y+=-(24*delta[i][j].y/pow(length[i][j],2))*epsilon*( ( (pow(sigma/length[i][j],6) )) - ( (2*pow(sigma/length[i][j],12))) );

			}
	}

	for (int i=0; i<number_of_particles; i++)
	{
		for (int j=i+1; j<number_of_particles; j++)
		{
			pot_energy[N]+=4*epsilon*(pow((sigma/length[i][j]),12) - pow((sigma/length[i][j]),6) );
		}
	}
// - VELOCITIES
	for (int i=0; i<number_of_particles; i++)
	{
		velocity[i][N].x=velocity[i][N-1].x+dt/2*(acceleration[i][N-1].x + acceleration[i][N].x);
		velocity[i][N].y=velocity[i][N-1].y+dt/2*(acceleration[i][N-1].y + acceleration[i][N].y);;
		//cout << "velocity_x_" << i << "_" << N << "= " <<  velocity[i][N].x << endl;			
		//cout << "velocity_y_" << i << "_" << N << "= " <<  velocity[i][N].y << endl;
	}
// KIN ENERGY AND TEMPETATURE
	for (int i=0; i<number_of_particles; i++)
	{
		kin_energy[N]+=(pow(velocity[i][N].x,2)+pow(velocity[i][N].y,2))/2;
	}	
	temperature_sum+=kin_energy[N]/number_of_particles;
	temperature_average=temperature_sum/N;
	temperature=kin_energy[N]/number_of_particles;
	temperature_desired=0.4;
	double coeff;
	coeff=temperature_desired/temperature_average;
	if (control_temp == 1 and (N==control_point_1 or N==control_point_2 or N==control_point_3))//and temperature_average-temperature_prev<0.00001)
	{
		for (int i=0; i<number_of_particles; i++)
		{
			velocity[i][N].x=sqrt(coeff)*velocity[i][N].x;
			velocity[i][N].y=sqrt(coeff)*velocity[i][N].y;
			//temperature_sum=0;
			//failback=1;
			//break;
		}
		o++;
	}
	
	for (int i=0; i<number_of_particles; i++)
	{
		coord[i][N+1].x=coord[i][N].x+velocity[i][N].x*dt+(1/2)*pow(dt,2)*acceleration[i][N].x;
		coord[i][N+1].y=coord[i][N].y+velocity[i][N].y*dt+(1/2)*pow(dt,2)*acceleration[i][N].y;
		if (coord[i][N+1].x>lattice_constant/2)
		{
			coord[i][N+1].x=coord[i][N+1].x-lattice_constant;
		}
		else if (coord[i][N+1].x<-lattice_constant/2)
		{
			coord[i][N+1].x=coord[i][N+1].x+lattice_constant;
		}
		if (coord[i][N+1].y>lattice_constant/2)
		{
			coord[i][N+1].y=coord[i][N+1].y-lattice_constant;
		}
		else if (coord[i][N+1].y<-lattice_constant/2)
		{
			coord[i][N+1].y=coord[i][N+1].y+lattice_constant;
		}		

		//cout << "i: " << coord[i][N+1].x << " , " << coord[i][N+1].y << endl;
		//cout << coord[i][N+1].x << " " << coord[i][N+1].y << endl;
		if (abs(coord[i][N+1].x)>lattice_constant/2 or abs(coord[i][N+1].y)>lattice_constant/2 or kin_energy[N]-kin_energy[N-1]>10 or acceleration[i][N].x-acceleration[i][N-1].x>1 or acceleration[i][N].y-acceleration[i][N-1].y>1 )
		{
			cout << "РАСХОДИМОСТЬ" << endl;
			stop = 1;
			cout << "velocity_x_" << i << "_" << N-1 << "= " <<  velocity[i][N-1].x << endl;			
			cout << "velocity_y_" << i << "_" << N-1 << "= " <<  velocity[i][N-1].y << endl;
			cout << "velocity_x_" << i << "_" << N << "= " <<  velocity[i][N].x << endl;			
			cout << "velocity_y_" << i << "_" << N << "= " <<  velocity[i][N].y << endl;		
			cout << "acceleration_x_" << i << "_" << N-1 <<  " = " << acceleration[i][N-1].x << endl;	
			cout << "acceleration_y_" << i << "_" << N-1 <<  " = " << acceleration[i][N-1].y << endl;
			cout << "acceleration_x_" << i << "_" << N <<  " = " << acceleration[i][N].x << endl;	
			cout << "acceleration_y_" << i << "_" << N <<  " = " << acceleration[i][N].y << endl;
			cout << "SUM VEL = " << k << endl;
			cout << "o = " << o << endl;
			break;
		}
			if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0)
			{
			write <<  coord[i][N].x << " " << coord[i][N].y << endl;
			if (i==14)
			{
				partwrite <<  coord[i][N].x << " " << coord[i][N].y << endl;
			}
		}
		//cout << coord[i][N].x << endl;
		//cout << coord[i][N].y << endl;
		//R+= pow((coord[i][N].x-coord[i][N-1].x),2)+pow((coord[i][N].y-coord[i][N-1].y),2);
		dx+=coord[i][N].x-coord[i][0].x;//velocity[i][N].x*dt+0.5*acceleration[i][N].x*pow(dt,2);
		dy+=coord[i][N].y-coord[i][0].y;//velocity[i][N].y*dt+0.5*acceleration[i][N].y*pow(dt,2);
		R+= pow(dx,2)+pow(dy,2);
		Z+=sqrt(pow(velocity[i][0].x+velocity[i][0].x,2))*sqrt(pow(velocity[i][N].y+velocity[i][N].y,2));
	}
		R=R/number_of_particles;
		Z=Z/number_of_particles;
		rwrite << N << " " << R << endl;
		zwrite << N << " " << Z << endl;
		//cout << "R MEOW =" << R << endl;
		R=0;
		//dx=0;
		//dy=0;
		//cout << "" << endl;
		//cout << "" << endl;	
		if (N % step == 0 or (N % step) % step ==0 or ((N % step) % step)%step == 0)
		{
		write << "" << endl;
		write << "" << endl;
		}


	//total energy = kin_energyх[0]+kin_energy[0]

	for (int i=0; i<number_of_particles; i++)
	{
		k+=velocity[i][N].x+velocity[i][N].y;
	}
	pulse << N << " " << k <<endl;
	k=0;
	cout << "" << endl;
	cout << "step " << N << ": kin energy = " << kin_energy[N] << endl;
	cout << "step " << N << ": pot energy = " << pot_energy[N] << endl;
	cout << "step " << N << ": total energy = " << kin_energy[N]+pot_energy[N] << endl;
	cout << "step " << N << ": temperature = " << temperature << endl;
	cout << "step " << N << ": pulse_sum = " << k << endl;
	//cout << "TOTAL energy " << N << " = " << pot_energy[N] + kin_energy[N] << endl;
	kinwrite << N << " " << kin_energy[N] << endl;
	potwrite <<  N << " " <<  pot_energy[N] << endl; 
	fullwrite << N  << " " << kin_energy[N]+pot_energy[N] << endl;  
	avtempwrite << N << " " << temperature_average << endl;
	tempwrite << N << " " << temperature << endl;
	//cout << "kin energy_" << 0 << " = " << kin_energy[0]  << endl;
	//cout << "pot energy_" << 0 << " = " << pot_energy[0] << endl;
	//cout << "TOTAL energy_" << 0 << " = " << pot_energy[0] + kin_energy[0]  << endl;
	//cout << "temperature_average" << " = " << temperature_average  << endl;

	
}
write.close();
kinwrite.close();
potwrite.close();
fullwrite.close();
tempwrite.close();
avtempwrite.close();
partwrite.close();
rwrite.close();
cout << "SUM VEL = " << k << endl;
cout << "o = " << o << endl;
	return 0;
}
