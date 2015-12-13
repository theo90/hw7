#include <cmath>
#include <fstream>

using namespace std;

void f ( double *ki, double x0,double y0,double x1,double y1, double t);
void RKstep(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,double *k7,  double* x, double t, double dt);
void RK4(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6, double *k7, double t,  double dt, double *x);
void RK5( double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,double *k7,double t,  double dt, double* y);
void max_bestimmen (double *x, double *y,  double &maximum);
int main()
{
	const int dim=4;
	double k1[dim], k2[dim], k3[dim], k4[dim], k5[dim], k6[dim], k7[dim];
	double dt=1e-3, tol=1e-5;
	double x[dim], y[dim];
	x[0]=0.994; x[1]=0.; x[2]=0. ;x[3]=-2.00158510637908;
	y[0]=0.994; y[1]=0.; y[2]=0. ;y[3]=-2.00158510637908;

	double T=17.065216560157, t=0.0, dt_new;
	double y_vor[dim];
	double maximum;
	ofstream out("RK_new.txt");
	while(t<T)
	{
		for(int i=0; i<dim; i++)
			y_vor[i]=y[i];
		RK4(k1, k2, k3, k4, k5, k6, k7, t, dt, y);
		RK5(k1, k2, k3, k4, k5, k6, k7, t, dt, x);
		max_bestimmen(x, y, maximum);
		for(int i=0; i<dim; i++)
			y[i]=y_vor[i];
		dt_new=dt*pow((tol/maximum), 0.2);
		dt=dt_new;
		RK4(k1, k2, k3, k4, k5, k6, k7, t, dt, y);
		for(int i=0; i<dim; i++) x[i]=y[i];
		t+=dt;
		out<<t<<"\t " <<y[0] <<"\t "<<y[1]<<"\t"<<dt <<endl;
	}
	out.close();
    
    return 0;
}
void f ( double *ki, double x0,double y0,double x1,double y1, double t)
{
	double const mu=0.012277471;
	double r=sqrt(pow((x0+mu),2)+pow(y0,2));
	double s=sqrt(pow((x0-1+mu),2)+pow(y0,2));
	ki[0]=x1;
	ki[1]=y1;
	ki[2]=x0+2*y1-((1-mu)*(x0+mu))/(pow(r,3))-mu*(x0-1+mu)/(pow(s,3));
	ki[3]=y0-2*x1-((1-mu)*y0/(pow(r,3)))-mu*y0/(pow(s,3));
}

void RKstep(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,double *k7,  double* x, double t, double dt)
{
	double x0, y0, x1, y1;
	x0=x[0]; y0=x[1]; x1=x[2]; y1=x[3];
	double const b11=3./40. , b12=9./40. , b10=3./10. ;
	double const b20=4./5. , b21=44./45. , b22=-56./15. , b23=32./9. ;
	double const b30=8./9. , b31=19372./6561. , b32=-25360./2187. , b33= 64448./6561. , b34=-212./729. ;
	double const b40=1. , b41=9017./3168. , b42=-355./33. , b43=46732./5247., b44=49./ 176. , b45=-5103./18656. ;
	double const b50=1. , b51=35./384. , b52=0, b53= 500. /1113. , b54=125./192. , b55=-2187./6784. , b56=11./84. ;
	
	
	f(k1, x0, y0, x1, y1, t);
	f(k2, x0+1./5.*dt*k1[0], y0+1./5.*dt*k1[1], x1+1./5.*dt*k1[2], y1+1./5.*dt*k1[3], t+dt*1./5.);
	f(k3, x0+dt*(b11*k1[0]+b12*k2[0]), y0+dt*(b11*k1[1]+b12*k2[1]), x1+dt*(b11*k1[2]+b12*k2[2]),
		y1+(b11*k1[3]+b12*k2[3])*dt, t+dt*b10);
	f(k4, x0+dt*(b21*k1[0]+b22*k2[0]+b23*k3[0]), y0+dt*(b21*k1[1]+b22*k2[1]+b23*k3[1]),
		x1+dt*(b21*k1[2]+b22*k2[2]+b23*k3[2]), y1+dt*(b21*k1[3]+b22*k2[3]+b23*k3[3]), t+dt*b20);

	f(k5, x0+dt*(b31*k1[0]+b32*k2[0]+b33*k3[0]+b34*k4[0]), y0+dt*(b31*k1[1]+b32*k2[1]+b33*k3[1]+b34*k4[1]),
		x1+dt*(b31*k1[2]+b32*k2[2]+b33*k3[2]+b34*k4[2]), y1+dt*(b31*k1[3]+b32*k2[3]+b33*k3[3]+b34*k4[3]),
		t+dt*b30);
	f(k6, x0+dt*(b41*k1[0]+b42*k2[0]+b43*k3[0]+b44*k4[0]+b45*k5[0]),
		y0+dt*(b41*k1[1]+b42*k2[1]+b43*k3[1]+b44*k4[1]+b45*k5[1]),
		x1+dt*(b41*k1[2]+b42*k2[2]+b43*k3[2]+b44*k4[2]+b45*k5[2]),
		y1+dt*(b41*k1[3]+b42*k2[3]+b43*k3[3]+b44*k4[3]+b45*k5[3]), t+dt*b40);
	f(k7, x0+dt*(b51*k1[0]+b52*k2[0]+b53*k3[0]+b54*k4[0]+b55*k5[0]+b56*k6[0]),
		y0+dt*(b51*k1[1]+b52*k2[1]+b53*k3[1]+b54*k4[1]+b55*k5[1]+b56*k6[1]),
		x1+dt*(b51*k1[2]+b52*k2[2]+b53*k3[2]+b54*k4[2]+b55*k5[2]+b56*k6[2]),
		y1+dt*(b51*k1[3]+b52*k2[3]+b53*k3[3]+b54*k4[3]+b55*k5[3]+b56*k6[3]),t+dt*b50);
}
void RK4(double *k1, double *k2, double *k3, double *k4, double *k5, double *k6, double *k7, double t,  double dt, double *x)
{
	const int dim=4;
	RKstep(k1, k2, k3,k4, k5, k6, k7, x,t,dt);
	for(int i=0; i<dim; i++)
		x[i]+=dt*(35./384.*k1[i]+500./1113*k3[i]+125./192*k4[i]-2187./6784.*k5[i]+11./84.*k6[i]);
}

void RK5( double *k1, double *k2, double *k3, double *k4, double *k5, double *k6,double *k7,double t,  double dt, double* y)
{
	const int dim=4;
	RKstep(k1, k2, k3,k4, k5, k6, k7, y, t, dt);
	for(int i=0; i<dim; i++)
		y[i]+=dt*(5179./57600*k1[i]+7571./16695.*k3[i]+393./640.*k4[i]-92097./339200.*k5[i]+187./2100.*k6[i]+1./40.*k7[i] );

}

void max_bestimmen (double *x, double *y,  double &maximum)
{
	
	const int dim=4;
	double norm[dim];
	for(int i=0; i<dim; i++)
		norm[i]=abs(y[i]-x[i]);

    maximum=norm[0];
	for(int i=0; i<dim; i++)
	{
		if(norm[i]>maximum)
			maximum=norm[i];
	}
}
