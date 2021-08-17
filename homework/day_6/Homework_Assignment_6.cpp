#include <iostream>
using namespace std;






double calcint(double x)
{

    return 1/(1 + x*x);

}

double calcappint(double a, double b, int npt) 
{
    double dx = (b-a) / npt;
    
    double totarea = 0;

    for(int i = 0; i < npt; i++) {
        double midpt = i*dx + dx/2;
        double fx = calcint(midpt);
        totarea += dx * fx;
    }
    return totarea;
}



int main() 
{
    
    int a = 0;
    int b = 1;
    int n = 900;

 
    


    double approx = calcappint(a,b,n);
    std::cout << "approx: " << approx << std::endl;
    std::cout << "pi: " << (approx*4) << std::endl;


}