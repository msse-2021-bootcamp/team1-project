#include <iostream>
using namespace std;
double calcint(double x)
{

    return 1/(1 + x*x);

}



int main() 
{
    
    int a = 0;
    int b = 1;
    int n = 50;
    


    double approx = calcint(a,b,n);
    cout << "approx: " << approx;
    cout << "pi: " << (approx/4);

}


double calcappint(int a, int b, int n); 
{
    double dx = (b-a) / n;
    
    double totarea = 0;

    for(int i = 0; i <= n; i++); {
        double midpt = idx + dx/2;
        double fx = calculate(midpt);
        totarea += dx * fx
    }
    return totarea
}
