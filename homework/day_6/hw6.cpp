#include <iostream>
using namespace std;

double f(double x);
int main()
{
	
    double a;
    double b;
    int n;
    double h;
    double area;
    double areab;
 
    cout << "What is the left x value? ";
    cin >> a;
    cout << "What is the right x value? ";
    cin >> b;
    cout << "How many pieces do you want to use? ";
    cin >> n;
 
    h = (b-a)/n;
    area = (f(a)+f(b))/2;
     
 
    for(int i=1;i<n;i++)
		{
		area+= f(a+i*h);
		}
    areab = area*h;
    cout << "The area under the curve of f(x) = 1/(1 + x*x) is ";
    cout.precision(4);
	cout << areab << endl;
	cout << a << "<= x <=" << b << endl; 

}



double f(double x) 
	   {
		   return 1/(1+ x*x);
	   }