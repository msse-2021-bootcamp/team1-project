#include <iostream>

double eval_integrand(double x)
{
	return 1 / (1 + x * x);
}

double integration_function(int a, int b, int points){
	
	double integeral_sum = 0;
	double rect_size = ((double)b-a)/points;

	for(int i=0; i < points; i++)
	{
		double step = eval_integrand((double)a+(rect_size/2)+(double)i*rect_size);
		double midpoint = (double)a+(rect_size/2)+(double)i*rect_size;
		std::cout << midpoint  << std::endl;
		integeral_sum += step;
	}

	return ((double)b-a)/points * integeral_sum;

}

int main(void)
{
	//std::cout << eval_integrand(1) << std::endl;
	std::cout << "(n=10) " << integration_function(0,1,10) << std::endl;
	//std::cout << "(n=50) " << integration_function(0,1,50) << std::endl;
	//std::cout << "(n=100) " << integration_function(0,1,100) << std::endl;
	//std::cout << "(n=500) " << integration_function(0,1,500) << std::endl;
	//std::cout << "(n=1000) " << integration_function(0,1,1000) << std::endl;
	//std::cout << "(n=5000) " << integration_function(0,1,5000) << std::endl;
	//std::cout << "(n=10000) " << integration_function(0,1,10000) << std::endl;
	//std::cout << "(n=50000) " << integration_function(0,1,50000) << std::endl;

}
