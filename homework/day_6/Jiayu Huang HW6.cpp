#include<iostream>

double F(double x)
{
    double a;
    a = 1 / (1 + x * x);
    return a;
}

void defIntegral(double a, double b, int n)
{
    double x = a;
    double dx = (b - a) / n;
    double s, S = 0;

    for (int i = 1; i <= n; i++)
    {
        s = F((x + x + dx) / 2) * dx;
        S += s;
        x += dx;
    }
    std::cout << s << std::endl;
}

int main(void)
{
    std::cout << "input between [0,1]" << std::endl;
    defIntegral(0, 1, 50);
    defIntegral(0, 1, 100);
    defIntegral(0, 1, 500);
    defIntegral(0, 1, 1000);
    defIntegral(0, 1, 10000);
    defIntegral(0, 1, 100000);
    return 0;
}
