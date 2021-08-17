
1. The function to evaluate the integrand will receieve a parameter of type double and will return a double as well.

The best way to do this assignment is for one person to share their screen and write the code while discussing the questions; the rest of the group should be copying the final code into their own files, and perhaps helping debug the code on their machine if issues arise.

2. The function to calculate the integral will include parameter list with bounds of integral a,b as double, als as part of the parameter list we will include npoints type int and for return type will be double after.

Each person should submit a pull request with their copy at the end of the day.

3. The width of the rectangles will evaluate to (b-a) / nPoints.

4. For every increment of dx, we need to find midpoint between the current and next dx value. Where dx is defined as (b-a) / npoints, for all points i (nPoints). From a to b, every idx + dx/2 . The other way would be as (idx +(1+i)dx)/2 

5. Write function integration function and to calculate total area, we summed all rectangle areas using formula the width we calculates in step 3 as dx = (b-a) / npoints

The integraiton function calculates total area, by summing all retangle areaas. For width we are using dx = (b-a) / npoints
for height fx = 1 /(1+x^2)

rectangle area = fx * dx

sum of all rectangleAreas

6. The main function handles the call to integration function and printing result. Between 188 and 189 for 5 dec places.



This is a group project, to be written in C++. The group should discuss the answers to the questions and write the README for the project repository.

1. The function to evaluate the integrand will receieve a parameter of type double and will return a double as well.

The best way to do this assignment is for one person to share their screen and write the code while discussing the questions; the rest of the group should be copying the final code into their own files, and perhaps helping debug the code on their machine if issues arise.

2. The function to calculate the integral will include parameter list with bounds of integral a,b as double, als as part of the parameter list we will include npoints type int and for return type will be double after.

Each person should submit a pull request with their copy at the end of the day.

3. The width of the rectangles will evaluate to (b-a) / nPoints.

4. Fore


Write a C++ function to evaluate the integrand 11+x2

. This function should take x as the only argument and return the value. Note C++ does not have the ** operator for exponentiation. Just multiply x by itself. Q: What type should be used for x and the return type?

You will be writing a function that takes the bounds of the integral a,b and the desired number of integration points, and returns the approximated value of the integral. Q: What should the types of a,b and npoints be? What type should the function return?

Given the bounds of integration a,b and the number of points, what is the width of the rectangles you will be using? Write this in terms of a,b, and npoints. We will call this dx

At what points will you need to evaluate the integrand? This should be written in terms of a,b, and dx.

Now you can write the integration function. You will not need to precompute the integration points - you can determine them on the fly. Write a for loop over the number of points, and evaluate the integrand using the function you wrote in Step 1 at the points determined in the previous question. For each point, calculate the area of a rectangle centered on that point, and sum these areas all the areas into a variable. This variable represents the approximated value of the integral.

Call this integration function from main and print the result. Compare, for example, with your first homework by calling your function with a=0, b=1, and some number of points (50 or so). The result should approach π4
as the number of points is increased. About how many points does it take to start converging to the proper result (4 decimal places)?


We have now seen the Monte-Carlo method and the Riemann Sum method for approximating integrals. You will notice that the number of integration points for the Riemann Sum method is much less than the Monte-Carlo method. Why would you ever use the Monte-Carlo method?