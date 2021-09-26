%% ODE Lab: Creating your own ODE solver in MATLAB
%
% In this lab, you will write your own ODE solver for the Improved Euler 
% method (also known as the Heun method), and compare its results to those 
% of |ode45|.
%
% You will also learn how to write a function in a separate m-file and 
% execute it.
% 
% Opening the m-file lab3.m in the MATLAB editor, step through each
% part using cell mode to see the results.  Compare the output with the
% PDF, which was generated from this m-file.
%
% There are six (6) exercises in this lab that are to be handed in on the
% due date. Write your solutions in the template, including
% appropriate descriptions in each step. Save the .m files and submit them 
% online on Quercus.
%
%% Student Information
%
% Student Name: Emily Traynor
%
% Student Number: 1005763372
%

%% Creating new functions using m-files.
%  
% Create a new function in a separate m-file:
%
% Specifics:  Create a text file with the file name f.m
% with the following lines of code (text):
%
%  function y = f(a,b,c) 
%  y = a+b+c;
%
% Now MATLAB can call the new function f (which simply accepts 3 numbers
% and adds them together).  
% To see how this works, type the following in the matlab command window:
% sum = f(1,2,3)

%% Exercise 1
%
% Objective: Write your own ODE solver (using the Heun/Improved Euler
% Method).
%
% Details: This m-file should be a function which accepts as variables 
% (t0,tN,y0,h), where t0 and tN are the start and end points of the 
% interval on which to solve the ODE, y0 is the initial condition of the
% ODE, and h is the stepsize.  You may also want to pass the function into
% the ODE the way |ode45| does (check lab 2).
%
% Note: you will need to use a loop to do this exercise.  
% You will also need to recall the Heun/Improved Euler algorithm learned in lectures.  
%


%NOTE: FUNCTIONS are at the bottom

%answerss = ImprovedE(0,0.5,1,0.01);

m = @(t,y) y^3 - t^2;
soln = ode45(m,[0,0.5],1);
disp(answerss);
disp(soln);
%plot(soln.x, soln.y, '--');





%% Exercise 2
%
% Objective: Compare Heun with |ode45|.
%
% Specifics:  For the following initial-value problems (from lab 2, 
% exercises 1, 4-6), approximate the solutions with your function from
% exercise 1 (Improved Euler Method).
% Plot the graphs of your Improved Euler Approximation with the |ode45| 
% approximation.
%
% (a) |y' = y tan t + sin t, y(0) = -1/2| from |t = 0| to |t = pi|
%
% (b) |y' = 1 / y^2 , y(1) = 1| from |t=1| to |t=10|
%
% (c) |y' =  1 - t y / 2, y(0) = -1| from |t=0| to |t=10|
%
% (d) |y' = y^3 - t^2, y(0) = 1| from |t=0| to |t=1|
%stops around 0.5
% Comment on any major differences, or the lack thereof. You do not need
% to reproduce all the code here. Simply make note of any differences for
% each of the four IVPs.

%for function a, the line follows closely except around t = 1.5 (because
% tan(pi/2) = INF), it causes the improved eulers method to spike. The
% larger the interval, the less likely a very steep gradient is taken but
% it makes the deviation from the actual value more obvious. The smaller
% the interval the larger the value aorund pi/2 becomes because each step
% generates a more and more negative number.

%for funciton b, both functions follow very closely, but it appears that
%the gap between them increases, suggesting that Euler's method accumulates
%error

%for function c, the improved eulers method produces a smoother line,
%however, both solutions follow very closely

%for funciton d, the Eulers method follows closely until t = 0.5, when
%ode45 reaches a vertical asymptote, but improved eulers method continues
%(in a kind of unpridictable manner)

%% Exercise 3
%
% Objective: Use Euler's method and verify an estimate for the global error.
%
% Details: 
%
% (a) Use Euler's method (you can use
% euler.m from iode) to solve the IVP
%
% |y' = 2 t sqrt( 1 - y^2 )  ,  y(0) = 0|
%
% from |t=0| to |t=0.5|.
%

E = @(t,y) 2*t*((1-y^2)^0.5);
t_at_tN = euler(E,0,0.5,0,25);
disp(t_at_tN)
%thus at t = 0.5, the euler value of y = 0.2427

%answer to initial value problem y = sin(x^2 + C1)
%input init values y = sin(t^2 + 0)   (could have different C1s)


% (b) Calculate the solution of the IVP and evaluate it at |t=0.5|.

%at t = 0.5,  y = 0.2499 degree

% (c) Read the attached derivation of an estimate of the global error for 
%     Euler's method. Type out the resulting bound for En here in
%     a comment. Define each variable.
%
%En <= (1+M)*(dt/2)*(e^(M*dt*n) - 1)
%M is the maximum value in the interval
%n is the step number
%dt is the step size

%first lets find M
%so basically M is bigger than f, the partial derivative of f with respect
%to time and y and M is bigger than 0
%f = 2t(1-y^2)^0.5
%df/dt = 2(11-y^2)^0.5
%df/dy = -2ty(1-y^2)^-0.5
%looking at f, it depends on t, lets see what value of t within [0,0.5]
%gives the largest M
%when t = 0.5, f is the largest so f = (1-y^2)^0.5
%but this is smaller than df/dt
%knowing the actual solution (y = sin(t^2)) we know y has the range [0,1]
%thus the max value of dy/dt = 2
%looking at df/dy, by subbing y = sin(t^2)and making it an absolute value
%the following eq is obtained:
%|df/dy| = -2tsin(t^2)/(1-(sin(t^2)^2))^0.5
%using desmos I find the largest value |df/dy| = 0.26
%df/dt produces a larger M
%thus M = 2

%thus En <= (3)*(dt/2)*(e^(2*dt*n) - 1)

% (d) Compute the error estimate for |t=0.5| and compare with the actual
% error.
%

%let the time step be dt = 0.01
%thus the step number is n = 50
%En <= (1 + 2)*(0.01/2)*(e^(2*0.01*50) - 1)
%En <= 0.0258

%looking at the actual value (0.2499) and the euler value (0.2427)
%the difference between them is 0.0072
%0.0072<0.0258
%thus the error estimate of 0.0258 holds true

% (e) Change the time step and compare the new error estimate with the
% actual error. Comment on how it confirms the order of Euler's method.

%im going to increase the step by 2 (dt=0.02)
%Eulers value = 0.2379
%actual value = 0.2499
%new En = (1 + 2)*(0.02/2)*(e^(2*0.02*25) - 1) = 0.0515
%difference in values = 0.012
%0.012<0.0515 thus En holds true
%and also 0.0515 is roughly double 0.0258 (the previous En)
%this confirms the order of Eulers method because the time is 
%proportional to the error bound


%% Adaptive Step Size
%
% As mentioned in lab 2, the step size in |ode45| is adapted to a
% specific error tolerance.
%
% The idea of adaptive step size is to change the step size |h| to a
% smaller number whenever the derivative of the solution changes quickly.
% This is done by evaluating f(t,y) and checking how it changes from one
% iteration to the next.

%% Exercise 4
%
% Objective: Create an Adaptive Euler method, with an adaptive step size |h|.
%
% Details: Create an m-file which accepts the variables |(t0,tN,y0,h)|, as 
% in exercise 1, where |h| is an initial step size. You may also want to 
% pass the function into the ODE the way |ode45| does.
%
% Create an implementation of Euler's method by modifying your solution to 
% exercise 1. Change it to include the following:
%
% (a) On each timestep, make two estimates of the value of the solution at
% the end of the timestep: |Y| from one Euler step of size |h| and |Z| 
% from two successive Euler steps of size |h/2|. The difference in these
% two values is an estimate for the error.
%
% (b) Let |tol=1e-8| and |D=Z-Y|. If |abs(D)<tol|, declare the step to be
% successful and set the new solution value to be |Z+D|. This value has
% local error |O(h^3)|. If |abs(D)>=tol|, reject this step and repeat it 
% with a new step size, from (c).
%
% (c) Update the step size as |h = 0.9*h*min(max(tol/abs(D),0.3),2)|.
%
% Comment on what the formula for updating the step size is attempting to
% achieve.

%it looks at how close the tolerance and difference are
%if they are similar, tol/dif creates a number larger than 0.3
%the max is picked and that number is compared to 2
%since D is larger than tol, tol/abs will never be more than 1
%it is then multiplied by 0.9 and h to get a slightly smaller h
%a reduced step size should reduce the error



%% Exercise 5
%
% Objective: Compare Euler to your Adaptive Euler method.
%
% Details: Consider the IVP from exercise 3.
%
% (a) Use Euler method to approximate the solution from |t=0| to |t=0.75|
% with |h=0.025|.
myeuler(0,0.75,0,0.025)
%
% (b) Use your Adaptive Euler method to approximate the solution from |t=0| 
% to |t=0.75| with initial |h=0.025|.
adaptive(0,0.75,0,0.025)
% (c) Plot both approximations together with the exact solution.

tt = linspace(0,0.75,30);
yy = sin(tt.^2);

plot(tt, yy, 'color', 'g');

%% Exercise 6
%
% Objective: Problems with Numerical Methods.
%
% Details: Consider the IVP from exercise 3 (and 5).
% 
% (a) From the two approximations calculated in exercise 5, which one is
% closer to the actual solution (done in 3.b)? Explain why.
% %the adaptive euler is much closer. this is because it can pick a step
% size most suitable based on an error detected
%as it reduces step size the error also reduces proporitonally
%this also lowers the accumulated error

% (b) Plot the exact solution (from exercise 3.b), the Euler's 
% approximation (from exercise 3.a) and the adaptive Euler's approximation 
% (from exercise 5) from |t=0| to |t=1.5|.

myeuler(0,3,0,0.025)
adaptive(0,3,0,0.025)
tt = linspace(0,3,120);
yy = sin(tt.^2);
plot(tt, yy, 'color', 'g');

% (c) Notice how the exact solution and the approximations become very
% different. Why is that? Write your answer as a comment.

%the euler and adaptive methods become unreliable near 1 and 1.2,
%respectively
%the adaptive method becomes unreliable sooner. 
%after looking at different h values for the adaptive method the faster the
%gradient of the exact solution changes, the lower the h has to be for the
%adaptive method to remain similar
%however both methods do not follow the exact after x=sqrt(pi/2)
%the gradient goes from positive to negative
%I also looked at the tolerance
%the reason adaptive becomes unreliable before eulers is since the
%tolerence is too restrictive

%my best guest for y the approximations are different is that the
%approximations are based on the actual first order derivative and they
%move by a finite dt, where as the exact solution is derived from a dt that
%is infinitesimally small dt
%dy =  2*t*(1-(y^2))^0.5;
%the approximates should start calculating complex numbers but matlab
%ignores the complex part
%this is what causes the deviaiton from the exact
%.....I tried......







%///////////////////////////////////////////////////////////////////////
%///////////////////////////////Functions///////////////////////////////
%///////////////////////////////////////////////////////////////////////

function estimate = ImprovedE(t0,tN,y0,h)
N = (tN-t0)/h;
N = round(N);
t = t0;
y = y0;
ts = zeros(N,1);
ys = zeros(N,1);
n=1;
ts(n) = t;
ys(n) = y;
while t <= tN
    ey = y + (h*g(t,y)); %euler's y
    a = h*ey; %area under graph (horizontal line/rectangle)
    ny = y+a; %new y
    Trap = (h/2)*(g(t,y) + g((t+h),ny));
    y = y + Trap;
    t = t + h;
    %estimate = t,y;
    n = n+1;
    ts(n) = t;
    ys(n) = y;
    
end
estimate = [ts,ys];
plot(ts,ys)
hold on
end

%///////////////////////////////////////////////////////////////////////

function [y t] = euler(f,a,b,ya,n)
h = (b - a) / n;
y(1,:) = ya;
t(1) = a;
for i = 1 : n
    y(i+1,:) = y(i,:) + h * f(t(i),y(i,:));
    t(i+1) = t(i) + h;
end
%plot(t,y)
end

function outputs = myeuler(t0,tN,y0,h)
N = (tN-t0)/h;
N = round(N);
ts = zeros(N:1);
ys = zeros(N:1);
n=1;
ts(n) = t0;
ys(n) = y0;
y = y0;
t = t0;
while t <= tN
    y = y + h*g(t,y);
    t = t+h;
    n = n+1;
    ts(n) = t;
    ys(n) = y;
end

plot(ts,ys, '--')
hold on
outputs = [ts,ys]
end




%////////////////////////////////////////////////////////////////////

function output = adaptive(t0,tN,y0,h)
N = (tN-t0)/h;
N = round(N);
ts = zeros(N,1);
ys = zeros(N,1);

t = t0;
y = y0;
dt = h;

n = 1;
ts(n) = t;
ys(n) = y;

tol = exp(-8)
while t <= tN
    fullstep = y + dt*g(t,y); % Y
    halfstep = y + (dt/2)*g(t,y);
    nexthalfstep = halfstep + (dt/2)*g((t+(dt/2)),halfstep); %Z
    difference = nexthalfstep - fullstep; %D
    if abs(difference) < tol
        y = nexthalfstep + difference;
        t = t + h;
        n = n+1;
        ts(n) = t;
        ys(n) = y;
        dt = h;
    else
        dt = 0.9*dt*min(max(tol/abs(difference),0.3),2);
    end
end

output = [ts,ys];
plot(ts,ys)
hold on
    
end


function dy = g(t,y)
dy =  2*t*(1-(y^2))^0.5;
end
