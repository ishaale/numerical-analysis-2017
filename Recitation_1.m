% Feb 3, 2017
% Recitation 1
% Ana Perez-Gea and Jiajun Tong

% We study the equation f(x) = exp(x)-2*x-1 = 0.
% We first check that there exists a solution on [1,2].
% Then we apply fixed point iteration and the Newton's
% method to find out the solution.

clc
clear
close all

%% Existence of root
% define function
f = @(x)(exp(x)-2*x-1);
    % This syntax is good for defining simple functions within an m-file.
    % The following syntax is OK for Matlab R2016b or later versions
            % function y = f(x)
            % 	y = exp(x)-2*x-1;
            % end
    % For more complicated function, you can define it in a separate
    % m-file with the function name = the file name, and store it in the
    % same path as the current m-file.
    
% there exists a root in [1,2]
f(1) % f(1) < 0
f(2) % f(2) > 0

%% Plot the function f on [1,2]
% generate equally spaced points from 1 to 2
help linspace
xx = linspace(1,2); % default of N=100 points 
    % linspace generates row vectors
    % It is equivalent to use
        % xx = linspace(1,2,100);
        
N = length(xx); % length returns the length of the largest array dimension
yy = zeros(N);
% Matlab runs faster if we pre-define a matrix/vector
% with known number of entries instead of 
% repeatedly changing its size in the loop

for i = 1:N
	yy(i) = f(xx(i));
end
    % A tricky point: in fact you can write yy = f(xx) in the place of the
    % whole for-loop, since the definition of f above is also well-defined
    % as a function of vectors. But one needs to be very cautious. 
    % It is NOT recommended for beginners.

% we plot
close all
figure(1)
plot(xx,yy,'b')
hold on 
% hold on means we do not erase what has been plotted when new plots are
% made. It can be used to a figure before we make any plot.
plot(xx,zeros(N),'k')

%% Fixed point method
g = @(x)(log(2*x+1));
% We know that solving f(x) = 0 is equivalent to solving the fixed point
% problem x = g(x).
    % Again, it is good to write the following in Matlab R2016b or later
    % versions.
        % function y = g(x)
        % 	y = log(2*x+1);
        % end

x0 = 1; % initial value
plot(x0,f(x0),'o')
iter = 10; % number of iterations to do

% % figure(1)
% % hold on

for i = 1:iter % a for-loop!
	x0 = g(x0);
	plot(x0,f(x0),'o')
end
plot(x0, f(x0), 'r*')

%     % If we would like to keep the history of the iteration, do the following.
%     iter_history = zeros(iter+1,1);
%     iter_history(1) = x0;
%     plot(iter_history(1),f(iter_history(1)),'o')
%     for i = 1:iter % a for-loop!
%     	iter_history(i+1)= g(iter_history(i));
%     	plot(iter_history(i+1),f(iter_history(i+1)),'o')
%     end
%     plot(iter_history(iter+1), f(iter_history(iter+1)), 'r*')

%% Newton's Method
% Let us start a new figure
figure(2)
hold on
plot(xx,yy,'b')
plot(xx,zeros(N),'k')

% we need the derivative h(x)=f'(x)
h = @(x)(exp(x)-2);
    % Again, it is good to write the following in Matlab R2016b or later
    % versions.
        % function y = h(x)
        % 	y = exp(x)-2;
        % end
        
% we take initial value random in [1,2]
help rand
x0 = 1+rand(1); 
plot(x0,f(x0),'k.', 'markersize', 15)
% if the solution changes less than tolerance, we stop
tol = 1e-6
% keep track of iterations
i = 0
x1 = x0-f(x0)/h(x0);
plot(x1,f(x1),'r.', 'markersize', 15)
i = i+1
while abs(x1-x0) > tol % another way of realizing loops
	x0 = x1;
	x1 = x0-f(x0)/h(x0);
	plot(x1,f(x1),'r.', 'markersize', 15)
	i = i+1
end
% Think about how to adapt the loop here to keep the history.