%% 1a
a = (1-1) + 1e-16
b = 1 - (1 + 1e-16)

%% 1b
H_5 = hilb(5);
InvH_5 = inv(H_5);
e_5 = ones(5,1);
H = norm(H_5*(InvH_5*e_5) - e_5)

H_10 = hilb(10);
InvH_10 = inv(H_10);
e_10 = ones(10,1);
H = norm(H_10*(InvH_10*e_10) - e_10)

H_20 = hilb(20);
InvH_20 = inv(H_20);
e_20 = ones(20,1);
H = norm(H_20*(InvH_20*e_20) - e_20)

%% 1c
f = @(x)(exp(x)/(power(cos(x), 3) + power(sin(x), 3)));

x = linspace(0,100);
y = zeros(length(x));
for i = 1:100
    y(i) = f(x(i));
end

x0 = pi/4;
f_prime = 3.101766393836051;
xy = zeros(16, 2);

for k = 1:16
    h = 10^-k;
    f_changed = f(x0 + 10^-k);
    f_(k) = (f_changed - f(x0))/(10^-k);
    diff = abs(f_prime - f_(k));
    fprintf('diff%2d = %5f\n', k, diff) 
    xy(k,2) = diff;
    xy(k,1) = h;
end

loglog(xy(:,1), xy(:,2))
grid on
%% 2a
A = [1 -2;3 2];
Inverse_A = inv(A);

norm(A,1);
norm(A,2);
norm(A,inf);
norm(Inverse_A,1);
norm(Inverse_A,2);
norm(Inverse_A,inf);

K1 = cond(A,1)
K2 = cond(A,2)
K_infinity = cond(A, inf)

%% 2c
n1 = 100;
A = rand(n1);
fprintf('Time for my functions for one and infinity norm for n = 100\n')
tic
    one_norm(A);
toc
tic
    infinity_norm(A);
toc
fprintf('Time for Matlab functions for one and infinity norms for n = 100\n')
tic
    norm(A,1);
toc
tic
    norm(A,inf);
toc

for k = 1:6
    n1 = 2*n1;
    A = rand(n1);
    fprintf('Time for my functions for one and infinity norm for n = %6d\n', n1)
    tic
    one_norm(A);
    toc
    tic
    infinity_norm(A);
    toc
    fprintf('Time for Matlab functions for one and infinity norms for n = %6d\n', n1)
    tic
    norm(A,1);
    toc
    tic
    norm(A,inf);
    toc
end
    
%calculate norm(A,one)
function y = one_norm(X)
    y = max(sum(abs(X)));
end

%calculate norm(A,infinity)
function z = infinity_norm(X)
    z = max(sum(abs(X')));
end

%% 2d

%% 7
A = [9,-6;12,-8;0,20]
[Q R] = qr(A,0)
b = [300, 600, 900];
b = b'
x = inv(R)*Q'*b

f = @(x)((-300 + 9*x)/6);
g = @(x)((-600 + 12*x)/8);
h = @(x)(900/20);

xx = linspace (0,100);
N = length(xx);
yy = zeros(N);
zz = zeros(N);
ww = zeros(N);
for i =1:N
    yy(i) = f(xx(i));
    zz(i) = g(xx(i));
    ww(i) = h(xx(i));
end

figure(1)
hold on
plot(xx, yy,'b', xx, zz, 'k', xx, ww, 'r')
plot(74, 45, '*')

%% 8
X = [0 .5 1 1.5 2 2.5];
X3 = X.^3;
X2 = X.^2;
X1 = X;
X0 = ones(1,6);
X3 = X3';
X2 = X2';
X1 = X1';
X0 = X0';

X = [X3 X2 X1 X0]
Y = [0 .2 .27 .3 .32 .33]'

[Q R] = qr(X,0)
x = inv(R)*Q'*Y

X'
A = X'*X
b = X'*Y
x = A\b
x = inv(A)*b