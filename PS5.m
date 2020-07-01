%% 2_c
for n = [5, 10, 20, 30]
    xvec = zeros(n,1);
    v = zeros(n+1);
    for i = 1:n+1
        x = -1 + (2*(i-1))/n;
        xvec(i) = x;
        for j = 1:n+1
            v(i,j) = x^(j-1);
        end
    end
    v;
    xvec;
    cond(v)
end

%% 3_a
X = [0 .5 1 1.5 2 2.5];
Y = [0 .2 .27 .3 .32 .33];
x = linspace(-1,4);
y = zeros(length(x));
X3 = X.^3;
X2 = X.^2;
X1 = X;
X0 = ones(1,6);
X3 = X3';
X2 = X2';
X1 = X1';
X0 = X0';

A1 = [X3 X2 X1 X0]
b = [0 .2 .27 .3 .32 .33]'
sol1 = A1\b
f = @(x) sol1(4) + sol1(3)*x + sol1(2)*x.^2 + sol1(1)*x.^3;

A1(6, :) = [];
A2 = A1
b = [0 .2 .27 .3 .32]'
sol2 = A2\b
g = @(x) sol2(4) + sol2(3)*x + sol2(2)*x.^2 + sol2(1)*x.^3;

A2(5, :) = [];
A3 = A2
b = [0 .2 .27 .3]'
sol3 = A3\b
h = @(x) sol3(4) + sol3(3)*x + sol3(2)*x.^2 + sol3(1)*x.^3;

figure(1)
hold on
plot(X, Y, '*', x, f(x), 'b', x, g(x), 'g', x, h(x), 'r')
legend('data points', 'y = .0033x^3 + .4970x^2 - .2738x +.0511', 'y = .0009x^3 + .5486x^2 - .3543x +.0800', 'y = .5900x^2 - .4400x +.1200')

%% 4_a
x = linspace(0,1);
f  =@(x)exp(3*x);
nodes = [0,.5,1];
c = polyfit(nodes, f(nodes), 2);
p2 = @(x)c(1)*x.^2 + c(2)*x + c(3);
plot(x, f(x), 'b', x, p2(x), 'g', nodes, f(nodes), 'ro');
legend('f(x)', 'p2(x)', 'Interpolation Points');

%% 4_b
x = linspace(0,1);
f  =@(x)exp(3*x);

nodes = [0,.5,1];
x0 = nodes(1);
y0 = f(x0);
x1 = nodes(2);
y1 = f(x1);
x2 = nodes(3);
y2 = f(x2);

L0 = @(x)((x-x1).*(x-x2)) / ((x0 - x1).*(x0 - x2));
L1 = @(x)((x-x0).*(x-x2)) / ((x1 - x0).*(x1 - x2));
L2 = @(x)((x-x0).*(x-x1)) / ((x2 - x0).*(x2 - x1));
p2 = @(x) y0*L0(x) + y1*L1(x) + y2*L2(x);

plot(x, f(x), 'k', x, L0(x), 'r', x, L1(x), 'g', x, L2(x), 'b', x, p2(x), 'm')
legend('f(x)', 'L0(x)', 'L1(x)', 'L2(x)', 'p2(x)');

%% 4_c
%error bound = sup [f^(3)(x)/3!]*(x-x0)*(x-x1)*(x-x2)
%since f(x) = exp(3x), f(3)(x) = 27exp(3x) which is maximized at x=1
%thus error bound = [(27e^3)/6]*(x-x0)*(x-x1)*(x-x2)
est_error = @(x)((27*exp(3))/6).*(x-x0).*(x-x1).*(x-x2);
error = @(x)f(x) - p2(x);
err = error(.75)
est_err = est_error(.75)
comparison = err - est_err

%% 4_d
x = linspace(0,1);
f  = @(x)exp(3*x);
g = @(x)3*exp(3*x);

nodes = [0,.5,1];
x0 = nodes(1);
y0 = f(x0);
z0 = g(x0);
x1 = nodes(2);
y1 = f(x1);
z1 = g(x1);

L0sq = @(x)((x-x1) / (x0-x1)).^2;
L1sq = @(x)((x-x0) / (x1-x0)).^2;
L0prime = 1/(x0-x1);
L1prime = 1/(x1-x0);

H0 = @(x)L0sq(x).*(1 - (2.*L0prime.*(x-x0)));
H1 = @(x)L1sq(x).*(1 - (2.*L1prime.*(x-x1)));
K0 = @(x)L0sq(x).*(x-x0);
K1 = @(x)L1sq(x).*(x-x1);

p3 = @(x) y0.*H0(x) + z0.*K0(x) + y1.*H1(x) + z1.*K1(x);

figure(1)
plot(x, f(x), 'k', x, p3(x), 'm')
legend('f(x)', 'p3(x)');
figure(2)
plot( x, H0(x), 'r', x, K0(x), 'g', x, H1(x), 'b', x, K1(x), 'y')
legend('H0(x)', 'K0(x)', 'H1(x)', 'K1(x)')

%% 5
figure(1)
xlabel('Number of Iterations')
ylabel('r')
hold on
for r = 0:.2:1
    A = [1 r;r 1]
    [U, cnt] = qr_algo(A)
    plot(cnt, r, '*')
end
   
function [A1, cnt] = qr_algo(A)
    if A == A'
        cnt = 0;
        tau = 1e-10;
        sortEig = sort(eig(A));
        [Q R] = qr(A, 0);
        A1 = R*Q;
        sortDiag = sort(diag(A1));
        while abs(sortEig - sortDiag) > tau
            cnt = cnt + 1;
            [Q R] = qr(A, 0);
            A1 = R*Q;
            sortDiag = sort(diag(A1));
            A = A1;
        end
    else
        A1 = 0;
        cnt = 0;
        fprintf('Not a symmetric matrix')
        return
    end
end