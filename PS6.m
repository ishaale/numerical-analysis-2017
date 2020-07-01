%% 1a
f = @(x)cos(3*pi*x).^4;
int_f = 3/8;
g = @(x)sqrt(x);
int_g = 2/3;

n = 100;
error_f = zeros(n,2);
error_g = zeros(n,2);
m = zeros(n,1);

for i = 1:n
    m(i,1) = 10*i;
    
    trap_f = trapez(f, 0, 1, 10*i);
    simp_f = simpson(f, 0, 1, 10*i);
    error_f(i,1) = abs(int_f - trap_f);
    error_f(i,2) = abs(int_f - simp_f);
    
    trap_g = trapez(g, 0, 1, 10*i);
    simp_g = simpson(g, 0, 1, 10*i);
    error_g(i,1) = abs(int_g - trap_g);
    error_g(i,2) = abs(int_g - simp_g);
end

figure(1)
loglog(m, error_f)
grid on;
title('f= cos(3*pi*x)^4');
xlabel('Varying values for m');
ylabel('Quadrature Errors Using Simpson''s Rule and Trapezium Rule');
legend('Trapezium Rule', 'Simpson''s Rule');

figure(2)
loglog(m, error_g)
grid on;
title('g = sqrt(x)');
xlabel('Varying values for m');
ylabel('Quadrature Errors Using Simpson''s Rule and Trapezium Rule');
legend('Trapezium Rule', 'Simpson''s Rule');

lef = log(error_f);
leg = log(error_g);
lm = [ones(n,1) log(m)];

logm = lm'*lm;
solTrapf = lm'*lef(:,1);
solSimpf = lm'*lef(:,2);
solTrapg = lm'*leg(:,1);
solSimpg = lm'*leg(:,2);

k_trap_f = inv(logm)*solTrapf
k_simp_f = inv(logm)*solSimpf
k_trap_g = inv(logm)*solTrapg
k_simp_g = inv(logm)*solSimpg

figure(3)
loglog(msimp(1,:),msimp(1,:).^kb(1)*exp(kb(2)))
grid on;
ylabel('Error')'
xlabel('Varying m');
title('Least Squares Estimate of the Simpson''s Quadrature Error for f=sqrt(x)');
legend('Quadrature Error' ,'Error = Cm^k')

function I = trapez(f,a,b,m)
    xx = linspace(a,b,m);
    h = (b-a)/m;
    y = h*.5*f(xx(1)) + h*.5*f(xx(m));
    for i = 2:m-1
        y = y + h*f(xx(i));
    end
    I = y;
end

function I = simpson(f,a,b,m)
    xx = linspace(a,b,2*m);
    h = (b-a)/(2*m);
    y = h/3*(f(xx(1)) + 4*f(xx(2*m-1))+ f(xx(2*m)));
 
    for i = 1:m-1
        y = y + h/3*(4*f(xx(2*i))+2*f(xx(2*i)));
    end
    I = y;
end

%% 3a
x = linspace(0,1);

phi0 = @(x)x.^0;
phi1 = @(x)x - (1/2);
phi2 = @(x)x.^2 - x +(1/6);

phi0num3 = @(x)(x.^3).*phi0(x);
innerProd0 = @(x)phi0(x).*phi0(x);
phi1num3 = @(x)(x.^3).*phi1(x);
innerProd1 = @(x)phi1(x).*phi1(x);
phi2num3 = @(x)(x.^3).*phi2(x);
innerProd2 = @(x)phi2(x).*phi2(x);

phi30 = integral(phi0num3, 0, 1)/integral(innerProd0, 0, 1)
phi31 = integral(phi1num3, 0, 1)/integral(innerProd1, 0, 1)
phi32 = integral(phi2num3, 0, 1)/integral(innerProd2, 0, 1)

phi3 = @(x)x.^3 - phi30.*phi0(x) - phi31.*phi1(x) - phi32.*phi2(x);

%% 3b
c = polyfit(x, phi3(x), 3)
r = roots(c)
L0 = @(x)((x - r(2))/(r(1) - r(2))).*((x - r(3))/(r(1) - r(3)));
L1 = @(x)((x - r(1))/(r(2) - r(1))).*((x - r(3))/(r(2) - r(3)));
L2 = @(x)((x - r(1))/(r(3) - r(1))).*((x - r(2))/(r(3) - r(2)));
L0sq = @(x)L0(x).*L0(x);
L1sq = @(x)L1(x).*L1(x);
L2sq = @(x)L2(x).*L2(x);
W = [integral(L0sq, 0, 1) integral(L1sq, 0, 1) integral(L2sq, 0, 1)]'

approx = W(1).*r(1) + W(2).*r(2) + W(3).*r(3)

%% 3c
for k = 0:7
    f = @(x) power(x, k);
    GaussQuad = (W(1) * f(r(1)) +  (W(2) * f(r(2))) + (W(3) * f(r(3))));
    Simp = simpson(f, 0, 1, 1000);
    errorsGQ(k + 1) =  abs((1/(k+1)) - GaussQuad);
    errorsSimp(k + 1) = abs((1/(k+1)) - Simp);
end
 
for i = 0 : 7
    plot(i, errorsGQ(i+1),'rx')
    hold on
    plot(i, errorsSimp(i+1),'bo')
    hold on
end
 
title('Errors Using Gaussian Quadrature formula and Simpson Rule');
xlabel('Varying values for k');
ylabel('Errors');
legend('Gaussian Quadrature Errors', 'Simpson''s Rule Errors');
grid on;

function I = simpson(f,a,b,m)
xx = linspace(a,b,2*m);
h = (b-a)/(2*m);
y = h/3*(f(xx(1)) + 4*f(xx(2*m-1))+ f(xx(2*m)));
    for i = 1:m-1
        y = y + h/3*(4*f(xx(2*i))+2*f(xx(2*i)));
    end
    I = y;
end

%% 4a
w = @(x) exp(-x);

l0 = @(x)x.^0;

l0num1 = @(x)w(x).*l0(x).*x;
l0num2 = @(x)w(x).*l0(x).*x.^2;
l0num3 = @(x)w(x).*l0(x).*x.^3;
innerProdl0 = @(x) l0(x).*l0(x).*w(x);

l11 = integral(l0num1, 0, inf)/integral(innerProdl0, 0, inf);
l1 = @(x)x - l11.*l0(x);

l1num2 = @(x)w(x).*l1(x).*x.^2;
l1num3 = @(x)w(x).*l1(x).*x.^3;
innerProdl1 = @(x)l1(x).*l1(x).*w(x);

l21 = integral(l0num2, 0, inf)/integral(innerProdl0, 0, inf);
l22 = integral(l1num2, 0, inf)/integral(innerProdl1, 0, inf);
l2 = @(x)x.^2 - l21.*l0(x) - l22.*l1(x);

l2num3 = @(x)w(x).*l2(x).*x.^3;
innerProdl2 = @(x)l2(x).*l2(x).*w(x);

l31 = integral(l0num3, 0, inf)/integral(innerProdl0, 0, inf);
l32 = integral(l1num3, 0, inf)/integral(innerProdl1, 0, inf);
l33 = integral(l2num3, 0, inf)/integral(innerProdl2, 0, inf);
l3 = @(x)x.^3 - l31.*l0(x) - l32.*l1(x) - l33.*l2(x);

x = linspace(0, 10);
c0 = polyfit(x,l0(x),0)
c1 = polyfit(x,l1(x),1)
c2 = polyfit(x,l2(x),2)
c3 = polyfit(x,l3(x),3)
plot(x, l0(x), 'r', x, l1(x), 'g', x, l2(x), 'b', x, l3(x), 'k')
legend('l0 = 1', 'l1 = x - 1', 'l2 = x^2 - 4x +2', 'l3 = x^3 -9x^2 + 18x - 6')

%% 4b
c2 = polyfit(x, l2(x), 2)
r2 = roots(c2)
c3 = polyfit(x, l3(x), 3)
r3 = roots(c3)

%% 4c
f = @(x)exp(-x);
g = @(x)exp(-(x.^2) + x);
approxf2 = .853553*f(.585786) + .146447*f(3.41421)
approxf3 = 0.711093*f(0.415775) + 0.278518*f(2.29428) + 0.0103893*f(6.28995)
approxf4 = 0.603154*f(0.322548) + 0.357419*f(1.74576) + 0.0388879*f(4.53662) + 0.000539295*f(9.39507)
approxg2 = .853553*g(.585786) + .146447*g(3.41421)
approxg3 = 0.711093*g(0.415775) + 0.278518*g(2.29428) + 0.0103893*g(6.28995)
approxg4 = 0.603154*g(0.322548) + 0.357419*g(1.74576) + 0.0388879*g(4.53662) + 0.000539295*g(9.39507)
errorf2 = (1/2) - approxf2;
errorg2 = (sqrt(pi)/2) - approxg2;
errorf3 = (1/2) - approxf3;
errorg3 = (sqrt(pi)/2) - approxg3;
errorf4 = (1/2) - approxf4;
errorg4 = (sqrt(pi)/2) - approxg4;
errorf = [errorf2 errorf3 errorf4]