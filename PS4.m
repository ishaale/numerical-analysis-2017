%% 2_a,b
A = [-2 1 4;1 1 1;4 1 -2];
x0 = [1 2 -1]';
v0 = [1 2 1]';

X1 = power_method(A,x0)
X2 = power_method(A,v0)
eig(A);

function vecSeq = power_method(A, x)
    vecSeq = [];
    eigSeq = [];
    
    n = length(A);
    %creating initial vector x0
    for i = 1:n
        vecSeq(i) = x(i,1);
    end
    %power method
    for i = 1:10
        for j = 1:n 
            y = A*vecSeq(i,:)';
            vecSeq(i+1,:) = y/(norm(y));
        end
        eigSeq(i+1,:) = vecSeq(i+1,:) * A * vecSeq(i+1,:)';
    end
    eig = eigSeq
end

%%2_c,d
A = [-2 1 4;1 1 1;4 1 -2];
x = [1 2 -1]';

X1 = inverse_power(A, x, 2)
X2 = inverse_power(A, x, -1)
X3 = inverse_power(A, x, -6.1)

function vecSeq = inverse_power(A, x, theta)
    vecSeq = [];
    eigSeq = [];
    
    n = length(A);
    I = eye(n);
    %creating initial vector
    for i = 1:n
        vecSeq(i) = x(i,1);
    end
    %Inverse Iteration Method
    for i = 1:20
        for j = 1:n
            y = inv(A - theta*I)*vecSeq(i,:)';
            vecSeq(i+1,:) = y/(norm(y));
        end
        eigSeq(i+1,:) = vecSeq(i+1,:) * A * vecSeq(i+1,:)';
    end
    eig = eigSeq
end

%%3_a,b
A = [2 1 2 2;1 -7 6 5;2 6 2 -5;2 5 -5 1];
T = h_trid(A)
function T = h_trid(A)

n = length(A);  % Preallocations. 
v = zeros(n,1);  
I = eye(n);  
A = A;  

for j=1:n-2  % Build each vector j and run the whole procedure.
    v(1:j) = 0;
    S = sqrt(sum(A(j+1:end,j).^2));
    v(j+1) = sqrt(.5*(1+abs(A(j+1,j))/(S+2*eps)));
    v(j+2:n) = A(j+2:n,j)*sign(A(j+1,j))...
                   /(2*v(j+1)*S+2*eps);
    H = I-2*v*v';
    T = H*A*H;
    A = T;
end
end

%%4b
n = 15;
r = 1:n;
c = poly(r);
c1 = c((n+1):-1:2)
A = [[zeros(1,n-1); eye(n-1)] -c1']

eig(A)
roots(c)

%%5a,b
n = 350;
A = (max(2,randn(n,n))-2);
A = A - diag(diag(A));
L = A * diag(1./(max(1e-10,sum(A,1))));

spy(L)
eig(L);

plot(real(eig(L)), imag(eig(L)), '*')
hold on;

n = 100;
angle = 0:2*pi/n:2*pi;            
R = 1;                   
x = R*cos(angle);  y = R*sin(angle);    

plot(x,y);                             
axis equal;
grid on;

%%5c
n = 350;
A = (max(2,randn(n,n))-2);
A = A - diag(diag(A));
L = A * diag(1./(max(1e-10,sum(A,1))));

k = .25;
E = zeros(size(A));
for i = 1:size(A)
    for j = 1:size(A);
        E(i,j) = 1/n;
    end
end

for k = .25:.5:2.25
    S = k*L+(1-k)*E;
    eig(S);
    figure();
    plot(real(eig(S)), imag(eig(S)), 'o')
    hold on;
end

