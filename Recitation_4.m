A = round(10*rand(5));
norm(A)
norm(A,1)
norm(A, 2)
norm(A, Inf)

Q = A*A';
cond(Q)
eig(Q)
max(eig(Q))/min(eig(Q))

A1 = [1,1;1,1.001]
b = [2;2]
K = cond(A1)
x = A1\b

ee = 1e-4
bb = b + [0; ee]
xx = A1\bb
dx = xx-x
db = bb-b

norm(dx)/norm(x)
norm(db)/norm(b)