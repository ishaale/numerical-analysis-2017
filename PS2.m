%% 6a

A = [6,2,1,-1;2,4,1,0;1,1,4,-1;-1,0,-1,3];

L = zeros(4, 4);
U = zeros(4,4);

 for i=1:4
     for k=1:i-1
         L(i,k)=A(i,k);
         for j=1:k-1
             L(i,k)= L(i,k)-L(i,j)*U(j,k);
         end
         L(i,k) = L(i,k)/U(k,k);
     end
     for k=i:4
         U(i,k) = A(i,k);
         for j=1:i-1
             U(i,k)= U(i,k)-L(i,j)*U(j,k);
         end
     end
 end
 for i = 1:4
     L(i, i) = 1;
 end
 
 L
 U
 
 %% 6b
 
A = [];
n = length(A);
 
L = zeros(n, n);
U = zeros(n, n);

 for i=1:n
     for k=1:i-1
         L(i,k)=A(i,k);
         for j=1:k-1
             L(i,k)= L(i,k)-L(i,j)*U(j,k);
         end
         if abs(U(k, k)) < 10e-8
             fprintf('Error: Leading Principle Submatrix is singular')
             break
         else
            L(i,k) = L(i,k)/U(k,k);
         end
     end
     for k=i:n
         U(i,k) = A(i,k);
         for j=1:i-1
             U(i,k)= U(i,k)-L(i,j)*U(j,k);
         end
     end
 end
 for i = 1:n
     L(i, i) = 1;
 end
 
 %% TA Answer
 
 function [L, U] = myLU(A) %does not work
 
 n = size(A);
 L = eye(n);
 U = zeros(n);
 U(1,:) = A(1,:);
 
    for i=1:n
       for k=1:i-1
           L(i,k)=A(i,k);
           for j=1:k-1
               L(i,k)= L(i,k)-L(i,j)*U(j,k);
           end
           if abs(U(k, k)) < 10e-8
               fprintf('Error: Leading Principle Submatrix is singular')
               break
           else
              L(i,k) = L(i,k)/U(k,k);
           end
       end
       for k=i:n
           U(i,k) = A(i,k);
           for j=1:i-1
               U(i,k)= U(i,k)-L(i,j)*U(j,k);
           end
       end
    end
    for i = 1:n
        L(i, i) = 1;
    end
 end