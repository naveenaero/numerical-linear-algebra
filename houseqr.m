function [ Q,R,U ] = houseqr( A )
%HOUSEQR This function computes the QR factorisation of 
%        input matrix A where Q is orthonormal and R is 
%        upper triangular   

[m,n]=size(A);
U = zeros(m,n);
Q = eye(m);
%% Function for computation of reflection vector v and coefficient $\beta$
    function [v, beta] = reflect(x)
        v = x;
        sigma = sign(x(1))*sqrt(x'*x);
        v(1) = x(1) + sigma;
        beta = 1/(sigma*v(1));
    end
%% Computation of the Q, R, U matrices
for j = 1:n
    [v,beta] = reflect(A(j:m,j));
    U(j:end,j) = v;
    temp = v'*A(j:m,j:n);
    A(j:m,j:n) = A(j:m,j:n) - beta*v*temp;
    H = eye(m);
    H(j:m,j:m) = H(j:m,j:m) - beta*(v*v');
    Q = Q*H';
end
R = A;
end

