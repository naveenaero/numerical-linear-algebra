function [ Q,R ] = gramschmidt( A, flag )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

r = size(A);
n = r(2);
Q = zeros(r);
R = zeros(n);

if flag
    for j=1:n
        w = A(:,j);
        for i=1:j-1
            R(i,j) = Q(:,i)'*A(:,j);
        end
        for i=1:j-1
            w = w - R(i,j)*Q(:,i);
        end
        R(j,j) = norm(w);
        Q(:,j) = w/R(j,j);
    end
    
else
    for j=1:n
        w = A(:,j);
        for i=1:j-1
            R(i,j) = Q(:,i)'*A(:,j);
            w = w - R(i,j)*Q(:,i);
        end
        R(j,j) = norm(w);
        Q(:,j) = w/R(j,j);
    end

end



