function [ Q,R ] = gramschmidt( A, flag )
%gramshmidt - This function computes the QR Factorisation
%             of the matrix A (full rank) using the Classical
%             Gram-Schmidt(flag=true) and Modified Gram-Schmidt(flag=flase)

%   Classical Gram-Schmidt is a numerically unstable algorithm
%   Whereas it stable counter part (arithmetically equivalent) 
%   algorithm, the modified Gram-Schmidt is numerically stable!

r = size(A);
n = r(2);
Q = zeros(r);
R = zeros(n);

% Classical Gram-Schmidt
if flag
    for j=1:n
        w = A(:,j);
        for i=1:j-1
            R(i,j) = Q(:,i)'*A(:,j);
            w = w - R(i,j)*Q(:,i);
        end
        R(j,j) = norm(w,2);
        Q(:,j) = w/R(j,j);
    end

% Modified Gram-Schmidt
else
    for j=1:n
        R(j,j) = norm(A(:,j),2);
        Q(:,j) = A(:,j)/R(j,j);
        for i = j+1:n
            R(j,i) = Q(:,j)'*A(:,i);
            A(:,i) = A(:,i) - R(j,i)*Q(:,j);
        end
    end
end



