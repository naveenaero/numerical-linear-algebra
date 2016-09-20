kappa = [1, 1e3, 1e6, 1e9];
error = zeros(3,4);
Qnorm = zeros(3,4);
for i=length(kappa)
    A = gallery('randsvd',100, kappa(i));
    [Q1,R1] = cgsmgs(A, true);
    [Q2, R2] = cgsmgs(A, false);
    [Q3,R3] = qr(A);
    normA = norm(A);
    error(1,i) = norm(A - Q1*R1)/normA;
    error(2,i) = norm(A - Q2*R2)/normA;
    error(3,i) = norm(A - Q3*R3)/normA;
    Qnorm(:,i) = [norm(Q1), norm(Q2), norm(Q3)]';
end
