%% Test Case for CGS, MGS and Matlab QR function
%% Defining Variables 
clear;
clc;
kappa = [1,1e3,1e6,1e9];
error = zeros(3,length(kappa));
I = eye(100);
%% Testing for different $$\kappa(A) $$ 
for i=1:length(kappa)
    A = gallery('randsvd',100, kappa(i));
    [Q1,R1] = gramschmidt(A, true);
    [Q2,R2] = gramschmidt(A, false);
    [Q3,R3] = qr(A);
    error(1,i) = norm(Q1'*Q1 - I);
    error(2,i) = norm(Q2'*Q2 - I);
    error(3,i) = norm(Q3'*Q3 - I);
end

%% Plotting the relative error norms
loglog( kappa,error(1,:), ...
        kappa,error(2,:), ...
        kappa,error(3,:));
xlabel('kappa');
ylabel('Relative l2 error norm');
legend('Classical GS', 'Modified GS', 'Matlab QR');
title('Comparison of Classical GS, Modified GS and Matlab QR factorisation');
