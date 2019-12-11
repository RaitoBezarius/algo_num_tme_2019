itermax = 10000;

[iters, B] = inv_matrix_newton(vander(1:20), 0.1, itermax);

disp(norm(vander(1:20), 'fro')^2);

plot(0:length(iters)-1, iters);
legend('InvNewton');
xlabel('Itérations');
ylabel('Erreur en norme (e_n)');
title('Inversion de matrice par Newton');

% Q2.
function [iters, B] = inv_matrix_newton(A, tol, itmax)
    B = A' / (norm(A, 'fro')^2); % A^T / Tr(A^T A) = A^T / ||A||^2_Frobenius
    iters = [iteration_error(B, A)];
    i = 0;
    while iteration_error(B, A) > tol && i < itmax
        B = 2*B - B*A*B;
        iters = [iters iteration_error(B, A)];
        i = i + 1;
    end
end

function e = iteration_error(X_n, A)
    n = length(A);
    e = norm(eye(n) - X_n*A);
end