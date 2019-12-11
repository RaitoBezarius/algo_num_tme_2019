% Q1.
N = 1000;
subplot(2,1,1);

m_timing = zeros(N, 1);
blas_timing = zeros(N, 1);
for n=1:N
    A = rand(n,n);
    
    tic;
    sqrt(manual_norm(A));
    m_timing(n) = toc;
    
    tic;
    norm(A, 1);
    blas_timing(n) = toc;
    
    fprintf("norm: %d/%d\n", n, N);
end

hold on;
plot(1:N, m_timing);
plot(1:N, blas_timing);
legend("MATLAB maison", "BLAS");
title("Comparaison des BLAS et de MATLAB maison (norme)");
ylabel("Temps en secondes");
xlabel("Taille de la matrice");

% Q2.
N = 350;
subplot(2,1,2);
m_timing = zeros(N, 1);
blas_timing = zeros(N, 1);
for n=1:N
    A = rand(n,n);
    B = rand(n,n);
    
    tic;
    mat_product(A,B);
    m_timing(n) = toc;
    
    tic;
    A*B;
    blas_timing(n) = toc;
    fprintf("matrix product: %d/%d\n", n, N);
end

hold on;
plot(1:N, m_timing);
plot(1:N, blas_timing);
legend("MATLAB maison", "BLAS");
title("Comparaison des BLAS et de MATLAB maison (multiplication matricielle)");
ylabel("Temps en secondes");
xlabel("Taille de la matrice");


function s = manual_norm(M)
    [m, n] = size(M);
    s = zeros(1, m);
    for i=1:m
        for j=1:n
            s(i) = s(i) + abs(M(i, j));
        end
    end
end

function R = mat_product(A, B)
    [m, n] = size(A);
    [p, q] = size(B);
    
    if p ~= n
        error('format invalid');
    end
    
    R = zeros([m, q]);
    
    for i=1:m
        for j=1:q
            for k=1:n
                R(i, j) = R(i, j) + A(i, k)*B(k, j);
            end
        end
    end
end    