MAX_N = 3000;
MIN_N = 100;
STEP = 100;
w0 = 0.4;
tol = 10e-10;

%figure(1);
%plot_performance_by_generation("I_n", @identity_matrix, MIN_N, MAX_N, STEP, w0, tol);
%figure(2);
%plot_performance_by_generation("Matrices à diagonale dominante", @generate_nice_matrix_for_jacobi, MIN_N, MAX_N, STEP, w0, tol);
%plot_performance_by_generation("Matrice très mal conditionnées", @generate_ill_conditionned_matrix, MIN_N, MAX_N, STEP, w0, tol);
plot_performance_by_generation("Matrices définies positives", @generate_sdp, MIN_N, MAX_N, STEP, w0, tol);

function plot_performance_by_generation(test_name, g, start, max_N, step, w0, tol)
    J_ = zeros((max_N - start)/step + 1, 2);
    GS_ = zeros((max_N - start)/step + 1, 2);
    SOR_ = zeros((max_N - start)/step + 1, 2);
    GCJ_ = zeros((max_N - start)/step + 1, 2);
    
    R = start:step:max_N;
    
    fprintf("Computing all the values...\n");
    for q = R
        [A, b] = g(q);
        x0 = zeros(q,1);
        %fprintf("Kappa(A) = %f\n", cond(A));
        [jacobi, gs, sor, gcj] = test_algorithms(A, b, x0, w0, tol);
        j = ((q - start)/step) + 1;
        J_(j,:) = jacobi;
        GS_(j,:) = gs;
        SOR_(j,:) = sor;
        GCJ_(j,:) = gcj;
        fprintf("%d/%d computed.\n", q, max_N);
    end
    
    disp(J_(:,1));
    disp(GS_(:,1));
    %disp(SOR_(:,1));
    disp(GCJ_(:,1));
    
    fprintf("Algorithms computed. Drawing now.\n")
    figure
    tiledlayout(2,1)
    
    for k = 1:2
        nexttile
        hold on
        plot(R, J_(:,k));
        plot(R, GS_(:,k));
        plot(R, SOR_(:,k));
        %plot(R, GCJ_(:,k));
        lgd = legend('Jacobi', 'Gauss Seidel', 'SOR', 'y=log(log(log(n)))', 'Conjugate gradient');
        xlabel('Taille de la matrice');
        if k == 1
            title(lgd, test_name + " - norme de l'erreur en fonction de la taille de la matrice");
            ylabel("Norme 2 de l'erreur entre la solution et l'approximation");
        else
            title(lgd, test_name + " - temps d'exécution en fonction de la taille de la matrice");
            ylabel('Temps en secondes');
        end
    end
end 



function [jacobi, gs, sor, gcj] = test_algorithms(A, b, x0, w0, tol)
    tic
    x_jacobi_approx = Jacobi(A, b, x0, tol);
    jacobi_time = toc;
    jacobi = [compute_approx_error(x_jacobi_approx) jacobi_time];
    
    tic
    x_gs_approx = GS(A, b, x0, tol);
    gs_time = toc;
    gs = [compute_approx_error(x_gs_approx) gs_time];
    
    tic
    x_sor_approx = SOR(A, b, x0, w0, tol);
    sor_time = toc;
    sor = [compute_approx_error(x_sor_approx) sor_time];
    
    tic
    x_gcj_approx = GCJ(A, b, x0);
    gcj_time = toc;
    
    gcj = [compute_approx_error(x_gcj_approx) gcj_time];
end

function e = compute_approx_error(x_approx)
    n = length(x_approx);
    e = norm(ones(n, 1) - x_approx);
end



% Matrices très mal conditionnés
function [A, b] = generate_ill_conditionned_matrix(n)
    A = vander(1:n); % Vandermonde([[1, n]])
    x = ones(n, 1);
    b = A*x;
end

% Diagonally dominant matrix (proven convergence for Jacobi)
function [A, b] = generate_nice_matrix_for_jacobi(n)
    x = 0.05; % part de zéros hors de la diagonale
    k = round(n*(n-1)*x);

    data = randn(n*(n-1)-k,1);
    data = [data;zeros(k,1)];
    data = data(randperm(length(data)));

    diag_index = 1:n+1:n*n;
    offd_index = setdiff(1:n*n,diag_index);
    A = zeros(n,n);
    A(offd_index) = data;
    A(diag_index) = sum(A,1);
    
    b = A*ones(n, 1);
end

% SDP matrices
function [A, b] = generate_sdp(n)
    A = rand(n, n);
    A = (1/2)*(A + A');
    A = A + n*eye(n);
    b = A*ones(n,1);
end

% Identity matrices
function [A, b] = identity_matrix(n)
    A = eye(n);
    b = A*ones(n, 1);
end

function x = Jacobi(A, y, x0, tol)
    n = length(y);
    x = x0;
    xold = x;
    itmax = 5000;
    divcounter = 1;
    it = 1;
    while (divcounter <= itmax && norm(xold - x) > tol && ~isnan(x)) || it == 1
        for i=1:n
            x(i) = (y(i) - A(i,[1:i - 1, i + 1:n]) * xold([1:i-1, i+1:n]))/A(i,i);
        end
        if norm(xold - x) < tol
            divcounter = 0;
        end
        xold = x;
        divcounter = divcounter + 1;
        it = it + 1;
    end
end

function x = GS(A, y, x0, tol)
    n = length(y);
    x = x0;
    
    itmax = 5000;
    divcounter = 1;
    it = 1;
    
    xold = x;
    
    while (divcounter <= itmax && norm(xold - x) > tol && ~isnan(x)) || it == 1
        for i=1:n
            x(i) = (y(i) - A(i, 1:i -1)*x(1:i-1) - A(i, i + 1:n) * x(i + 1:n)) / A(i, i);
        end
        if norm(xold - x) < tol
            divcounter = 0;
        end
        xold = x;
        divcounter = divcounter + 1;
        it = it + 1;
    end
end


function x = SOR(A, y, x0, w, tol)
    n = length(y);
    x = x0;
    
    itmax = 5000;
    divcounter = 1;
    it = 1;
    
    xold = x;
    
    while (divcounter <= itmax && norm(xold - x) > tol && ~isnan(x)) || it > 1
        for i=1:n
            x(i) = (w/A(i,i))*(y(i) - A(i, 1:i-1)*x(1:i-1) - A(i, i+1:n).x(i + 1:n)) + (1 - w)*x(i);
        end
        
        if norm(xold - x) < tol
            divcounter = 0;
        end
        xold = x;
        divcounter = divcounter + 1;
        it = it + 1;
    end
end

function x = GCJ(A, y, x0)
    n = length(y);
    r = y - A*x0;
    p = r;    
    tol = 1e-10;
    
    rsold = r' * r;
    x = x0;
    for k=0:n-1
        Ap = A*p;
        a = rsold / (p' * Ap);
        x = x + a*p;
        r = r - a*Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < tol
            break;
        end
        p = r + (rsnew / rsold) * p;
        if isnan(x)
            break;
        end
    end
end
