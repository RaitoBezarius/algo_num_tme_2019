% (0 + 2X + 3X^2 + 5X^3)(0 + 1X + 2X^2 + 3X^3)
MAX_N = 3000;
fftprod_timing = zeros(MAX_N, 1);
naive_timing = zeros(MAX_N, 1);
conv_timing = zeros(MAX_N, 1);

for n=1:MAX_N
    P = rand(1, n);
    Q = rand(1, n);
    
    tic;
    fft_prod(P, Q);
    fftprod_timing(n) = toc;
    
    %tic;
    %naive_product(P, Q);
    %naive_timing(n) = toc;
    
    tic;
    conv(P, Q);
    conv_timing(n) = toc;
    
    fprintf("%d/%d\n", n, MAX_N);
end

R = 1:MAX_N;
hold on;
plot(R, fftprod_timing);
%plot(R, naive_timing);
plot(R, conv_timing);
legend('Produit par FFT', 'Naïf', 'Convolution MATLAB');
ylabel('Temps en secondes');
xlabel('Degré maximal des polynômes');
title('Efficacité des méthodes de multiplication polynomiales');


function res = fft_prod(P, Q)
    N = 2*max(length(P), length(Q));
    P_ = fft([P zeros(1, N - length(P))]);
    Q_ = fft([Q zeros(1, N - length(Q))]);
    
    R_ = P_ .* Q_;
    
    res = ifft(R_);
end

function res = naive_product(p,q)
    n = length(p) ;
    R = zeros(1,2*n-1) ;
    for i=1:n
        for j=1:n
            R(i+j-1) = R(i+j-1) + p(i)*q(j) ;
        end
    end
    res = R;
end