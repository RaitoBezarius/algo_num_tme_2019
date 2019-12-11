MAX_N = floor(log2(1000000)) + 1;

fft_rec_time = zeros(MAX_N, 1);
fft_iter_time = zeros(MAX_N, 1);
matlab_fft_time = zeros(MAX_N, 1);

for n=1:MAX_N
    P = rand(2^n,1);
    W = exp(-2*pi*1i/n);
    tic
    fft_rec(P, W);
    fft_rec_time(n) = toc;    
    tic;
    fft_iter(P);
    fft_iter_time(n) = toc;
    tic;
    fft(P);
    matlab_fft_time(n) = toc;
end

hold on;
range = arrayfun(@(n) 2^n, 1:MAX_N);
plot(range, fft_rec_time);
plot(range, fft_iter_time);
plot(range, matlab_fft_time);
legend('FFT récursive', 'FFT itérative', 'MATLAB FFT');
ylabel('Temps en secondes');
xlabel('Taille du vecteur');
title('Vitesse de différentes implémentations de la FFT');

function res = fft_rec(A, w)
    n = length(A);
    res = A;
    if w ~= 1 && n > 1
        B = fft_rec(A(1:2:n), w^2);
        C = fft_rec(A(2:2:n), w^2);
        
        res = zeros(n, 1);
        for j=1:n/2
            res(j) = B(j) + w^(j - 1) * C(j);
            res(j + n/2) = B(j) - w^(j - 1) * C(j);
        end
    end
end

function res = fft_iter(A)
    n = length(A);
    q = log2(n);
    res = A;
    
    if n > 1
        res = zeros(n, 1);
        [~, indexes] = bitrevorder(0:(n-1));
        for j=1:n
            res(indexes(j)) = A(j);
        end
        bsize = 1;
        
        for i=1:q
            bsize = bsize * 2;
            w = exp(1i*2*pi/bsize);

            for j=1:bsize:n
                k = j + bsize/2;
                W = 1;

                for p=0:(bsize/2 - 1)
                    c_even = res(j + p);
                    c_odd = res(k + p);

                    res(j + p) = c_even + W*c_odd;
                    res(k + p) = c_even - W*c_odd;
                end
            end
        end
    end
end

unction res = fft_iter(A)
    n = length(A);
    q = log2(n);
    res = A;
    
    if n > 1
        res = zeros(n, 1);
        [~, indexes] = bitrevorder(0:(n-1));
        for j=1:n
            res(indexes(j)) = A(j);
        end
        bsize = 1;
        
        for i=1:q
            bsize = bsize * 2;
            w = exp(1i*2*pi/bsize);

            for j=1:bsize:n
                k = j + bsize/2;
                W = 1;

                for p=0:(bsize/2 - 1)
                    c_even = res(j + p);
                    c_odd = res(k + p);

                    res(j + p) = c_even + W*c_odd;
                    res(k + p) = c_even - W*c_odd;
                end
            end
        end
    end
end

function z = rootsOfUnity(n)
 
    assert(n >= 1,'n >= 1');
    z = roots([1 zeros(1,n-1) -1]);
 
end
