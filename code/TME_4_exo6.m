load clown.mat;
% valeur intéressantes: eps=5*10e1, eps=5*10e2, eps=5*10e3.

for eps=[2*10e2]
    disp(eps);
    figure;
    colormap('gray');
    hold on;

    subplot(1, 2, 1);
    Y=compresser_image(X, eps);
    spy(Y);
    subplot(1, 2, 2);
    image(ifft2(Y));
end

function res = compresser_image(X, eps)
    Y = fft2(X);
    for i=1:size(Y,1)
        for j=1:size(Y,2)
            if norm(Y(i,j)) <= eps
                Y(i,j) = 0;
            end
        end
    end
    fprintf("Taux de compression: %f %% (%d)\n", taux_compression(X, Y)*100, nnz(Y));
    res = Y;
end

function taux = taux_compression(X_orig, X_compresse)
    taux = 1 - nnz(X_compresse) / nnz(X_orig);
end