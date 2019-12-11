% Q2. D'après WA, on attend x ~0.865474.
x_ = 0.865474;
test_f = @(x)(x^3 - cos(x));

[iters, x] = newton(test_f, 0, 0.01, 0.01);
fprintf("Newton solution: %f\n", x);

error = abs(iters - x_);
plot(0:(length(iters) - 1), error);
legend('x_k');
ylabel('Erreur en valeur absolue');
xlabel('Itération');
title('Itérations de la méthode de Newton appliqué au problème x^3 = cos(x)');

% Q1.
function df = approx_d(f, a, stepsize)
    df = (f(a + stepsize) - f(a)) / stepsize;
end

function [iters, x] = newton(f, x0, tol, eps)
    oldx = x0;
    x = oldx - f(oldx)/approx_d(f, oldx, eps); 
    iters = [x0 x];
    
    while abs(oldx - x) > tol
        oldx = x;
        x = oldx - f(oldx)/approx_d(f, oldx, eps);
        iters = [iters x];
    end
 end