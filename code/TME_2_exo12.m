syms f(x, y);
f(x, y) = 100*(y - x^2)^2 + (1 - x)^2;

% Q1.
G = gradient(f);
H = hessian(f);

% Q2.
fprintf("grad f(x^*) = (%g, %g), det H_f(x^*) = %i\n", G(1, 1), det(H(1, 1)));
% si grad f(x*) = 0 & det H_f(x*) > 0 alors x* est minimum local.

% Q3.
figure;
subplot(2,1,1);
iters = newton(f, [-1 -2], 5)
hold on
for i=1:5
    plot(iters(i, 1), iters(i, 2), 'ro');
end
for i=1:4
    p1 = iters(i, :);
    p2 = iters(i + 1,:);
    dp = p2 - p1;
    quiver(p1(1), p1(2), dp(1), dp(2), 0);
end
fcontour(f, [-1.5 2 -3 3]);
title("Contour de la fonction muni des itérations de Newton dans l'ordre");

% Q4.
error = zeros(5, 1);

for i=1:6
    error(i) = norm(iters(i,:) - [1 1]);
end

subplot(2,1,2);
plot(0:5, error);
ylabel('Erreur en norme 2');
xlabel('Itérations');
title('Erreurs de la méthode de Newton');


% Newton impl
function iters = newton(f, x0, itmax)
    H = hessian(f);
    G = gradient(f);
    x = x0;
    iters = zeros(itmax + 1, 2);
    iters(1,:) = x0;
    for i=1:itmax
        s = linsolve(H(x(1), x(2)), -G(x(1), x(2)))';
        x = x + s;
        iters(i + 1,:) = x;
    end
end