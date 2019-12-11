% Q1
syms f(x)
x = sym('x', [1 2]);
f(x) = x(1)^2 - 2*x(2)^2;

x0 = [-1 -1];

x_ = [0 0];

subplot(2,1,1);
g_range = 0.1:0.2:1.3;
g_f_errors = zeros(length(g_range), 1);
for n=1:numel(g_range)
    alpha = g_range(n);
    y = gradient_method(f, alpha, x0, 35);
    g_f_errors(n) = norm(y - x_);
    fprintf("[Initial function] With alpha = %f, gradient method: [%0.5e %0.5e] (error: %0.5e)\n", alpha, y, norm(y - x_));
end
semilogy(g_range, g_f_errors);
legend('(x, y) → x^2 - 2y^2');
title("Erreur de la méthode du gradient en fonction d'alpha");
ylabel('Erreur en norme 2');
xlabel('Alpha');

% Q2
syms rosenbrock(x)
x = sym('x', [1 2]);
rosenbrock(x) = 100*(x(2) - x(1)^2)^2 + (1 - x(1))^2;
x0 = [-1 1.2];

x_ = [1 1];

rosenbrock_range = 0.0001:0.0004:0.0043;
rosenbrock_g_errors = zeros(length(rosenbrock_range), 1);
for n=1:numel(rosenbrock_range)
    alpha = rosenbrock_range(n);
    y = gradient_method(rosenbrock, alpha, x0, 10);
    fprintf("[Rosenbrock] With alpha = %f, gradient method: [%0.5e %0.5e] (error: %0.5e)\n", alpha, y, norm(y - x_));
    rosenbrock_g_errors(n) = norm(y - x_);
end

subplot(2,1,2);
plot(rosenbrock_range, rosenbrock_g_errors);
legend('Fonction de Rosenbrock');
xlabel('Alpha');
ylabel('Erreur en norme 2');
title("Erreur de la méthode du gradient en fonction d'alpha");

% Q3
y = gradient_with_wolfe(rosenbrock, x0, 12);
fprintf("[Rosenbrock] Wolfe's Gradient method: [%0.5e %0.5e] (error: %0.5e)\n", y, norm(y - x_));


function y = gradient_method(f, alpha, x0, itmax)
    y = x0;
    G = gradient(f);
    for i=1:itmax
        y = y - alpha*G(y(1), y(2))';
    end
end

function t = wolfe_linear_search(f, G, y, d)
    t_g = 0;
    t_d = +inf;
    t = 1;
    m_1 = 0.1;
    m_2 = 0.9;
    d_g0 = dot(-d, d);
    
    
    go = true;
    while go
        % dérivées directionnelles
        % h : R → R^p
        % h(t) = x + th
        % g = f rond h
        % g'(t) = (grad f(x + th) | h)
        p = y + t*d;
        d_gt = dot(G(p(1), p(2))', d);
        
        c1 = m_1 * t * d_g0;
        c2 = m_2 * d_g0;
        
        % g(t) <= g(0) + m_1 t g'(0) && g'(t) >= m_2 g'(0)
        if f(p(1), p(2)) <= f(y(1), y(2)) + c1 && d_gt >= c2
            go = false;
        % g(t) >= g(0) + m_1 t g'(0)
        elseif f(p(1), p(2)) > f(y(1), y(2)) + c1
            t_d = t;

            if isinf(t_d)
                t = 10*t_g;
            else
                t = (t_d + t_g) / 2;
            end
        % g(t) <= g(t) + m_1 t g'(0) && g'(t) < m_2 g'(0)
        elseif f(p(1), p(2)) <= f(y(1), y(2)) + c1 && d_gt < c2
            t_g = t;
            if isinf(t_d)
                t = 10*t_g;
            else
                t = (t_d + t_g)/2;
            end
        end
    end
end

function y = gradient_with_wolfe(f, x0, itmax)
    y = x0;
    G = gradient(f);
    
    for i=1:itmax
        fprintf("[Wolfe] %d/%d\n", i, itmax);
        u = G(y(1), y(2))';
        a = wolfe_linear_search(f, G, y, -u);
        y = y - a*u;
    end
end