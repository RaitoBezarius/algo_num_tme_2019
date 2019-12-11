f = cell(1, 5);
f{1} = @(x)(0.5 - x.*exp(-x.^2));
f{2} = @sin;
f{3} = @(x)atan(x)^2;
f{4} = @(x)abs(log(x));
f{5} = @abs;

g = cell(1, 5);
g{1} = "x → 0.5 - xexp(-x^2)";
g{2} = "sin";
g{3} = "arctan(x)^2";
g{4} = "|log(x)|";
g{5} = "|⋅|";

ranges = [0 2; 
    0 pi/2;
    -1 1;
    1/2 4;
    -1 1];

for i = 1:5
    disp(f{i})
    x_goldensearch = golden_search(ranges(i, 1), ranges(i, 2), f{i}, 0.0001, false);
    x_newton = newton((ranges(i, 1) + ranges(i, 2))/2, f{i}, 100, 0.00001);
    x_minbnd = fminbnd(f{i}, range(i, 1), ranges(i, 2));
    fprintf("golden search=%f\n", golden_search(ranges(i, 1), ranges(i, 2), f{i}, 0.0001, false));
    fprintf("newton=%f\n", x_newton);
    fprintf("fminbd=%f\n", x_minbnd);
    figure;
    Y = [x_goldensearch x_newton x_minbnd];
    X = categorical({'Section dorée', 'Newton', 'MATLAB'});
    X = reordercats(X, {'Section dorée', 'Newton', 'MATLAB'});
    bar(X, Y);
    ylabel('x^* optimal trouvé');
    xlabel('Méthodes');
    title("Optimization de la fonction: " + g{i});
end

function res = num_deriv(f, x, h)
    res = (f(x + h) - f(x))/h;
end

function res = num_deriv2(f, x, h)
    res = (f(x + h) - 2*f(x) + f(x - h))/h^2;
end

function res = newton(x0, f, nb_iterations, step_size)
    x_k = x0;
    
    df = num_deriv(f, x_k, step_size);
    df2 = num_deriv2(f, x_k, step_size);
    
    for i=1:nb_iterations
        x_k = x_k - df/df2;
        df = num_deriv(f, x_k, step_size);
        df2 = num_deriv2(f, x_k, step_size);
    end
    
    res = x_k;
end

function res = golden_search(a, b, f, tol, debug)
    tau = (sqrt(5) - 1)/2;
    x1 = a + (1 - tau)*(b - a);
    f1 = f(x1);
    x2 = a + tau*(b - a);
    f2 = f(x2);
    i = 0;
    
    while (b - a) > tol && i < 100
        if f1 > f2
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + tau*(b - a);
            f2 = f(x2);
        else
            b = x2;
            x2 = x1;
            f2 = f1;
            x1 = a + (1 - tau)*(b - a);
            f1 = f(x1);
        end
        i = i + 1;
    end
    
    if debug
        fprintf("golden search debug info\n");
        fprintf("f1 = %f, f2 = %f\n", f1, f2);
        fprintf("x1 = %f, x2 = %f\n", x1, x2);
        fprintf("nb iterations = %i\n", i);
        fprintf("error = %f\n", (b - a));
    end
    res = x1;
end
