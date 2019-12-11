x = linspace(-1, 1);
y = Taylored(x);
y2 = BadFunction(x);
plot(x, y, 'k.', x, y2, '--')
xlim([-1 1])
grid on

function res = BadFunction(x)
    res = zeros(length(x), 1);
    for i=1:length(x)
        res(i) = (exp(x(i)) - 1 - x(i)) / x(i)^2;
    end
end

function res = Taylored(x)
    acc = 0;
    for i=1:15
        acc = acc + x.^(i - 1)/factorial(i+1);
    end
    res = acc;
end