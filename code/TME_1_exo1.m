format longE;
z = Higham(4)

x = logspace(0, 1, 2013);
y = Higham(x);
plot(x, y, 'k.', x, x, '--')

function res = Higham(x)
    for i=1:52
        x=sqrt(x);
    end
    for i=1:52
        x=x.^2;
    end
    res=x;
end