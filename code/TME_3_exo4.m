syms f(x)
f(x) = x^2 - 2;

% Q2.
% converge dans Z_7
newton_padique(f, 7, 3, 30)

% Q1.
function x = newton_padique(f, p, a, itmax)
    g = diff(f);
    i = mulinv(g(a), p);
    x = sym(a);
    r = sym(p);
    for j=1:itmax
        x = sym(x + r*mod(-f(x)*i/r, p));
        fprintf("x_%d = %s, f(x_%d) = %d mod %s\n", j, char(x), j, mod(f(x), p*r), char(p*r));
        r = p*r;
    end
    fprintf("f(⋅) = %d\n", mod(f(x), r));
end

function y = mulinv(x,p)
    if ~isprime(p)
        disp('The field order is not a prime number');
        return
    elseif x>=p
        disp('All or some of the numbers do not belong to the field');
        return
    elseif x == 0
        disp('0 does not have a multiplicative inverse');
        return
    end
    k = 0;
    m=mod(k*p+1,x);
    while m       
        k=k+sign(m);
        m=mod(k*p+1,x);
    end
    y=(k*p+1)./x;
end