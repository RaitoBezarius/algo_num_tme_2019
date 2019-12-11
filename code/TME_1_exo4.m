x=0:1:2000;
y=test_sinus(x);
z=test_sinus_taylor(x);
err=zeros(size(z));
for i=1:size(err)
    err(i) = abs(y(i) - z(i))/y(i);
end
errorbar(x, y, err, '-');
title('Sinus avec erreurs relatives');
legend('sin(-10 + k/100)');



function res = sinus_taylor1(x)
    acc = 0;
    for i=1:15
        acc = acc + (-1)^i * (x.^(2*i + 1)) / factorial(2*i + 1);
    end
    res = acc;
end

function res = sinus_taylor2(x)
    r = mod(x, 2*pi);
    acc = 0;
    i = 0;
    limit_reached = false;
    while ~limit_reached
        new_acc = acc + (-1)^i * (r.^(2*i + 1)) / factorial(2*i + 1);
        if new_acc == acc
            limit_reached = true;
        end
        acc = new_acc;
        i = i + 1;
    end
    res = acc;
end
        
function res = test_sinus_taylor(k)
    res = sinus_taylor2(-10 + k/100);
end

function res = test_sinus(k)
    res = sin(-10 + k/100);
end

        