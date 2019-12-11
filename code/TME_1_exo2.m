format longE

fprintf("Precision machine:\n")
fprintf("\tCalculee manuellement: \t%.16g\n", prec_machine)
fprintf("\tCalculée théoriquement:\t%.16g\n", 2^(-52))
fprintf("\tConstante MATLAB: \t%.16g\n\n", eps)

fprintf("Plus petit nombre normalise\n")
fprintf("\tCalculée manuellement: \t%.16g\n", normalized_min)
fprintf("\tCalculée théoriquement:\t%.16g\n", 2^(-1022))
fprintf("\tConstante MATLAB: \t%.16g\n\n", realmin)

fprintf("Plus petit flottant denormalise\n")
fprintf("\tCalculée manuellement: \t%.16g\n", denormalized_min)
fprintf("\tCalculée théoriquement:\t%.16g\n", 2^(-1074))
fprintf("\tConstante MATLAB: \t%.16g\n\n", realmin*eps)

fprintf("Plus grand nombre flottant\n")
fprintf("\tCalculée manuellement: \t%.16g\n", float_max)
fprintf("\tCalculée théoriquement:\t%.16g\n", 2^(1023)*(2 - eps))
fprintf("\tConstante MATLAB: \t%.16g\n\n", realmax)

function res = prec_machine()
    e = 1;
    while e + 1 > 1
        res = e;
        e = e/2;
    end
end

function res = normalized_min()
    e = eps;
    var = 1;
    while var + e ~= var
        res = var;
        e = e / 2;
        var = var / 2;
    end
end

function res = denormalized_min()
    var = 1;
    while var > 0
        res = var;
        var = var / 2;
    end
end

function res = max_float()
    x = 2 - eps;
    while x < inf
        res = x;
        x = x*2;
    end
end
