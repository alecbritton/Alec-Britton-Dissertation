function result = differentiate(inp,order)
    syms f(x);
    lambda = 4.370041;
    alpha_n = (sinh(lambda) - sin(lambda))/(cos(lambda) - cosh(lambda));
    f(x) = (sin(lambda*x) - sinh(lambda*x)) + alpha_n * (cos(lambda*x) - cosh(lambda*x));
    diffMode = diff(f,x,order);
    result = double(diffMode(inp));
end