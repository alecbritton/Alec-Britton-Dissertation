function StatDC(L,b,E,h,g0)
    ws = 0;
    vDCi = 0;
    I = (1/12)*b*h^3;
    e0 = 8.854*(10^-12);
    alpha = ((L^4)*e0*b)/(2*E*I*(g0^3));
    lambda = 4.370041;
    alpha_n = (sinh(lambda) - sin(lambda))/(cos(lambda) - cosh(lambda));
    syms x;
    modeFun = @(x) (sin(lambda.*x) - sinh(lambda.*x)) + alpha_n * (cos(lambda.*x) - cosh(lambda.*x));
    modeFunSq = @(x) ((sin(lambda.*x) - sinh(lambda.*x)) + alpha_n * (cos(lambda.*x) - cosh(lambda.*x))).^2;
    km = integral(matlabFunction(diff(modeFun(x),4)),0,1);
    
    for vDCi1 = 0.1:0.1:20
        ke = alpha* vDCi^2 * (1-ws)^-3 * integral((modeFunSq),0,1);
        F = alpha*(vDCi1^2-vDCi^2)*(1-ws)^-2*integral(modeFun,0,1);
        if(km-ke ~= 0)
            a = F/(km-ke);
            ws = ws + a * modeFun(0.5);
            plot(vDCi1,ws,'.');
            hold on
        end
        vDCi = vDCi1;
    end
end