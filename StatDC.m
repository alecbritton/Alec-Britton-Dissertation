function StatDC(L,b,E,h,g0)
    ws = 0;
    vDCi = 0;
    I = (1/12)*b*h^3;
    e0 = 8.854*(10^-12);
    alpha = ((L^4)*e0*b)/(2*E*I*(g0^3));
    lambda = 4.730041;
    alpha_n = (sin(lambda) - sinh(lambda))/(cosh(lambda) - cos(lambda));
    
    dx = 0.01;
    x = 0:dx:1;
    
    modeFun = (sin(lambda.*x) - sinh(lambda.*x)) + alpha_n * (cos(lambda.*x) - cosh(lambda.*x));
    modeFunSq = ((sin(lambda.*x) - sinh(lambda.*x)) + alpha_n * (cos(lambda.*x) - cosh(lambda.*x))).^2;
    km = sum(modeFun.*differentiate(x,4).*dx);
        
    for vDCi1 = 0.1:0.1:25
        disp(vDCi1);
        ke = alpha* (vDCi^2) * sum(((1-ws)^-3) .* modeFunSq.*dx);
        F = alpha*(vDCi1^2-vDCi^2)*sum(((1-ws)^-2).*modeFun.*dx);
        a = F/(km-ke);
        ws = (ws + a * modeFun(51));
        plot(vDCi1,ws,'b--.');
        title('Static analysis of electrostatically actuated microbeam');
        ylabel('Beam deflection');
        xlabel('Voltage / V');
        hold on
        vDCi = vDCi1;
    end
end

%L = 350*(10^-6);
%b = 50*10^-6;
%E = 169.9*10^9;
%h = 3*10^-6;
%g0 = 1*10^-6;
%StatDC(L,b,E,h,g0);