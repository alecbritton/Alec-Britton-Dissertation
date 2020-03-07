function DynamDC(e0,b,h,c,rho,E,w,Vdc)
I = (1/12)*b*h^3;
dx = 0.01;
x = 0:dx:1;
A = b*h;
lambda = 4.730041;
alpha_n = (sin(lambda) - sinh(lambda))/(cosh(lambda) - cos(lambda));

modeFun = (sin(lambda.*x) - sinh(lambda.*x)) + alpha_n * (cos(lambda.*x) - cosh(lambda.*x));
modeFunSq = ((sin(lambda.*x) - sinh(lambda.*x)) + alpha_n * (cos(lambda.*x) - cosh(lambda.*x))).^2;

a1 = E*I*sum(differentiate(modeFun,4).*modeFun.*dx);
a2 = rho*A*sum(modeFunSq*dx);
a3 = c*sum(modeFunSq*dx);
a4 = e0*b*Vdc^2*0.5*sum((1/(1-w)^2).*modeFun.*dx);

t = 0:0.5:3;
u0 = 0;
uDot0 = 0;

[t,u] = ode45(@ut, t, [u0,uDot0]);

plot(t,u(:,1));
xlabel('t');
ylabel('u');

    function dudt = ut(t,u)
        dudt1 = u(2);
        dudt2 = (1/a2)*(a4-a1*u(1)-a3*u(2));
        
        dudt = [dudt1,dudt2];
        dudt = dudt(:);
    end
end