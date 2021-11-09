function xdot = nonlinear_ode(t,x)

if t<=2
    u=25;
else
    u=35;
end

xdot = zeros(3,1);
xdot(1) = x(2);
xdot(2) = x(3);
xdot(3) = (3/5)*u - 8*x(3)- 30*x(2) - 3*x(1)^4;