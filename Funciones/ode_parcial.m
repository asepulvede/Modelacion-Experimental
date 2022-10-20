function dxdt = ode_parcial(t, x, u, tiempo, a, b, c, d)
    u_int = interp1(tiempo, u, t);
    dxdt = [-a*x(1) + b*x(1)*x(2); c*x(1)^2 - d*x(2) + u_int];
end

