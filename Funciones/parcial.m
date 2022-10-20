function dxdt = parcial(t, x, u, tiempo,odes,a,b)
    u_int = interp1(tiempo, u, t);
    dxdt = odes(a,b,u_int,x(1),x(2));
end