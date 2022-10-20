syms x1 x2 u a b
vars = [x1, x2];
odes = [-x1+a*x2*(1+x1); -b*x1*(1+x1)+u];
% Revisamos que el Jacobiano esté bien:
J = jacobian(odes, vars);
syms T
A = eye(2) + J*T;
%A(T,x1,x2): función de la matriz del modelo linealizado y discretizado
A = matlabFunction(A); 
Ts = 0.1; % tiempo de muestreo T.
tiempo = 0:Ts:50; % periodo de tiempo a usar.
x0 = [1, 1]; % Condición inicial
% Probamos dos entradas:
u = sin(5*tiempo) + randn(1,length(tiempo));%ones(1, length(tiempo));%sin(5*tiempo) + randn(1,length(tiempo));
plot(tiempo, u)
title('Entrada');
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$u(k)$', 'Interpreter', 'Latex', 'FontSize', 15);
grid on
odes= matlabFunction(odes)
a = 0.01; b = 0.02;
[t, x] = ode45(@(t,x) parcial(t, x, u, tiempo, odes,a,b), tiempo, x0);
y = x(:,1);
% Condiciones iniciales dados por el usuario:
M = 0.01*eye(2); N = 0.5; P0 = M; x0 = [1 1]'; Ts = 0.1;
% Parámetros conocidos:
C = [1 0];
n = length(y); x_est = zeros(n,2); P_trace = zeros(n,1); K_norm = zeros(n,1);
x_act = x0; P_act = P0;
for k = 1:n
    A_ = A(Ts,a,b,x_act(1), x_act(2));
    % Predicción:
    x_pred = [x_act(1) + Ts*(a*x_act(2)*(1+x_act(1))-x_act(1));
              x_act(2) + Ts*(u(k)-b*x_act(1)*(1+x_act(1)))];
    P_pred = A_*P_act*A_' + M;
    % Corrección:
    K = P_pred*C'*inv(C*P_pred*C' + N);
    x_act = x_pred + K*(y(k) - C*x_pred);
    P_act = (eye(2) - K*C)*P_pred;
    % Estimación:
    x_est(k,:) = x_act';
    % Medidas:
    P_trace(k) = trace(P_act);
    K_norm(k) = norm(K);
end
plot(t, x(:,2), t, x_est(:,2))
legend({'Modelo original', 'Estimación Kalman Extendido'}, 'Location', 'best');
xlabel('$k$', 'Interpreter', 'Latex', 'FontSize', 15);
ylabel('$x_{1}(k)$', 'Interpreter', 'Latex', 'FontSize', 15);
title('Comparación modelo original con estimación kalman extendido', 'FontSize', 12);
grid on