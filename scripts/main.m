clc;
clear;

%Trabajo basado en:
%https://www.sciencedirect.com/science/article/pii/S0029801819304767

%Parametros del modelo
nx = 12;
nu = 6;
ny = 12;

% Puntos de equilibrio
thetas = [0; 0; 0];
omegas = [0; 0; 0];

%% Ecuacion de estado
Ts = 1; %Tiempo de muestreo

sys = ecuacion_estado(equils, Ts);

A_matriz_dt = sys.A;
B_matriz_dt = sys.B;
C_matriz_dt = sys.C;

%% Simulacion en Tiempo
T = 100; %Tiempo de simulacion. En Segundos.

p0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
dtao = [0.1; 0; 0; 0; 0; 0];
tao = equils(13:end) + dtao;

[t_nolin, x_nolin] = ode45(@(t,y) auv_system(t,y,tao), [0, T], ...
    p0);

x_dt = [];
x_dt(1, :) = p0;

for k=1:T/Ts
   x_dt(k+1, :) = A_matriz_dt*x_dt(k, :)' + B_matriz_dt*dtao;
end

x=zeros(size(x_dt));

figure;
hold on
plot(t_nolin, x_nolin(:, 1)', 'r', 'LineWidth', 1);
for i=1:nx
    x(:,i) = x_dt(:,i) + equils(i);
end
plot(0:Ts:T, x(:,1),'k--', 'LineWidth', 2);
legend('No Lin', 'Lin');
xlabel('Tiempo (s)');
ylabel('X (m)');


