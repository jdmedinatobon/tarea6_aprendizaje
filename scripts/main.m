clc;
clear;

%Trabajo basado en:
%https://www.sciencedirect.com/science/article/pii/S0029801819304767
%El paper del ifac tambien

%Creo que lo que puedo hacer es utilizar lo de dinamicas poblacionales
%distribuidas discretas para generar las trayectorias de cada submarino y
%despues usar el MPC para lograr que siguan esas trayectorias.

%Iteraciones maximas.
max_iter = 30;

%Voy a hacer esa parte primero.
%Por ahora con 4 robots.
n = 4;

%Numero de dimensiones.
dim = 6;

%Dimensiones de una piscina olimpica
global beta%beta_x beta_y beta_z

beta = [50; 25; 2; 2*pi; 2*pi; 2*pi];

% beta_x = 50;
% beta_y = 25;
% beta_z = 2;

%Con Smith dynamics saturadas
global deltas pos_actual pos_lider_actual

%Variable global con las desviaciones con respecto al lider de cada
%seguidor. El delta en la primera posicion siempre es 0, ya que corresponde
%al lider.
deltas = [0 ,   5 , 10 ,   15  ;
          0 ,  2.5,  5 ,   7.5 ;
          0 , 0.25, 0.5,  0.75 ;
          0 , pi/2,  pi, 3*pi/2;
          0 ,   0 ,  0 ,   0   ;
          0 ,   0 ,  0 ,   0  ];

%Variable global con las posiciones de cada robot y de la variable 
%adicional.
%Esto corresponde al r en el paper. La primera posicion es la variable
%adicional.
pos_inicial(:,:,1) = [ 140,  10 ,  20, 30 ;
                        70,   5 ,  10, 15 ;
                        7 ,  0.5,  1 , 1.5;
                      n*2*pi,   0 ,  0 ,  0 ;
                      n*2*pi,   0 ,  0 ,  0 ;
                      n*2*pi,   0 ,  0 ,  0 ;];
            
pos_actual = pos_inicial(:,:,1);

pos = cat(3, pos_inicial, zeros(dim, n, max_iter));

pos_lider_inicial = [0; 0; 0; 0; 0; 0];
pos_lider_actual = pos_lider_inicial;

pos_lider = [pos_lider_inicial, zeros(dim, max_iter)];

%Cell con los vecinos de cada robot.
neighbors = {[2 4], [1 3], [2 4], [1 4]};

epsilon = calcular_epsilon(neighbors);

aux = pos_actual;

%Hacer esto para todas las dimensiones y depronto tambien los angulos.
for k=2:max_iter+1
   
%     if k < 150
%         pos_lider_actual = (5/150)*k + 2*sin((2*pi/75)*k);
%     else
%         pos_lider_actual = 2*sin((2*pi/75)*(k-150)) + 5;
%     end
    
    pos_lider(:, k) = pos_lider_actual;
    
    for i=1:n
        suma = zeros(dim, 1);
        robot_neighbors = neighbors{i};
        fi = fitness(i);
        
        for r=1:length(robot_neighbors)
            j = robot_neighbors(r);
            
            fj = fitness(j);
            suma = suma + (fi-fj).*theta(i, j, fi, fj).*phi(fi, fj);
        end

    aux(:, i) = pos_actual(:, i) + epsilon(:, i).*suma;
    pos(:, i, k) = aux(:, i);
    
    end
    
    pos_actual = aux;
end

disp("Termino");

%Graficas
labels = string(2:1:n);
ejes_y = ["x", "y", "z", "\phi_x", "\phi_y", "\phi_z"];
titulos = ["Posicion en", "Posicion en", "Posicion en",...
           "Orientacion en", "Orientacion en", "Orientacion en "];
leyenda = string(repmat("Submarino", 1, n-1));
leyenda = join([leyenda; labels], " ", 1);
leyenda = ["Submarino Lider", leyenda];
% leyenda = ["Submarino Lider", leyenda , "Variable Adicional"];

figure(1);
set(1, "defaultAxesFontSize", 12);

for i=1:dim
    subplot(2,3,i);
    hold on
    
    plot(pos_lider(i,:));

    for j=2:n
        plot(squeeze(pos(i, j, :)));
    end

    %plot(pos_x(:,1));
    legend(leyenda);
    xlabel("Iteraciones");
    ylabel("Posicion en " + ejes_y(i));
    grid on
    title(titulos(i) + " " + ejes_y(i) + " de cada submarino");
end

%% MPC

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