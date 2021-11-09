%ENERGY 295 Spring 2021
%HW 1 Solution of Problem 1 Part 6 and Part 7


clear all; close all; clc

%% Initialization

t = [0:0.1:25]; % define time vector
 
% define linearized state space matrices
A_lin = [0,1,0;0,0,1;-12*(5)^0.75,-30,-8];
B_lin = [0,0,3/5]';
C_lin = eye(3); 
D_lin = [0,0,0]';
sys_model_lin = ss(A_lin,B_lin,C_lin,D_lin);

% At equilibrium
u_e  = 25; % input at equilibrium
y_e = [5^(0.25),0, 0]'; % output at equilibrium 
x_e = [5^(0.25),0, 0]; % states at equilibrium

%% Transfer Function

G_s=tf(sys_model_lin)


%% Simulation of Linearized System

% Assuming the step input starts at t=2 seconds,i.e. the input \delta u is zero at
% time zero (i.e. u=25) until time =2sec and then it goes to 10 (i.e. u=35)

delta_u = 10*[zeros(1,20),ones(1,length(t)-20)]; 
delta_y = lsim(sys_model_lin,delta_u,t);

% For model output, add the equilirium y to delta_y
% Parse outputs (y1, y2, y3)
y_linear_1 = y_e(1) + delta_y(:,1);
y_linear_2 = y_e(2) + delta_y(:,2);
y_linear_3 = y_e(3) + delta_y(:,3);
% 

%% Simulation of Nonlinear System

[T,x_nonlin] = ode45(@nonlinear_ode, t, x_e);

% Parse outputs (y1, y2, y3)
y_nonlin_1 = x_nonlin(:,1);
y_nonlin_2 = x_nonlin(:,2);
y_nonlin_3 = x_nonlin(:,3);

%% Comparison plots

figure(1)
plot(t,y_linear_1,T,y_nonlin_1,'--','LineWidth',2);
xlabel('Time [s]','FontSize',16);ylabel('Output - y_{1}','FontSize',16);
legend('linearized','nonlinear'); set(gca,'FontSize',16); grid on

figure(2)
plot(t,y_linear_2,T,y_nonlin_2,'--','LineWidth',2);
xlabel('Time [s]','FontSize',16);ylabel('Output - y_{2}','FontSize',16);
legend('linearized','nonlinear'); set(gca,'FontSize',16); grid on

figure(3)
plot(t,y_linear_3,T,y_nonlin_3,'--','LineWidth',2);
xlabel('Time [s]','FontSize',16);ylabel('Output - y_{3}','FontSize',16);
legend('linearized','nonlinear'); set(gca,'FontSize',16); grid on
