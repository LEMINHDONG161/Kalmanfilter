clear all;
close all;
clc;
load signal.mat;
% *********** Measurement data ***********
T = 0.1; %( step size)
Time = [0:T:100];
N = length(t_s);

%R_noise = sqrt(5) * randn(1,N);
%z = 10 + 5*cos(2*pi*Time/25) + 5*cos(2*pi*Time/50) + R_noise; % To generate a N array of random numbers with a zero mean that are distributed with a variance of 1
z=s_n;
%z_true = 10 + 5*cos(2*pi*Time/25) + 5*cos(2*pi*Time/50);
% *********** Initial Parameter for Kalman Filtering ***********

F = [1 T; 0 1]; % transition matrix
H = [1 0];      % measurement matrix 

% measurement noise covariance 
R = 13; 

% process noise covariance 
Q = 0.7;

P = [0 0; 0 0];


x_initial = [1.5 ; 0]; % initial state vector
x_hat = x_initial; % initial state estimate

% ****************** Kalman Filtering ********************
for k = 1:N
    % Update the most recent state estimate to the present time.
    x_hat = (F * x_hat) ;

    % initial estimation covariance
    P_k = F*P*F' + Q;

    %  covariance of the Correction Vector
    Re = (H * P_k * H') + R ;

    %  Kalman Gain matrix.
    K_k = P_k * H'* inv(Re);

    % Update the state estimate.
    x_hat = x_hat + (K_k * (z(k) - (H * x_hat)));

    % Compute the covariance of the estimation error.
    I = eye(2,2);
    P_k = (I - (K_k * H)) * P_k;

    xest(:,k)=x_hat;
end

% *******************   Error ***********************

%error_Kalman = mean(sqrt((z_true-xest(1,:)).^2))     % Kalman 

% ********************** Figure *****************
figure;
plot (t_s,z,'-')
hold on 
plot ( t_s,xest(1,:),'r-', 'linewidth',3); % title('Kalman Filter')  % Kalman 
xlabel('Time [sec]', 'fontsize',16)
ylabel('Voltage [V]', 'fontsize',16)
legend('Measured data','Kalman')
set(gca, 'fontsize', 16);
hold off

%savefile = 'measure_data.mat';
%save(savefile, 'Time', 'z');