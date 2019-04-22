%% Home work1 
% To get the data of the messages.mat
%---> Find the pdf from the message.
%---> Find the likelihood ratio & MLD rule.
%---> Find the error probability.
clear all
close all
load messages.mat; 

figure(1)
plot(t_m, m0, '--', t_m, m1, '-')

xlabel('Time [sec]')
ylabel('Voltage [v]')

figure(2) 
hist(m0);
xlabel('Voltage [v]')
ylabel('Number')

figure(3) 
hist(m1);
xlabel('Voltage [v]')
ylabel('Number')
%% Find the pdf from message m0 & m1
figure(4)
set(gcf,'numbertitle','off','name', 'Probability Distribution messasge m0 & m1');
[a,c]=hist(m0); 
b=a./20001;
plot(c,b,'r-');
xlabel('Voltage [v]')
ylabel('Probability')
hold on
[mean_m0,var_m0,muci,sigmaci] = normfit(m0)
[bincounts,binpositions] = hist(m0);
binwidth = binpositions(2) - binpositions(1);
histarea = binwidth*sum(bincounts)/20001;
x = binpositions(1):0.001:binpositions(end);
y = normpdf(x,mean_m0,var_m0);
plot(x,histarea*y,'r','LineWidth',2)
hold on
[a,c]=hist(m1); 
b=a./20001;
plot(c,b,'b-');
xlabel('Voltage [v]')
ylabel('Probability')
hold on
[mean_m1,var_m1,muci1,sigmaci1] = normfit(m1)
[bincounts,binpositions] = hist(m1);
binwidth = binpositions(2) - binpositions(1);
histarea = binwidth*sum(bincounts)/20001;
x = binpositions(1):0.001:binpositions(end);
y = normpdf(x,mean_m1,var_m1);
plot(x,histarea*y,'b','LineWidth',2)
%% Find the error m0 & m1
x = -5:0.1:5;

error_d1_m0 = quad('4/3*exp((-x.^2)/0.18)',0.65,5);
error_d0_m1 =quad('4/3*exp((-(x-1).^2)/0.18)',-5,0.5);




%% Signal & Noise
% To get the data of the signal.mat
%---> Decide the signal from the signal+noise data.
%---> How many errors will be the result?
load signal.mat;

figure,
plot(t_s, s_n, '-')
xlabel('Time [sec]')
ylabel('Voltage [v]')
figure, 
hist(s_n);
xlabel('Voltage [v]')
ylabel('Number')

%% KALMAN FILTER
T = 0.1; %( step size)
Time = [0:T:200];
N = length(t_s);
z=s_n;
%z_true = 10 + 5*cos(2*pi*Time/25) + 5*cos(2*pi*Time/50);
% *********** Initial Parameter for Kalman Filtering ***********

F = [1 T; 0 1]; % transition matrix
H = [1 0];      % measurement matrix 

% measurement noise covariance 
R = 15; 

% process noise covariance 
Q = 0.15;

P = [0 0; 0 0];


x_initial = [0 ; 0]; % initial state vector
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

error_Kalman = mean(sqrt((z-xest(1,:)).^2))     % Kalman 

% ********************** Figure *****************
figure,
plot (t_s,z,'-')
hold on 
plot ( t_s,xest(1,:),'r-', 'linewidth',3); % title('Kalman Filter')  % Kalman 
xlabel('Time [sec]', 'fontsize',16)
ylabel('Voltage [V]', 'fontsize',16)
legend('Measured data','Kalman')
set(gca, 'fontsize', 16);
hold off
