% To get the data of the messages.mat
%---> Find the pdf from the message.
%---> Find the likelihood ratio & MLD rule.
%---> Find the error probability.
clear
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

% To get the data of the signal.mat
%---> Decide the signal from the signal+noise data.
%---> How many errors will be the result?
load signal.mat;

figure(4)
plot(t_s, s_n, '-')
xlabel('Time [sec]')
ylabel('Voltage [v]')