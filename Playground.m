clear all
clc
close all

%% Definitions and parameters
hbar = 6.62606876e-34 / (2*pi); %Planck constant
S_x = (1/2) * [0 1; 1 0];
S_y = (1/2) * [0 -1i; 1i 0];  % 1i is the imaginary unit
S_z = (1/2) * [1 0; 0 -1];      
plus = [1; 0];
minus = [0; 1];
B_0 = 1;     % Static magnetic field = 1T
B_1 = 2e-4;  % Oscillating magnetic field, at frequency omega_0
gamma = 2*pi*28e9;  % electron gyromagnetic ratio = 28 GHz/T
omega_0 = gamma * (B_0);   % Larmor frequency of a spin in the field B_0
omega_1 = gamma * B_1;   % Rabi frequency. This is the frequency at which a 
                         % spin would precess if placed in a static field B_1
%% Time evolution
dt = 2e-12;  % time sampling interval. How do you choose this value properly?


H0 = (omega_0) * S_z;

half_period = 8.93e-8;
psi0 = minus;   %initial state.
time = 0:dt:half_period;


noisePwr = 100e-6;

avgSz=zeros(length(time),1);
avgSx=zeros(length(time),1);
avgSy=zeros(length(time),1);
B_n = noisePwr*wgn(1,length(time),0);
for t = 1:length(time)
   % Hamiltonian describing the rotating magnetic field
   %H1 = omega_1 * (S_x * cos(omega_0 * time(t)) + S_y * sin(omega_0 * time(t))); 
   H1 = omega_1 * (S_x * cos(omega_0 * time(t))); 
   U = expm(-1i * (H0 + H1 + B_n(t)*gamma*S_z) * dt); 
   psi0 = U * psi0; % time-evolved state
end


tau = 4e-11;
dt2 = 3e-13;
rotateTime = 0:dt2:tau;
%Free Evolution
B_n = noisePwr*wgn(1,length(rotateTime),0);
avgSz = zeros(length(rotateTime),1);
avgSx=zeros(length(rotateTime),1);
avgSy=zeros(length(rotateTime),1);




for t = 1:length(rotateTime)
   % Hamiltonian describing the rotating magnetic field
   %H1 = omega_1 * (S_x * cos(omega_0 * time(t)) + S_y * sin(omega_0 * time(t))); 
   U = expm(-1i * (H0 + B_n(t)*gamma*S_z) * dt2); 
   psi0 = U * psi0; % time-evolved state
   avgSx(t) = psi0' * S_x * psi0;
   avgSy(t) = psi0' * S_y * psi0;
   avgSz(t) = psi0' * S_z * psi0;
end

figure(1)
subplot(1,3,1)
plot(rotateTime,avgSx)
xlabel('seconds(s)')
ylabel('<Sx>')
subplot(1,3,2)
plot(rotateTime,avgSy)
xlabel('seconds(s)')
ylabel('<Sy>')
subplot(1,3,3)
xlabel('seconds(s)')
ylabel('<Sz>')
plot(rotateTime,avgSz)
