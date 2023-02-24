%% Auth: Arian Velayati, PhD
% The script can be used to find the ultrasonic elastic constants
...of the rock from s and p-wave traveling times

clc; close; clear;

%% Input

L = [15.24]; % Traveling Distance, mm
P = []; % Bulk Density, g/cc
Tp = [12]; % P-wave traveling time, us
Ts = [7.23]; % S-wave traveling time, us

% Face to Face Arrival
Tp0 = 24.374; % us
Ts10 = 33.29; % us
Ts20 = 38.42; % us

%% Calculations

P = P*1000; % Density in kg/m3
Vp = (L/1000)/((Tp-Tp0)/1e6); % Compression wave velocity, m/s
Vs = (L/1000)/((Ts-Ts0)/1e6); % Shear wave velocity, m/s

E = (P*Vs^2*(3*Vp^2-4*Vs^2))/(Vp^2-Vs^2); % Young's Modulus
G = P*Vs^2; % Shear Modulus
v = (Vp^2-2*Vs^2)/(2*(Vp^2-Vs^2)); % Poisson ratio
k = P*(3*Vp^2-4*Vs^2)/3; % Bulk Modulus
E = E/1e9; G = G/1e9; k = k/1e9; % Units in GPa

%% Output

   T = table(E,G,v,k,'VariableNames',{'YM_GPa','SM_GPa','PR','BM_GPa'})
  
    writetable(T,'ULT_Moduli.csv')
%% Signals

clc;clear;
% Waves = load('Allwaveforms.csv');
% Waves2 = Waves(32:end,:);
% Waves2 = table2array(Waves2);

filename = 'C:\Users\velayati\Documents\MATLAB\Allwaveforms.csv';
delimiter = ',';
startRow = 33;
% Format string for each line of text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Create output variable
Allwaveforms = [dataArray{1:end-1}];
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

Time = Allwaveforms(:,1);
[x,y] = size(Allwaveforms);

% Signal = Allwaveforms(:,5);
% plot(Time,Signal)

for i=2:y
Signal = Allwaveforms(:,i);
figure(i)
plot(Time,Signal)
hold on
end



