function [frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, dt, snapshot, snapshotName, S1,S2,S3,movie,C_p,C_s, decibel_txt] = input1()
%Input parameters for point source
frequency = 2;
x0 = 5;
y0 = 5;

%Defining number of gridpoints
m_x = 60;
m_y = 15;

%Defining lenght of domain
W_x = 40;
W_y = 10;

roomLength = 0.5*W_y;
roomHeight = 0.5*W_y;
roomX = 0.8*W_x;
roomY = 0.25*W_y;

C_p = 1500; %P_wave
C_s = 800; % S-wave
rho = 2000; %Density

%Material 'constants', will vary by function
mu = rho*C_s^2;      %Second lamé parameter
lambda = rho*C_p^2-2*mu;  %First lamé parameter

%Which BC,[NBC, ABC, NABC]
BC = 'ABS';

%If eigenvalue plot is wanted [Y,N]
eig_answer = 'N';

%Whats going to be plotted
% 1: "V_x"
% 2: "V_y"
% 3: "V"
% 4: "Sigma_xx"
% 5: "Sigma_xy"
% 6:  "Sigma_xx"
% 7: "P"
plot_answer = 3;



%timestep
dt = 0.00001;

%Simulation lenght
T = 0.1;
T = T - mod(T, dt);


%%%SAVING DATA STUFF%%%
%Snapshot wanted?
snapshot = "Y";
snapshotName = "input0";

S1 = (x0+2)/C_p - mod((x0+2)/C_p, dt);
S2 = sqrt((roomY - y0)^2+(roomX-x0)^2)/C_p - mod(sqrt((roomY - y0)^2+(roomX-x0)^2)/C_p, dt);
S3 = T - mod(T, dt);

%movie?
movie = "Y";

%Saving decibel to txt file?
decibel_txt = "N";

end