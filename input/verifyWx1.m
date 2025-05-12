function [frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, snapshot, snapshotName,movie,C_p,C_s,decibel_txt, save_movie] = verifyWx1()
%Input parameters for point source
frequency = 50; %1,5,10,15
x0 = 5; %Låst
y0 = 5; %Låst



roomLength = 5; %Låst
roomHeight = 3; %Låst
roomX = x0 + 100; % (2,5)


%Börja högt och jobba ner
%Börja med låg W_y för när vi varierar W_x
%Defining lenght of domain
W_x = roomX + 30; %Låst, VERIFIERA rum_pos + (30, 60, 90, 120)
W_y = 40; %Låst, VERIFIERA (10,20,30,40)


roomY = ( W_y - roomHeight) - 2; %testa 2 olika (2,12)


%Defining number of gridpoints, TESTA SIST MED OLIKA
m_x = (10*W_x*frequency)/332; %ska bero på våglängde
m_y = (10*W_y*frequency)/332; %ska bero på våglängden (minst 40)
if m_y < 40
    m_y = 40;
end


m_x = ceil(m_x);
m_y = ceil(m_y);

C_p = 1500; %P_wave,  OM VI HINNER (granit)
C_s = 800; %S-wave, OM VI HINNER (granit om vi hinner)
rho = 2000; %Density, OM VI HINNER

%Material 'constants', will vary by function
mu = rho*C_s^2;      %Second lamé parameter
lambda = rho*C_p^2-2*mu;  %First lamé parameter

%Which BC,[NBC, ABC, NABC]
BC = 'NABC';

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
plot_answer = 7;

%Simulation lenght, ändra dt, 
T = 0.3;
%T = T - mod(T, dt);

%Snapshot wanted?
snapshot = "N";
snapshotName = "verifyWx1";

%movie?
movie = "Y";
save_movie = "N";

%save decibel data to txt?
decibel_txt = "Y";

end