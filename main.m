% m: number of gridpoints
% W_x: length of domain
% W_y: height of domain
% L: Lenght of room
% H: Height of room

u = [v_x, v_y, sigma_xx, sigma_xy, sigma_yy];

m_x = 51;
my = 51;

W_x = 10;
W_y = 20;

%gridspacing for height and length
hx = W_x/(m_x-1);
hy = W_y/(m_y-1);

%Material 'constants', will vary by function
lambda = 1;  %First lamé parameter
mu = 1;      %Second lamé parameter
rho = 1;     %Density

%Defining matrixes for equation 11 in report
C = [rho, 0,   0, 0, 0;
     0,   rho, 0, 0, 0;
     0,   0,   1, 0, 0;
     0,   0,   0, 1, 0;
     0,   0,   0, 0, 1];

A = [0,           0, 1, 0, 0;
     0,           0, 0, 1, 0;
     lambda+2*mu, 0, 0, 0, 0;
     0,           mu, 0, 0, 0;
     lambda,      0, 0, 0, 0];

B = [0, 0,          0, 1, 0;
     0, 0,          0, 0, 1;
     0, lambda,     0, 0, 0;
     mu, 0,         0, 0, 0;
     0, lambda+2*mu, 0, 0, 0];

%Characteristic Boundary operator 
%NEEDS TO BE CORRECTLY DEFINED
L_x = [L_w; L_e];
L_y = [L_s; L_n];

%Discrete innerproduct
HII = inv(kron(H,H)*C);

%Defining projectionoperator, one for each direction is needed
P_x = I - HII * L_x' * inv(L_x * HII * L_x') * L_x; %#ok<MINV> 
P_y = I - HII * L_y' * inv(L_y * HII * L_y') * L_y; %#ok<MINV> 

P = P_x*P_y;

D_x = []; %NEEDS TO BE DEFINED
D_y = []; %NEEDS TO BE DEFINED


M = C\(D_x+D_y);

B = P*M*P;

%Kör med små m
%GÖR ALLT SPARSE FORMAT
%No positive real-parts for stability
%CFL relevant
figure;
ee=eig(h*B);
plot(real(ee),imag(ee),'*','MarkerSize',8);

%CHANGE TITLE FOR UNDERSTANDING
title('Eigenvalues to h\cdot B')

xlabel('\Re(h\cdotB)');ylabel('\Im(h\cdotB)');
ax = gca;          % current axes
ax.FontSize = 10;





