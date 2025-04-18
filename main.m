% m: number of gridpoints
% W_x: length of domain
% W_y: height of domain
% L: Lenght of room
% H: Height of room

u = [v_x, v_y, sigma_xx, sigma_xy, sigma_yy];

m_x = 51;
m_y = 51;

W_x = 10;
W_y = 20;

%gridspacing for height and length
hx = W_x/(m_x-1);
hy = W_y/(m_y-1);

d1_upwind_2;

%Material 'constants', will vary by function
lambda = 1;  %First lamé parameter
mu = 1;      %Second lamé parameter
rho = 1;     %Density

LA = lambdaDiagonal(lambda,m_x,m_y);
MU = lambdaDiagonal(mu,m_x,m_y);
RH = lambdaDiagonal(rho,m_x,m_y);

Id = speye(m_x*m_y);

Ze = sparse(m_x,m_y);

%Defining matrixes for equation 11 in report
C = [RH, Ze,   Ze, Ze, Ze;
     Ze,   RH, Ze, Ze, Ze;
     Ze,   Ze,   Id, Ze, Ze;
     Ze,   Ze,   Ze, Id, Ze;
     Ze,   Ze,   Ze, Ze, Id];

A = [Ze,           Ze, Id, Ze, Ze;
     Ze,           Ze, Ze, Id, Ze;
     LA+2*MU, Ze, Ze, Ze, Ze;
     Ze,           MU, Ze, Ze, Ze;
     LA,      Ze, Ze, Ze, Ze];

B = [Ze, Ze,          Ze, Id, Ze;
     Ze, Ze,          Ze, Ze, Id;
     Ze, LA,     Ze, Ze, Ze;
     MU, Ze,         Ze, Ze, Ze;
     Ze, LA+2*MU, Ze, Ze, Ze];

D = [Dp, Ze, Ze, Ze, Ze;
     Ze, Dp, Ze, Ze, Ze;
     Ze, Ze, Dm, Ze, Ze;
     Ze, Ze, Ze, Dm, Ze;
     Ze, Ze, Ze, Ze, Dm];


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

D_x = A*D; %NEEDS TO BE DEFINED
D_y = B*D; %NEEDS TO BE DEFINED


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





