% m: number of gridpoints
% W_x: length of domain
% W_y: height of domain
% L: Lenght of room
% H: Height of room

m_x = 51;
m_y = 51;

W_x = 10;
W_y = 20;

Id = speye(m_x*m_y);

Ze = sparse(m_x,m_y);

e1 = [Id,Ze,Ze,Ze,Ze];
e2 = [Ze,Id,Ze,Ze,Ze];
e3 = [Ze,Ze,Id,Ze,Ze];
e4 = [Ze,Ze,Ze,Id,Ze];
e5 = [Ze,Ze,Ze,Ze,Id];

%u = [v_x, v_y, sigma_xx, sigma_xy, sigma_yy];



%gridspacing for height and length
hx = W_x/(m_x-1);
hy = W_y/(m_y-1);

x = 0 : hx : W_x;
y = 0 : hy : W_y;
[X, Y] = ndgrid(x, y);

[H_x, HI_x, Dp_x, Dm_x, e1_x, em_x] = d1_upwind_2(m_x, hx);
[H_y, HI_y, Dp_y, Dm_y, e1_y, em_y] = d1_upwind_2(m_y, hy);

%Material 'constants', will vary by function
lambda = 1;  %First lamé parameter
mu = 1;      %Second lamé parameter
rho = 1;     %Density

LA = lambdaDiagonal(lambda,m_x,m_y);
MU = lambdaDiagonal(mu,m_x,m_y);
RH = lambdaDiagonal(rho,m_x,m_y);


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
BC_fun;

L_x = [L_w; L_e];
L_y = [L_s; L_n];

%Discrete innerproduct
HII = inv(kron(H,H)*C);

%Defining projectionoperator, one for each direction is needed
P_x = I - HII * L_x' * inv(L_x * HII * L_x') * L_x; %#ok<MINV> 
P_y = I - HII * L_y' * inv(L_y * HII * L_y') * L_y; %#ok<MINV> 

P = P_x*P_y;

D_x = A*D; 
D_y = B*D; 


M = C\(D_x+D_y);

B = P*M*P;

u_t = B*U;

%RK4

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



%Animation for magintude if velocity
% Kvar att göra:
% - Anpassa kod för vårt fall
% - IC
% - Randvillkor, fixat
% - Hitta dt enligt cfl
% - Skriva kod för simulation



%INSERT IC



%Simulation stuff
t = 0; T = 2;

%ny vx och vy för varje steg.
%Beräkna magnituden för varje punkt



