% m: number of gridpoints
% W_x: length of domain
% W_y: height of domain
% L: Lenght of room
% H: Height of room

m_x = 10;
m_y = 15;

W_x = 20;
W_y = 20;

%gridspacing for height and length
hx = W_x/(m_x-1);
hy = W_y/(m_y-1);

x0 = 10;
y0 = 10;

sig = 1;

Id = speye(m_x*m_y);

Ze = sparse(m_x*m_y,m_x*m_y);

Id0 = 1+sparse(m_x*m_y,1);

Ze0 = sparse(m_x*m_y,1);

e1 = [Id,Ze,Ze,Ze,Ze];
e2 = [Ze,Id,Ze,Ze,Ze];
e3 = [Ze,Ze,Id,Ze,Ze];
e4 = [Ze,Ze,Ze,Id,Ze];
e5 = [Ze,Ze,Ze,Ze,Id];

e01 = [Id0,Ze0,Ze0,Ze0,Ze0];
e02 = [Ze0,Id0,Ze0,Ze0,Ze0];
e03 = [Ze0,Ze0,Id0,Ze0,Ze0];
e04 = [Ze0,Ze0,Ze0,Id0,Ze0];
e05 = [Ze0,Ze0,Ze0,Ze0,Id0];

x = 0 : hx : W_x;
y = 0 : hy : W_y;



[X, Y] = ndgrid(x, y);

Z = exp(-0.7*((X-x0).^2+(Y-y0).^2));
v_x0 = sparse(Z(:));
v_y0 = sparse(Z(:));

sigma_xx0 = sparse(m_x*m_y,1);
sigma_xy0 = sparse(m_x*m_y,1);
sigma_yy0 = sparse(m_x*m_y,1);


u = e1'*v_x0+e2'*v_y0+e3'*sigma_xx0+e4'*sigma_xy0+e5'*sigma_yy0;

TEST = e01'*v_x0;


[H_x, HI_x, Dp_x, Dm_x, e1_x, em_x] = d1_upwind_2(m_x, hx);
[H_y, HI_y, Dp_y, Dm_y, e1_y, em_y] = d1_upwind_2(m_y, hy);

%Flippa ordning
Dp_y2 = kron(Dp_y,speye(m_x) );
Dm_y2 = kron(Dm_y, speye(m_x));
Dp_x2 = kron(speye(m_y), Dp_x);
Dm_x2 = kron(speye(m_y),Dm_x );

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

A = [Ze,           Ze, Id*Dp_x2, Ze, Ze;
     Ze,           Ze, Ze, Id*Dp_x2, Ze;
     (LA+2*MU)*Dm_x2, Ze, Ze, Ze, Ze;
     Ze,           MU*Dm_x2, Ze, Ze, Ze;
     LA*Dm_x2,      Ze, Ze, Ze, Ze];

B = [Ze, Ze,          Ze, Id*Dp_y2, Ze;
     Ze, Ze,          Ze, Ze, Id*Dp_y2;
     Ze, LA*Dm_y2,     Ze, Ze, Ze;
     MU*Dm_y2, Ze,         Ze, Ze, Ze;
     Ze, (LA+2*MU)*Dm_y2, Ze, Ze, Ze];


%Characteristic Boundary operator
[L_w, L_e, L_s, L_n] = BC_fun(e1, e2, e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y);

L_x = [L_w; L_e];
L_y = [L_s; L_n];

H = kron(H_y,H_x);
H5 = kron( speye(5), H);
%Discrete innerproduct
HII = inv(H5*C);

I = speye(5*m_x*m_y);

%Defining projectionoperator, one for each direction is needed
P_x = I - HII * L_x' * inv(L_x * HII * L_x') * L_x; %#ok<MINV> 
P_y = I - HII * L_y' * inv(L_y * HII * L_y') * L_y; %#ok<MINV> 



P = P_x*P_y;

D_x = A; 
D_y = B; 

disp(em_y)

M = C\(D_x+D_y);


OP = P*M*P;

%u_t = B*U;

%RK4

%Kör med små m
%GÖR ALLT SPARSE FORMAT
%No positive real-parts for stability
%CFL relevant
%testa köra längre simulering

figure;
ee=eig(full(H5*C*OP));
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
t = 0; T = 20;
dt = 0.001; %ändra efter CFL
count = 0;

plotDatax = reshape(e1*u, m_x, m_y);

plotDatay = reshape(e2*u, m_x, m_y);

plotData = sqrt(plotDatax.^2 + plotDatay.^2);

%IC plot
figure;
hSurf = surf(X,Y,plotData);
colorbar;
caxis([0 0.14])
axis tight;
xlim([0 W_x])
ylim([0 W_y])
zlim([-0.3 1.5])
xlabel('x'); ylabel('y'); zlabel('x-hastighet')

while t < T
    u_next = RK4(u,dt,OP);
    t = dt + t;
    count = count +1;

    plotDatax = reshape(e1*u, m_x, m_y);
    plotDatay = reshape(e2*u, m_x, m_y);

    plotData = sqrt(plotDatax.^2 + plotDatay.^2);

    if mod(count,10) == 0
        set(hSurf, 'ZData', plotData);
        drawnow;
    end
    u = u_next;
end

%Vad kan vara fel?
% - disp med input vad man vill plotta
% - Snygga till koden
% - Få in ett nytt rum i domänen
% - Testa olika komponenter till 0, är nån stabil? Mental sjukhus?


