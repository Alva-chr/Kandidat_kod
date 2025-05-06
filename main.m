

%Defining number of gridpoints
m_x = 60;
m_y = 40;

%Defining lenght of domain
W_x = 30;
W_y = 20;

roomLength = 0.5*W_y;
roomHeight = 0.5*W_y;
roomX = 0.8*W_x;
roomY = 0.25*W_y;

%gridspacing for height and length
hx = W_x/(m_x-1);
hy = W_y/(m_y-1);

%Defining startpoint of gausspoint
x0 = W_y/2;
y0 = W_y/2;

sig = 1;

Id = speye(m_x*m_y);

Ze = sparse(m_x*m_y,m_x*m_y);

Id0 = 1+sparse(m_x*m_y,1);

Ze0 = sparse(m_x*m_y,1);

%Block matrix for extracting one specifik component
e1 = [Id,Ze,Ze,Ze,Ze];
e2 = [Ze,Id,Ze,Ze,Ze];
e3 = [Ze,Ze,Id,Ze,Ze];
e4 = [Ze,Ze,Ze,Id,Ze];
e5 = [Ze,Ze,Ze,Ze,Id];


%Defining gridspace
x = 0 : hx : W_x;
y = 0 : hy : W_y;


[X, Y] = ndgrid(x, y);

%defining gausspoint IC
Z = exp(-0.7*((X-x0).^2+(Y-y0).^2));
v_x0 = sparse(Z(:));
v_y0 = sparse(Z(:));

%Defining IC for sigmas
sigma_xx0 = sparse(m_x*m_y,1);
sigma_xy0 = sparse(m_x*m_y,1);
sigma_yy0 = sparse(m_x*m_y,1);

%Contstructing u for IC
u = e1'*v_x0+e2'*v_y0+e3'*sigma_xx0+e4'*sigma_xy0+e5'*sigma_yy0;



%Constructing SBP operators
[H_x, HI_x, Dp_x, Dm_x, e1_x, em_x] = d1_upwind_3(m_x, hx);
[H_y, HI_y, Dp_y, Dm_y, e1_y, em_y] = d1_upwind_3(m_y, hy);

%
Dp_y2 = kron(Dp_y,speye(m_x) );
Dm_y2 = kron(Dm_y, speye(m_x));
Dp_x2 = kron(speye(m_y), Dp_x);
Dm_x2 = kron(speye(m_y),Dm_x );

C_p = 1500;
C_s = 800;
rho = 2;     %Density

%Material 'constants', will vary by function
% mu = rho*C_s^2;      %Second lamé parameter
% lambda = rho*C_p^2-2*mu;  %First lamé parameter
mu = 1;
lambda = 1;




[LA, MU, RH] = materialValues(roomX,roomY,roomLength,roomHeight,m_x,m_y,X,Y,lambda,mu,rho);



% LA = spdiags(lambda*ones(m_x*m_y), 0, m_x*m_y, m_x*m_y);
% MU = spdiags(lambda*ones(m_x*m_y), 0, m_x*m_y, m_x*m_y);
% RH = spdiags(lambda*ones(m_x*m_y), 0, m_x*m_y, m_x*m_y);

%Defining matrixes as blockmatrix for equation 11 in report
C = [RH, Ze,   Ze, Ze, Ze;
     Ze,   RH, Ze, Ze, Ze;
     Ze,   Ze,  ((LA+2*MU)/(4*MU))/(LA+MU), Ze, ((-LA)/(4*MU))/(LA+MU);
     Ze,   Ze,   Ze, inv(MU), Ze;
     Ze,   Ze,   ((-LA)/(4*MU))/(LA+MU), Ze, ((LA+2*MU)/(4*MU))/(LA+MU)];

A = [Ze,           Ze, Dp_x2, Ze, Ze;
     Ze,           Ze, Ze, Dp_x2, Ze;
     Dm_x2, Ze, Ze, Ze, Ze;
     Ze,           Dm_x2, Ze, Ze, Ze;
     Ze,      Ze, Ze, Ze, Ze];

B = [Ze, Ze,          Ze, Dp_y2, Ze;
     Ze, Ze,          Ze, Ze, Dp_y2;
     Ze, Ze,     Ze, Ze, Ze;
     Dm_y2, Ze,         Ze, Ze, Ze;
     Ze, Dm_y2, Ze, Ze, Ze];


%Defining and creating BC for west, east, south and north boundary
[L_w, L_e, L_s, L_n] = boundaryNeumann(e1, e2, e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y);

%[L_w, L_e, L_s, L_n] = boundaryAbsorbing(e1, e2, e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y, C_p, C_s, rho);


%Defining Characteristic Boundary operator
L_x = [L_w; L_e];
L_y = [L_s; L_n];


%Discrete inner product
H = kron(H_y,H_x);
H5 = kron( speye(5), H);
HII = H5*C;

I = speye(5*m_x*m_y);

L = [L_w; L_e; L_s; L_n];

P = I - (HII\L')* ((L* (HII\ L'))\ L);

%Renaming A and B to be easier read
D_x = A; 
D_y = B; 

M = C\(D_x+D_y);

%Projection to hopefully make sure our BC can be handled
OP = P*M*P;

%%% Collecting inputs from users regarding what needs to be plot %%%

%Eigenvalues plot or not
eig_answer = upper(input ("Do you want to plot the eigenvalues of h*M (Y/N)?  ", "s"));

%What parameters user wants to plot
plot_parameters = ["V_x", "V_y", "V", "Sigma_xx", "Sigma_xy", "Sigma_xx"];
for i = 1:length(plot_parameters)
    disp([num2str(i) ': ' plot_parameters(i)])
end

plot_answer = int8(input ("Which parameter do you wanna plot from above choices? Enter number from 1 to 6 "));

%Needed simulation defintions
T = input("Input wanted length of simulation (integer): ");
if ~( isnumeric(T) && isscalar(T) )
    T = 10;
end
dt = input("Input wanted timestep: ");
if ~( isnumeric(dt) && isscalar(dt) )
    dt = 0.001;
end

%Plotting eigenvalues 
if eig_answer == "Y"
    ee=eig(full(H5*C*OP));
    plot(real(ee),imag(ee),'*','MarkerSize',8);
    
    % CHANGE TITLE FOR UNDERSTANDING
    title('Eigenvalues to h\cdot B')
    
    xlabel('\Re(h\cdotB)');ylabel('\Im(h\cdotB)');
    ax = gca;          % current axes
    ax.FontSize = 10;
else
    disp("No eigenplot wanted")
end


%Calculating IC data
if plot_answer == 1
    plotData = reshape(e1*u, m_x, m_y);

elseif plot_answer == 2
    plotData = reshape(e2*u, m_x, m_y);

elseif plot_answer == 3
    plotDatay = reshape(e2*u, m_x, m_y);
    plotDatax = reshape(e1*u, m_x, m_y);
    plotData = sqrt(plotDatax.^2 + plotDatay.^2);

elseif plot_answer == 4
%Change for Sigma_xx

elseif plot_answer == 5
%Change for Sigma_xy

elseif plot_answer == 6
%Change for Sigma_yy

else
    plot_answer = 3;
    plotDatay = reshape(e2*u, m_x, m_y);
    plotDatax = reshape(e1*u, m_x, m_y);
    plotData = sqrt(plotDatax.^2 + plotDatay.^2);
end


%Starting values
t = 0;
count = 0;

%IC plot
figure;
hSurf = surf(X,Y,plotData);


%Style changes for graph
colorbar;
%clim([0 0.14])  
axis tight;
xlim([0 W_x])
ylim([0 W_y])
zlim([-0.3 max(e3*u)])
xlabel('x'); ylabel('y'); zlabel(plot_parameters(plot_answer))

axis tight
daspect([1 1 1])

while t < T
    u_next = RK4(u,dt,OP);
    t = dt + t;
    count = count +1;

    %Plotting what we want
    if plot_answer == 1
        plotData = reshape(e1*u, m_x, m_y);
    
    elseif plot_answer == 2
        plotData = reshape(e2*u, m_x, m_y);
    
    elseif plot_answer == 3
        plotDatay = reshape(e2*u, m_x, m_y);
        plotDatax = reshape(e1*u, m_x, m_y);
        plotData = sqrt(plotDatax.^2 + plotDatay.^2);
    
    elseif plot_answer == 4
        %Change for Sigma_xx
    
    elseif plot_answer == 5
        %Change for Sigma_xy
    
    elseif plot_answer == 6
        %Change for Sigma_yy
    end

    if mod(count,15) == 0
        set(hSurf, 'ZData', plotData);
        drawnow;
        disp(["displaying " count])
    end
    u = u_next;
end

%Vad kan vara fel?
% - disp med input vad man vill plotta
% - Snygga till koden
% - Få in ett nytt rum i domänen
% - Testa olika komponenter till 0, är nån stabil? Mental sjukhus?


