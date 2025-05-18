function [] = simulation(frequency, x0, y0, m_x, m_y, W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, snapshot, snapshotName, movie,C_p,C_s,decibel_txt, save_movie)


if W_y/m_y > 2
    m_y = 2*W_y;
end

if W_x/m_x > 2
    m_x = 2*W_x;
end



m_x = ceil(m_x);
m_y = ceil(m_y);

T = sqrt(W_y^2+W_x^2)/C_s+2/frequency;

%gridspacing for height and length
hx = W_x/(m_x-1);
hy = W_y/(m_y-1);


Id = speye(m_x*m_y);

Ze = sparse(m_x*m_y,m_x*m_y);

%Block matrix for extracting one specifik component
e1 = [Id,Ze,Ze,Ze,Ze];
e2 = [Ze,Id,Ze,Ze,Ze];
e3 = [Ze,Ze,Id,Ze,Ze];
e4 = [Ze,Ze,Ze,Id,Ze];
e5 = [Ze,Ze,Ze,Ze,Id];


%Defining gridspace and grid
x = 0 : hx : W_x;
y = 0 : hy : W_y;
[X, Y] = ndgrid(x, y);



%Contstructing u for IC
u = e1'*sparse(m_x*m_y,1)+e2'*sparse(m_x*m_y,1)+e3'*sparse(m_x*m_y,1)+e4'*sparse(m_x*m_y,1)+e5'*sparse(m_x*m_y,1);
%Constructing SBP operators

[H_x, HI_x, Dp_x, Dm_x, e1_x, em_x] = d1_upwind_3(m_x, hx);
[H_y, HI_y, Dp_y, Dm_y, e1_y, em_y] = d1_upwind_3(m_y, hy);

Dp_y2 = kron(Dp_y,speye(m_x) );
Dm_y2 = kron(Dm_y, speye(m_x));
Dp_x2 = kron(speye(m_y), Dp_x);
Dm_x2 = kron(speye(m_y),Dm_x );



%Creating room
[LA, MU, RH, inside] = materialValues(roomX,roomY,roomLength,roomHeight,m_x,m_y,X,Y,lambda,mu,rho);

assert(isdiag(MU));

%Defining matrixes as blockmatrix for equation 11 in report
C = [RH, Ze,   Ze, Ze, Ze;
     Ze,   RH, Ze, Ze, Ze;
     Ze,   Ze,  ((LA+2*MU)/((4*MU)))/(LA+MU), Ze, ((-LA)/(4*MU))/(LA+MU);
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
if BC == "NBC"
    [L_w, L_e, L_s, L_n] = boundaryNeumann(e1, e2, e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y);
elseif BC == "ABC" 
    [L_w, L_e, L_s, L_n] = boundaryAbsorbing(e1, e2, e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y, C_p, C_s, rho);
else 
    [L_w, L_e, L_s, L_n] = boundaryAbsorbingNeumann(e1, e2, e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y, C_p, C_s, rho);
end

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

dt = 2.5/abs(eigs(OP,1));


dt = T/(ceil(T/dt));


%Setting up point source location
%Findong the nearest grid‐indicess to (p_x, p_y)
[~, ix] = min(abs(x - x0));
[~, iy] = min(abs(y - (W_y-y0)));

%Building a mask with a single one
PointSource = sparse(m_x,m_y);
PointSource(ix, iy) = 1;
PointSource = sparse(PointSource(:));

%Setting up vectors for point sources
%f = e1'*PointSource+e2'*sparse(m_x*m_y,1)+e3'*sparse(m_x*m_y,1)+e4'*sparse(m_x*m_y,1)+e5'*sparse(m_x*m_y,1);
f = e1'*sparse(m_x*m_y,1)+e2'*PointSource+e3'*sparse(m_x*m_y,1)+e4'*sparse(m_x*m_y,1)+e5'*sparse(m_x*m_y,1);
%f = e1'*sparse(m_x*m_y,1)+e2'*sparse(m_x*m_y,1)+e3'*sparse(m_x*m_y,1)+e4'*sparse(m_x*m_y,1)+e5'*sparse(m_x*m_y,1);
f = H5\f;

disp(f)


plot_parameters = ["V_x","V_y","V","Sigma_xx","Sigma_xy","Sigma_xx", "P"];

%Plotting eigenvalues 
if eig_answer == "Y"
    ee=eig(full(H5*C*OP));
    plot(real(ee),imag(ee),'*','MarkerSize',8);
    
    % CHANGE TITLE FOR UNDERSTANDING
    title('Eigenvalues to h\cdot B')
    
    xlabel('\Re(h\cdotB)');ylabel('\Im(h\cdotB)');
    ax = gca;          % current axes
    ax.FontSize = 10;
    
    set(gcf, 'Units','pixels', 'Position', [0 0 1920 1080]);
    takeSnapshot("eig", snapshotName)


else
    disp("No eigenplot wanted")
end


%Calculating IC data
plotData = plot_u(plot_answer,u,m_x,m_y);


%Starting values
t = 0;
count = 0;



m_y = ceil(m_y);
m_x = ceil(m_x);


%Snapchot timings
S1 = (x0+2)/C_p - mod((x0+2)/C_p, dt);
S2 = sqrt((roomY - y0)^2+(roomX-x0)^2)/C_p - mod(sqrt((roomY - y0)^2+(roomX-x0)^2)/C_p, dt);
S3 = T;

%IC plot
figure;
set(gcf, 'Units','pixels', 'Position', [0 0 1920 1080]);
hSurf = surf(X,Y,plotData, "edgecolor", "None");
view(2)

%Style changes for graph
colorbar;
clim([-15 15])  
axis tight;
xlim([0 W_x])
ylim([0 W_y])
zlim([-0.3 max(e3*u)])
xlabel('x'); ylabel('y'); zlabel(plot_parameters(plot_answer))

axis tight
daspect([1 1 1])

meanPressureTime = sparse(m_x,m_y);

if snapshot == "Y" || save_movie == "Y"
    if ~exist(snapshotName, 'dir')
        mkdir(snapshotName);
    end
end

fid = fopen("decibel_saved.txt", "a");

if save_movie == "Y"
    filename = snapshotName +"_movie.mp4";
    fullpath = fullfile(snapshotName, filename);
    v = VideoWriter(fullpath, 'MPEG-4');
    v.Quality = 100;
    open(v)
end

while t < T
    u_next = RK4(u,dt,OP,t,f,frequency);
    t = dt + t;
    count = count +1;

    %Plotting what we want
    plotData = plot_u(plot_answer,u,m_x,m_y);

    if snapshot == "Y"
        if (t <= S1+dt/2 && t >= S1-dt/2)
            set(hSurf, 'ZData', plotData);
            drawnow;

            %Takes snapshot
            takeSnapshot("S1", snapshotName)          

        elseif (t <= S2+dt/2 && t >= S2-dt/2)
            %disp("displaying S2")
            set(hSurf, 'ZData', plotData);
            drawnow;

            takeSnapshot("S2", snapshotName)  

        elseif (t <= S3+dt/2 && t >= S3-dt/2)
            %disp(S3)
            set(hSurf, 'ZData', plotData);
            drawnow;

            takeSnapshot("S3", snapshotName)  
        end
    

    end
    if (T-t)<=2/(f)
        meanPressureTime = meanPressureTime + reshape((0.5*e3*(u)+0.5*e5*(u)).^2,m_x,m_y);
    end



    if mod(count, 15) == 0
        if save_movie == "Y"
            set(hSurf, 'ZData', plotData);
            frame = getframe(gcf);
            writeVideo(v, frame)
        end

        if movie == "Y"
            set(hSurf, 'ZData', plotData);
            drawnow;
            %disp(["displaying " count t])
        end
    end
    

    u = u_next;
end

meanPressureTime = frequency*dt*meanPressureTime/2;

meanPressureTime(~inside) = 0;

meanPressureRoom = sqrt( sum( meanPressureTime(:).^2 )/(roomHeight*roomLength) );

disp(meanPressureRoom)

decibel = 20*log10(meanPressureRoom/(20*10^-6));

disp(decibel)

disp(datetime)

if decibel_txt == "Y" 
    fprintf(fid,'%s\t%f\n', snapshotName, full(decibel));
end

fclose(fid);
end

