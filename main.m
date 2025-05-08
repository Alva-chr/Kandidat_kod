paramFuncs = {@testing};

for k = 1:numel(paramFuncs)
    %taking from inputfile
    [frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T , snapshot, snapshotName, movie,C_p,C_s] = paramFuncs{k}();

    %Simulates based on inputfile
    simulation(frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, snapshot, snapshotName, movie,C_p,C_s);
end

