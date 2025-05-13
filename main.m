paramFuncs = {@verifyWx1_30,@verifyWx1_60,@verifyWx1_90,@verifyWx1_120,@verifyWx1_150};



for k = 1:numel(paramFuncs)
    %taking from inputfile
    [frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T , snapshot, snapshotName, movie,C_p,C_s, decibel_txt,save_movie] = paramFuncs{k}();

    %Simulates based on inputfile
    simulation(frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, snapshot, snapshotName, movie,C_p,C_s,decibel_txt, save_movie);
end

