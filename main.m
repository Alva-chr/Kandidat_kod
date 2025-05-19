paramFuncs = {@simx_17,@simx_18,@simx_19,@simx_20,@simx_21,@simx_22,@simx_23,@simx_24,@simx_25,@simx_26,@simx_27,@simx_28,@simx_29,@simx_30};
%paramFuncs = {@simy_15};



for k = 1:numel(paramFuncs)
    %taking from inputfile
    [frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T , snapshot, snapshotName, movie,C_p,C_s, decibel_txt,save_movie] = paramFuncs{k}();

    %Simulates based on inputfile
    simulation(frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, snapshot, snapshotName, movie,C_p,C_s,decibel_txt, save_movie);
end

