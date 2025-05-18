paramFuncs = {@simy_1,@simy_2,@simy_3,@simy_4,@simy_5,@simy_6,@simy_7,@simy_8,@simy_9,@simy_10,@simy_11,@simy_12,@simy_13,@simy_14,@simy_16,@simy_17,@simy_18,@simy_19,@simy_20,@simy_21,@simy_22,@simy_23,@simy_24,@simy_25,@simy_26,@simy_27,@simy_28,@simy_29,@simy_30,@simx_3,@simx_5,@simx_7,@simx_9,@simx_11,@simx_13,@simx_15,@simx_17,@simx_19,@simx_21,@simx_23,@simx_25,@simx_27,@simx_29};
%paramFuncs = {@simy_15};



for k = 1:numel(paramFuncs)
    %taking from inputfile
    [frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T , snapshot, snapshotName, movie,C_p,C_s, decibel_txt,save_movie] = paramFuncs{k}();

    %Simulates based on inputfile
    simulation(frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, snapshot, snapshotName, movie,C_p,C_s,decibel_txt, save_movie);
end

