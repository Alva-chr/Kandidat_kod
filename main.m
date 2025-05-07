paramFuncs = {@input1};

for k = 1:numel(paramFuncs)
    % call the k-th input-function
    [frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, dt ] = paramFuncs{k}();

    % now feed all of those into your solver
    simulation(frequency, x0, y0, m_x, m_y,W_x, W_y, roomLength, roomHeight, roomX, roomY, mu, lambda,rho, BC, eig_answer, plot_answer, T, dt);
end

