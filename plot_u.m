function [plotData] = plot_u(plot_answer, u, m_x, m_y)
    Id = speye(m_x*m_y);
    Ze = sparse(m_x*m_y,m_x*m_y);
    
    e1 = [Id,Ze,Ze,Ze,Ze];
    e2 = [Ze,Id,Ze,Ze,Ze];
    e3 = [Ze,Ze,Id,Ze,Ze];
    e4 = [Ze,Ze,Ze,Id,Ze];
    e5 = [Ze,Ze,Ze,Ze,Id];
    
    if plot_answer == 1
        plotData = reshape(e1*u, m_x, m_y);
    
    elseif plot_answer == 2
        plotData = reshape(e2*u, m_x, m_y);
    
    elseif plot_answer == 3
        plotDatay = reshape(e2*u, m_x, m_y);
        plotDatax = reshape(e1*u, m_x, m_y);
        plotData = sqrt(plotDatax.^2 + plotDatay.^2);
    
    elseif plot_answer == 4
        plotData = reshape(e3*u, m_x, m_y);
    
    elseif plot_answer == 5
        plotData = reshape(e4*u, m_x, m_y);
    
    elseif plot_answer == 6
        plotData = reshape(e5*u, m_x, m_y);

    elseif plot_answer == 7
        plotDatax = reshape(e3*u, m_x, m_y);
        plotDatay = reshape(e5*u, m_x, m_y);
        plotData = 0.5*(plotDatax + plotDatay);
    end
end