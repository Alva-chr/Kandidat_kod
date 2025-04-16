function u_next = RK4 (u, dt, B)
    % RK4 calculations
    k1 = dt * B * u;
    k2 = dt * B * (u + 0.5 * k1);
    k3 = dt * B * (u + 0.5 * k2);
    k4 = dt * B * (u + k3);
    
    %updating the next step
    u_next = u + (k1 + 2*k2 + 2*k3 + k4) / 6;

end