function u_next = RK4 (u, dt, B, t, f, frequency)
    
    g = @(t) sin(2*pi*frequency*t);

    % RK4 calculations
    k1 = B * u + g(t)*f;
    k2 = B * (u + 0.5 * dt * k1) + g(t + 0.5 * dt) * f;
    k3 = B * (u + 0.5 * dt * k2) + g(t + 0.5 * dt) * f;
    k4 = B * (u + dt * k3) + g(t + dt) * f;


    %updating the next step
    u_next = u + dt * (k1 + 2*k2 + 2*k3 + k4) / 6;

end