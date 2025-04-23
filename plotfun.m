function plotfun(u,t,plt,txt)
    plt.ZData = u;
    txt.String = "$u(x,y,t ="+ t + ")$";
    drawnow;
end