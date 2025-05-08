function [LA, MU, RH, inside] = materialValues(roomX,roomY,roomLength,roomHeight,m_x,m_y,X,Y,lambda,mu,rho)

domain = ones(m_x,m_y);
LA = lambda*domain;
MU = mu*domain;
RH = rho*domain;

inside = (X >= roomX) & (X <= roomX+roomLength) & (Y >= roomY) & (Y <= roomY+roomHeight);

LA(inside) = 14.632;
MU(inside) = 0.01;
RH(inside) = 1.293;


LA = spdiags(LA(:),0,m_x*m_y,m_x*m_y);



MU = spdiags(MU(:),0,m_x*m_y,m_x*m_y);


RH = spdiags(RH(:),0,m_x*m_y,m_x*m_y);





end