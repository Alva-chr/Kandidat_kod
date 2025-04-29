function [LA, MU, RH] = materialValues(roomX,roomY,roomLength,roomHeight,m_x,m_y,X,Y,lambda,mu,rho);

domain = ones(m_x,m_y);


inside = (X >= roomX) & (X <= roomX+roomLength) & (Y >= roomY) & (Y <= roomY+roomHeight);

domain(inside) = 0.0001;

LA = lambda*domain;

LA = spdiags(LA(:),0,m_x*m_y,m_x*m_y);


MU = mu*domain;
MU = spdiags(MU(:),0,m_x*m_y,m_x*m_y);

RH = rho*domain;
RH = spdiags(RH(:),0,m_x*m_y,m_x*m_y);





end