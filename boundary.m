function [L_w, L_e, L_s, L_n] = boundary(e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y)
%BOUNDARY  Construct SAT boundary operators for traction-free BCs
%   [L_w, L_e, L_s, L_n] = boundary(m_x, m_y, e1_x, em_x, e1_y, em_y)
%   returns the west, east, south, and north boundary operators
%   for a 2D SBP grid of size m_x-by-m_y.
%   e1_x, em_x: m_x-by-1 selectors for x=1 and x=m_x
%   e1_y, em_y: m_y-by-1 selectors for y=1 and y=m_y

% build 1D-to-2D pick-out rows:
KxW = kron(e1_x', speye(m_y));   % selects all y at x=1 (west face)
KxE = kron(em_x', speye(m_y));   % selects all y at x=m_x (east face)
KyS = kron(speye(m_x), e1_y');   % selects all x at y=1 (south face)
KyN = kron(speye(m_x), em_y');   % selects all x at y=m_y (north face)

L_w = [ kron(e3, KxW);
        kron(e4, KxW) ];
L_e = [ kron(e3, KxE);
        kron(e4, KxE) ];

L_s = [ kron(e4, KyS);
        kron(e5, KyS) ];
L_n = [ kron(e4, KyN);
        kron(e5, KyN) ];
end