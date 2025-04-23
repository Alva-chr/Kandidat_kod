function [L_w, L_e, L_s, L_n] = BC_fun(e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y)

NxW = kron(e1_x', speye(m_y));  
NxE = kron(em_x', speye(m_y));  
NyS = kron(speye(m_x), e1_y');   
NyN = kron(speye(m_x), em_y');  

L_w = [ kron(e3, NxW);
        kron(e4, NxW) ];
L_e = [ kron(e3, NxE);
        kron(e4, NxE) ];

L_s = [ kron(e4, NyS);
        kron(e5, NyS) ];
L_n = [ kron(e4, NyN);
        kron(e5, NyN) ];
end