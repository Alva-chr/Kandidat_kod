function [L_w, L_e, L_s, L_n] = BC_fun(e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y)

KxW = kron(e1_x', speye(m_y));  
KxE = kron(em_x', speye(m_y));  
KyS = kron(speye(m_x), e1_y');   
KyN = kron(speye(m_x), em_y');  

L_w = [ kron(e3, KxW);
        kron(e4, KxW) ];
L_e = [ kron(e3, KxE);
        kron(e4, KxE) ];

L_s = [ kron(e4, KyS);
        kron(e5, KyS) ];
L_n = [ kron(e4, KyN);
        kron(e5, KyN) ];
end