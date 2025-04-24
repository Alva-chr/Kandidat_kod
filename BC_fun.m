function [L_w, L_e, L_s, L_n] = BC_fun(e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y)

Im_x = speye(m_x);
Im_y = speye(m_y);



NxW = kron(e1_x', Im_y);  
NxE = kron(em_x', Im_y);  
NyS = kron(Im_x, e1_y');   
NyN = kron(Im_x, em_y');  

L_w = [ e3*NxW;
        e4*NxW ];
L_e = [ e3*NxE;
        e4*NxE ];

L_s = [ e4*NyS;
        e5*NyS ];
L_n = [ e4*NyN;
        e5*NyN ];
end