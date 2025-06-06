function [L_w, L_e, L_s, L_n] = boundaryNeumann(e1, e2, e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y)

Im_x = speye(m_x);
Im_y = speye(m_y);

e_S = kron(e1_y', Im_x);  
e_N = kron(em_y', Im_x);  
e_W = kron(Im_y, e1_x');   
e_E = kron(Im_y, em_x');  

e_W_no_corner = e_W(2:end-1,:);
e_E_no_corner = e_E(2:end-1,:);

L_w = [ e_W*e3;
        e_W_no_corner*e4 ];
L_e = [ e_E*e3;
        e_E_no_corner*e4];

L_s = [ e_S*e4;
        e_S*e5 ];
L_n = [ e_N*e4;
        e_N*e5 ];
end