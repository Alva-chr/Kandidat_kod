function [L_w, L_e, L_s, L_n] = boundaryAbsorbingNeumann(e1, e2, e3, e4, e5, m_x, m_y, e1_x, em_x, e1_y, em_y, C_p, C_s, rho)

Im_x = speye(m_x);
Im_y = speye(m_y);

Z_s = rho*C_s;
Z_p = rho*C_p;

e_S = kron(e1_y', Im_x);  
e_N = kron(em_y', Im_x);  
e_W = kron(Im_y, e1_x');   
e_E = kron(Im_y, em_x');  

e_W_no_corner = e_W(2:end-1,:);
e_E_no_corner = e_E(2:end-1,:);

L_w = [ e_W*(e3-Z_p*e1);
        e_W_no_corner*(-e4+Z_s*e2)];
L_e = [ e_E*(e3+Z_p*e1);
        e_E_no_corner*(e4+Z_s*e2)];

L_s = [ e_S*(e5-Z_p*e2);
        e_S*(-e4+Z_p*e1) ];
L_n = [ e_N*e5;
        e_N*e4];
end