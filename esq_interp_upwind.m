function [Ap,Anb_uab] = esq_interp_upwind(Fu_w,u_W,Ap,Anb_uab,Fu_x)
%Upwind interpolation squeme
if Fu_w>0

    U_up=u_W;
    %r_z=-1 + 2*(dot(grad_u_p,rp_w))/(U_do-U_up);
    %yf_r=0;
    Ap=[Ap 0];
    Anb_uab =[Anb_uab Fu_x*U_up];

else 
    %r_z=-1 + 2*(dot(grad_u_p,rp_w))/(U_do-U_up);
    %yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
    Ap=[Ap -Fu_x];
    Anb_uab =[Anb_uab 0];

end

end