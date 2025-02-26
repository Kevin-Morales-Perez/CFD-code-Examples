function [U_i_j, d_Ij] = x_Momentum_U(x, y, u_N, u_E, u_S, u_W, u_P, u_N_w, u_N_e, u_S_w, u_S_e, v_n_w, v_n_e, v_s_w, v_s_e, rho, nu, p_P, p_E, p_P_n, p_E_n, p_P_s, p_E_s, dx)
% Solution of the momentum equation in X
% Integration of the momentum equation by finite volumes 
% Control volume ui,j
% Coordinates of points 1,2,3,4,5,6 cells pi,j, pi,j+1  Entry x,y
% Adjacent velocities, Entry
% Component U
% u_N Ui-1,j
% u_E Ui,j+1
% u_S Ui+1,j
% u_W Ui,j-1
% u_P Ui,j to calculate the gradient in P
% Additional U velocities to calculate non-orthogonal diffusion terms,
% Entry 
% u_N_w Ui-1,j-1
% u_N_e Ui-1,j+1
% u_S_w Ui+1,j-1
% u_S_e Ui+1,j-1
% Component V, entry
% v_n_w Vi-1,j
% v_n_e Vi-1,j+1
% v_s_w Vi,j
% v_s_e Vi,j+1
% Density, input rho;
% Viscosity, input nu;
% Adjacent cell pressures, Input 
% p_P Pij
% p_E Pi,j+1
% p_P_n Pi-1,j
% p_E_n Pi-1,j+1
% p_P_s Pi+1,j
% p_E_s Pi+1,j+1
% Step in X, Input, Constant dx;

% Tangents angles 1,2,3,4
tan_u = [];
tan_u(1) = (y(5) - y(6)) / dx;
tan_u(2) = (y(2) - y(1)) / dx;
tan_u(3) = (y(3) - y(2)) / dx;
tan_u(4) = (y(4) - y(5)) / dx;

% Length of face W cell i,j-1
dy_1 = y(1) - y(6);
% Length of face E cell i,j-1
dy_2 = y(2) - y(5);
% Length of face E cell i,j
dy_3 = y(3) - y(4);

% Control volume coordinates for Uij A, B, C, D
vc_x = [];
vc_y = [];

x_c_1 = dy_1 * dx / (dy_2 + dy_1); % Centroid of cell i,j-1
x_c_2 = dy_2 * dx / (dy_3 + dy_2); % Centroid of cell i,j
vc_x(1) = x(1) + x_c_1;
vc_x(2) = vc_x(1);
vc_x(3) = x(2) + x_c_2;
vc_x(4) = vc_x(3);

vc_y(1) = y(6) + tan_u(1) * x_c_1;
vc_y(2) = y(1) + tan_u(2) * x_c_1;
vc_y(3) = y(2) + tan_u(3) * x_c_2;
vc_y(4) = y(5) + tan_u(4) * x_c_2;

% Centroid coordinates 
c_x = [];
c_y = [];
c_x(1) = vc_x(1);
c_x(2) = vc_x(3);
c_y(1) = y(6) + (dy_2 + (y(5) - y(6))) * dy_1 / (dy_2 + dy_1);
c_y(2) = y(5) + (dy_3 + y(4) - y(5)) * dy_2 / (dy_2 + dy_3);

% Coordinates of Ui,j-1, Ui,j, Ui,j+1
cU_x = [];
cU_y = [];
cU_x(1) = x(6);
cU_y(1) = y(6) + 0.5 * dy_1;
cU_x(2) = x(5);
cU_y(2) = y(5) + 0.5 * dy_2;
cU_x(3) = x(4);
cU_y(3) = y(4) + 0.5 * dy_3;

%----------Initializing coefficients Ap - Anb - Scd ---------

Ap = []; % Coefficient that multiplies Up, left side of the equation, Outlet

Anb_uab = []; % Product of coefficients Anb and velocities Ap

S_dc = []; % Non-orthogonal diffusion terms

%----------Diffusion Terms---------
%-----Face W -----
% Tangent and cosine of the angle normal to face W and the line that connects the closest U velocities
tanW_n = (cU_y(1) - cU_y(2)) / (cU_x(1) - cU_x(2));
cosW_n = 1 / sqrt(tanW_n^2 + 1);
% Distance between velocities uij, ui,j-1, delta_psi
delta_ep_w = sqrt((cU_x(2) - cU_x(1))^2 + (cU_y(2) - cU_y(1))^2);
% Length of face W, Aw, delta_n
dA_w = vc_y(2) - vc_y(1);

% Non-orthogonal diffusion term 
s_CD_w = nu * tanW_n * 0.25 * (u_N_w + u_N - u_S_w - u_S);

% Diffusion terms Pressure Gradient P_Ai
D_w = nu * dA_w / (cosW_n * delta_ep_w);

Ap(1) = D_w;
Anb_uab(1) = D_w * u_W;
S_dc(1) = s_CD_w;

%----- Face N -----
% Angles 
tanN_n = -(vc_y(3) - vc_y(2)) / (vc_x(3) - vc_x(2));
cosN_n = 1 / sqrt(tanN_n^2 + 1);
% Distance between velocities ui-1,j, uij, delta_psi
delta_ep_n = dy_2;
% Length of face N, An, delta_n
dA_n = sqrt((vc_y(3) - vc_y(2))^2 + (vc_x(3) - vc_x(2))^2);

% Non-orthogonal diffusion term 
s_CD_n = nu * tanN_n * 0.25 * (u_N_e + u_E - u_N_w - u_W);

% Diffusion terms Pressure Gradient P_Ai
D_n = nu * dA_n / (cosN_n * delta_ep_n);

Ap(2) = D_n;
Anb_uab(2) = D_n * u_N;
S_dc(2) = s_CD_n;

%----- Face E ------
% Angles 
tanE_n = (cU_y(3) - cU_y(2)) / (cU_x(3) - cU_x(2));
cosE_n = 1 / sqrt(tanE_n^2 + 1);
% Distance between velocities, Ui,j+1, Uij, delta_psi
delta_ep_E = sqrt((cU_y(3) - cU_y(2))^2 + (cU_x(3) - cU_x(2))^2);
% Length of face E, delta_n
dA_e = vc_y(3) - vc_y(4);

% Non-orthogonal diffusion term 
s_CD_e = nu * tanE_n * 0.25 * (u_S_e + u_S - u_N_e - u_N);

% Diffusion terms Pressure Gradient P_Ai
D_e = nu * dA_e / (cosE_n * delta_ep_E);

Ap(3) = D_e;
Anb_uab(3) = D_e * u_E;
S_dc(3) = s_CD_e;

% Face S
% Angles
tanS_n = -(vc_y(4) - vc_y(1)) / (vc_x(4) - vc_x(1));
cosS_n = 1 / sqrt(tanS_n^2 + 1);
% Distance between velocities, ui,j, Ui+1,j, delta_psi
delta_ep_S = dy_2;
% Length of face S, delta_n
dA_s = sqrt((vc_y(4) - vc_y(1))^2 + (vc_x(4) - vc_x(1))^2);

% Non-orthogonal diffusion term 
s_CD_s = nu * tanS_n * 0.25 * (u_S_w + u_W - u_S_e - u_E);

% Diffusion terms Pressure Gradient P_Ai
D_s = nu * dA_s / (cosS_n * delta_ep_S);

Ap(4) = D_s;
Anb_uab(4) = D_s * u_S;
S_dc(4) = s_CD_s;

%----------------------------- ----------

%--------Convection terms----------


%Mass Flux velocity U

%Face W
%Upwind 
Fu_w=0.5*(u_W+u_P)*dA_w;
%rp_w=[(cU_x(1)-cU_x(2)) (cU_y(1)-cU_y(2))];%Vector que une los nodos , P y W en vel_U

%if Fu_w>0

%    U_up=u_W;
%    U_do=u_P;
%    r_z=-1 + 2*(dot(grad_u_p,rp_w))/(U_do-U_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap Fu_w*yf_r*0.5];
%    Anb_uab =[Anb_uab Fu_w*U_up*(1-0.5*yf_r)];

%else 
%    U_up=u_P;
%    U_do=u_W;
%    r_z=-1 + 2*(dot(grad_u_p,rp_w))/(U_do-U_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap -Fu_w*(1-0.5*yf_r)];
%    Anb_uab =[Anb_uab -Fu_w*U_do*yf_r*0.5];0

%end

%[Ap,Anb_uab]=esq_interp(Fu_w,u_W,u_P,grad_u_p,rp_w,Ap,Anb_uab,Fu_w);
[Ap,Anb_uab]=esq_interp_upwind(Fu_w,u_W,Ap,Anb_uab,Fu_w);

%u_w=U_up +0.5*yf_r*(U_do -U_up);

%Face  N
%Diferencias centrales
Fu_n =0.5*(u_N + u_P)*cosN_n*tanN_n*rho;

%Face  E
%Upwind 
Fu_e = -0.5*(u_P + u_E)*dA_e*rho;
%rp_e=[(cU_x(3)-cU_x(2)) (cU_y(3)-cU_y(2))];%Vector que une los nodos , P y E en vel_U

%if Fu_e>0
%    U_up=u_E;
%    U_do=u_P;
%    r_z=-1 + 2*(dot(grad_u_p,rp_e))/(U_do-U_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap Fu_e*yf_r*0.5];
%    Anb_uab =[Anb_uab Fu_e*U_up*(1-0.5*yf_r)];

%else
%    U_up=u_P;
%    U_do=u_E;
%    r_z=-1 + 2*(dot(grad_u_p,rp_e))/(U_do-U_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap -Fu_e*(1-0.5*yf_r)];
%    Anb_uab =[Anb_uab -Fu_e*U_do*yf_r*0.5];

%end

%[Ap,Anb_uab]=esq_interp(Fu_e,u_E,u_P,grad_u_p,rp_e,Ap,Anb_uab,Fu_e);
[Ap,Anb_uab]=esq_interp_upwind(Fu_e,u_E,Ap,Anb_uab,Fu_e);

%u_e=U_up + 0.5*yf_r*(U_do -U_up);

%Face S
Fu_s =-0.5*(u_P + u_S)*dA_s*cosS_n*tanS_n*rho;

%Mass Flux for  V
%Face  W
%Fv_w=0 velocity vector perpendicular to face W

%Face N
%Upwind TVD
Fv_n = -0.5*(v_n_w*v_n_e)*dA_n*cosN_n*rho;
%rp_e=[0 dy_2];%Vector que une los nodos  P y N en vel_v

%if Fv_n>0
%    U_up=u_N;
%    U_do=u_P;
%    r_z=-1 + 2*(dot(grad_u_p,rp_e))/(U_do-U_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap Fv_n*yf_r*0.5];
%    Anb_uab =[Anb_uab Fv_n*U_up*(1-0.5*yf_r)];
%[Ap,Anb_uab]=esq_interp(Fv_n,u_N,u_P,grad_u_p,rp_e,Ap,Anb_uab,Fv_n);
[Ap,Anb_uab]=esq_interp_upwind(Fv_n,u_N,Ap,Anb_uab,Fv_n);

    %Face N , U component
%    Ap=[Ap Fu_n*yf_r*0.5];
%    Anb_uab =[Anb_uab Fu_n*U_up*(1-0.5*yf_r)];
%[Ap,Anb_uab]=esq_interp(Fv_n,u_N,u_P,grad_u_p,rp_e,Ap,Anb_uab,Fu_n);
[Ap,Anb_uab]=esq_interp_upwind(Fv_n,u_N,Ap,Anb_uab,Fu_n);

%else
%   U_up=u_P;
%   U_do=u_N;
%    r_z=-1 + 2*(dot(grad_u_p,rp_e))/(U_do-U_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap -Fv_n*(1-0.5*yf_r)];
%    Anb_uab =[Anb_uab -Fv_n*U_do*yf_r*0.5];
    %Contribución de la cara N ,flujo másico en U
%    Ap=[Ap -Fu_n*(1-0.5*yf_r)];
%    Anb_uab =[Anb_uab -Fu_n*U_do*yf_r*0.5];


%end
%u_n=U_up +0.5*yf_r*(U_do -U_up);
 
%Cara E
%Fv_e = 0 Vector velocidad perpendicular a la cara W

%Face S
%Upwind TVD
Fv_s= 0.5*(v_s_w*v_s_e)*dA_n*cosS_n*rho;
%rp_s=[0 -dy_2];%Vector que une los nodos  P y N en vel_v
%if Fv_s>0
%    U_up=u_S;
%    U_do=u_P;
%    r_z=-1 + 2*(dot(grad_u_p,rp_s))/(U_do-U_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap Fv_s*yf_r*0.5];
%    Anb_uab =[Anb_uab Fv_s*U_up*(1-0.5*yf_r)];
    %Cara S , contribución de flujo masico de U
%[Ap,Anb_uab]=esq_interp(Fv_s,u_S,u_P,grad_u_p,rp_s,Ap,Anb_uab,Fv_s);
[Ap,Anb_uab]=esq_interp_upwind(Fv_s,u_S,Ap,Anb_uab,Fv_s);

    %Ap=[Ap Fu_s*yf_r*0.5];
    %Anb_uab =[Anb_uab Fu_s*U_up*(1-0.5*yf_r)];
%[Ap,Anb_uab]=esq_interp(Fv_s,u_S,u_P,grad_u_p,rp_s,Ap,Anb_uab,Fu_s);
[Ap,Anb_uab]=esq_interp_upwind(Fv_s,u_S,Ap,Anb_uab,Fu_s);

%else
%    U_up=u_P;
%    U_do=u_S;
%    r_z=-1 + 2*(dot(grad_u_p,rp_s))/(U_do-U_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap -Fv_s*(1-0.5*yf_r)];
%    Anb_uab =[Anb_uab -Fv_s*U_do*yf_r*0.5];
    %Contribución de la cara S ,flujo másico en U
%    Ap=[Ap -Fu_s*(1-0.5*yf_r)];
%    Anb_uab =[Anb_uab -Fu_s*U_do*yf_r*0.5];

%end
%u_s=Us_up +0.5*yf_r*(U_do -U_up);


%DP/DX 
dp_dx1=(p_E-p_P)/sqrt((c_x(1)-c_x(2))^2 + (c_y(1)-c_y(2))^2); %pressure gradient , line that joins cell nodes
dp_dy=0.25*(p_P_n + p_E_n -p_P_s - p_E_s)/dy_2;%pressure gradient respect to Y
tan_cen_p = (c_y(2)-c_y(1))/(c_x(2)-c_x(1)); % tangent to the angle , line that joins pressure cell nodes
theta_cent = atan(tan_cen_p); %interior line angle
dp_dx= dp_dx1*cos(theta_cent) - (dp_dy/cos(theta_cent) - dp_dx1*tan_cen_p)*sin(theta_cent);
delta_V = 0.5*(dA_e + dA_w)*dx; %vol of the cell 



%--------    Final Equation    ----------

% Final coefficient AP, left side of the equality;
Ap = Ap(~isnan(Ap));
Anb_uab = Anb_uab(~isnan(Anb_uab));
S_dc = S_dc(~isnan(S_dc));

A_P_i = sum(Ap); % Coefficient needed to calculate U_p and pressure correction;
d_Ij = -delta_V / (A_P_i * dx);
U_i_j = (1 / A_P_i) * (sum(Anb_uab) + sum(S_dc) - dp_dx * delta_V); % New U_P, main output
end
