function [V_i_j, d_iJ] = y_Momentum_V(x,y,v_N,v_E,v_S,v_W,v_P,v_N_w,v_N_e,v_S_w,v_S_e,u_n_w,u_n_e,u_s_w,u_s_e,rho,nu,p_P,p_S,dx)
% Solution of the momentum equation in Y
% Integration of the momentum equation by finite volumes.
% Control volume for velocities V (Transported variable is V)
% Vij
% Coordinates of points 1,2,3,4,5,6 cells pi,j, pi+1,j  Entry x and y
% Adjacent velocities, Entry
% Component V,
% v_N = Vi-1,j
% v_E = Vi,j+1
% v_S = Vi+1,j
% v_W = Vi,j-1
% v_P = Vi,j to calculate the gradient in P

% Additional velocities of v to calculate non-orthogonal diffusion terms
% v_N_w = Vi-1,j-1
% v_N_e = Vi-1,j+1
% v_S_w = Vi+1,j-1
% v_S_e = Vi+1,j-1

% Component U, entry
% u_n_w = Ui,j-1
% u_n_e = ui,j
% u_s_w = ui+1,j-1
% u_s_e = u+1,j

% Density, input
% rho
% Viscosity, input 
% nu

% Pressures for dp_dy
% p_P = Pij
% p_S = Pi+1,j

%-------- Geometric parameters ------
% Step in X, Input, Constant
% dx

% Tangent angles 1,2,3,4

tan_v=[];
tan_v(1)=(y(2)-y(1))/dx;
tan_v(2) = (y(3)-y(6))/dx;
tan_v(3) = (y(4)-y(5))/dx;

%Length of face W i,j-1
dy_1=y(1)-y(6);
%Length of face E  i,j-1
dy_2=y(2)-y(3);

% Control volume coordinates for Uij A, B, C, D
vc_x=[];
vc_y=[];
vc_x(1)=x(6);
vc_x(2)=x(6);
vc_x(3)=x(3);
vc_x(4)=x(3);

vc_y(1)=y(6) - 0.5*dy_1;
vc_y(2)=y(6) + 0.5*dy_1;
vc_y(3)=y(3) + 0.5*dy_2;
vc_y(4)=y(3) - 0.5*dy_2;

% Centroid of cell i,j-1
x_c_1=dy_1*dx/(dy_2+dy_1);

% Centroid of cell i,j
c_x=[];
c_y =[];
c_x(1)=x_c_1;
c_x(2)=x_c_1;
c_y(1)=vc_y(2) + x_c_1*(vc_y(3)-vc_y(2))/dx;
c_y(2)=vc_y(1) + x_c_1*(vc_y(4)-vc_y(1))/dx;

%coordinates of  Vi-1,j , Vi,j ,Vi,j+1
cV_x=[];
cV_y=[];

cV_x(1)=x_c_1;
cV_x(2)=x_c_1;
cV_x(3)=x_c_1;

cV_y(1)=y(1) + tan_v(1)*x_c_1;
cV_y(2)=y(6) + tan_v(2)*x_c_1;
cV_y(3)=y(5) + tan_v(3)*x_c_1;


%----------Initializing coefficients Ap - Anb - Scd ---------

Ap = []; % Coefficient that multiplies Up, left side of the equation, Outlet

Anb_vab = []; % Product of coefficients Anb and velocities Ap

S_dc = []; % Non-orthogonal diffusion terms

%----------Diffusion Terms---------
%-----Face W -----
%tangente y coseno de angulo normal a cara W y linea que une velocidades U mas
%cercanas
tanW_n=-tan_v(2);
cosW_n=1/sqrt(tanW_n^2 + 1);
% Tangent and cosine of the angle normal to face W and the line that connects the closest U velocities
delta_ep_w= 2*sqrt((cV_y(2)-y(6))^2 + (cV_x(2)-x(6))^2);
% Length of face W, Aw, delta_n
dA_w=dy_1;
% Non-orthogonal diffusion term 
s_CD_w =nu*tanW_n*0.25*(v_N_w+ v_N-v_S_w- v_S);
% Diffusion terms Pressure Gradient P_Ai
D_w = nu*dA_w/(cosW_n*delta_ep_w);

Ap(1)=D_w;
Anb_vab(1)=D_w*v_W;
S_dc(1)=s_CD_w;


%----- face  N -----
%Angles
tanN_n = (vc_y(3)-vc_y(2))/dx;
cosN_n=1/sqrt(tanN_n^2 + 1);
% Distance between velocities ui-1,j, uij, delta_psi
delta_ep_n=cV_y(1)-cV_y(2);
% Length of face N, An, delta_n
dA_n=sqrt((vc_y(3)-vc_y(2))^2 + (vc_x(3)-vc_x(2))^2);

% Non-orthogonal diffusion term 
s_CD_n = nu*tanN_n*0.25*(v_N_e+ v_E-v_N_w- v_W);

% Diffusion terms Pressure Gradient P_Ai
D_n = nu*dA_n/(cosN_n*delta_ep_n);

Ap(2)=D_n;
Anb_vab(2)=D_n*v_N;
S_dc(2)= s_CD_n ;

%----- Face E ------
% Angles  
tanE_n=-tan_v(2);
cosE_n=1/sqrt(tanE_n^2 + 1);
% Distance between velocities, Ui,j+1, Uij, delta_psi
delta_ep_E=2*sqrt((y(3)-cV_y(2))^2 + (x(3)-cV_x(2))^2);
% Length of face E, delta_n
dA_e=vc_y(3)-vc_y(4);

% Non-orthogonal diffusion term 
s_CD_e = nu*tanE_n*0.25*(v_S_e+ v_S-v_N_e- v_N);

% Diffusion terms Pressure Gradient P_Ai
D_e = nu*dA_e/(cosE_n*delta_ep_E);

Ap(3)=D_e;
Anb_vab(3)=D_e*v_E;
S_dc(3)= s_CD_e;

% Face S
% Angles
tanS_n=(vc_y(4)-vc_y(1))/dx;
cosS_n=1/sqrt(tanS_n^2 + 1);
% Distance between velocities, ui,j, Ui+1,j, delta_psi
delta_ep_S=  sqrt((vc_y(3)-vc_y(2))^2+(vc_x(3)-vc_x(2))^2);
% Length of face S, delta_n
dA_s= sqrt((vc_y(4)-vc_y(1))^2 + (vc_x(4)-vc_x(1))^2);

% Non-orthogonal diffusion term 
s_CD_s = nu*tanS_n*0.25*(v_S_w+ v_W-v_S_e- v_E);

% Diffusion terms Pressure Gradient P_Ai
D_s = nu*dA_s/(cosS_n*delta_ep_S);

Ap(4)=D_s;
Anb_vab(4)=D_s*v_S;
S_dc(4)= s_CD_s;

%--------Terminos de convección----------
%Gradiente de V en P(centro del volumen de control)
%diferencias de distancia entre nodo P y nodos adyacentes
%dif_G_P = [2*(x(6)-cV_x(2)) 0 2*(x(3)-cV_x(2)) 0;2*(y(6)-cV_y(2)) cV_y(1)-cV_y(2) 2*(y(3)-cV_y(2)) cV_y(3)-cV_y(2)]' ;
%dif_fhi =[v_N-v_P;v_E-v_P;v_S-v_P;v_W-v_P];
%Gradiente en P pro minimos cuadrados 
%grad_v_p=inv(dif_G_P.'*dif_G_P)*dif_G_P.'*dif_fhi;

%Flujos másicos para velocidad U , variable transportada v

%cara W
Fu_w=0.5*(u_n_w+u_s_w)*dA_w*rho;
%rp_w=[2*(x(6)-cV_x(2)) 2*(y(6)-cV_y(2))];%Vector que une los nodos , P y W en vel_U

%if Fu_w>0

%    V_up=v_W;
%    V_do=v_P;
%    r_z=-1 + 2*(dot(grad_v_p,rp_w))/(V_do-V_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap Fu_w*yf_r*0.5];
%    Anb_vab =[Anb_vab Fu_w*V_up*(1-0.5*yf_r)];

%else 
%    V_up=v_P;
%    V_do=v_W;
%    r_z=-1 + 2*(dot(grad_v_p,rp_w))/(V_do-V_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap -Fu_w*(1-0.5*yf_r)];
%    Anb_vab =[Anb_vab -Fu_w*V_do*yf_r*0.5];

%end

%[Ap,Anb_vab]=esq_interp(Fu_w,v_W,v_P,grad_v_p,rp_w,Ap,Anb_vab,Fu_w);
[Ap,Anb_vab]=esq_interp_upwind(Fu_w,v_W,Ap,Anb_vab,Fu_w);
%Cara N
Fu_n =0.5*(u_n_w+u_n_e)*dA_n*cosN_n*tanN_n*rho;

%Cara E
Fu_e = -0.5*(u_n_e+u_s_e)*dA_e*rho;
%rp_e=[2*(x(3)-cV_x(2)) 2*(y(3)-cV_y(2))];%Vector que une los nodos , P y E en vel_V
%if Fu_e>0
%    V_up=v_E;
%    V_do=v_P;
%    r_z = -1 + 2*(dot(grad_v_p,rp_e))/(V_do-V_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap Fu_e*yf_r*0.5];
%    Anb_vab =[Anb_vab Fu_e*V_up*(1-0.5*yf_r)];

%else
%    V_up=v_P;
%    V_do=v_E;
%    r_z = -1 + 2*(dot(grad_v_p,rp_e))/(V_do-V_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap -Fu_e*(1-0.5*yf_r)];
%    Anb_vab =[Anb_vab -Fu_e*V_do*yf_r*0.5];

%end
%[Ap,Anb_vab]=esq_interp(Fu_e,v_E,v_P,grad_v_p,rp_e,Ap,Anb_vab,Fu_e);
[Ap,Anb_vab]=esq_interp_upwind(Fu_e,v_E,Ap,Anb_vab,Fu_e);

%Cara S
Fu_s = -0.5*(u_s_w+u_s_e)*dA_s*tanS_n*cosS_n*rho;

%Flujos másicos para velocidad V

%Cara W
%Fv_w=0;

%Cara N
Fv_n =-0.5*(v_N + v_P)*dA_n*cosN_n*rho;
%rp_e=[0 cV_y(1)-cV_y(2)];%Vector que une los nodos  P y N en vel_v

%if Fv_n>0
%    V_up=v_N;
%    V_do=v_P;
%    r_z = -1 + 2*(dot(grad_v_p,rp_e))/(V_do-V_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap Fv_n*yf_r*0.5];
%    Anb_vab =[Anb_vab Fv_n*V_up*(1-0.5*yf_r)];
%[Ap,Anb_vab]=esq_interp(Fv_n,v_N,v_P,grad_v_p,rp_e,Ap,Anb_vab,Fv_n);
[Ap,Anb_vab]=esq_interp_upwind(Fv_n,v_N,Ap,Anb_vab,Fv_n);

    %Cara N , contribución de flujo masico de U
%    Ap=[Ap Fu_n*yf_r*0.5];
%    Anb_vab =[Anb_vab Fu_n*V_up*(1-0.5*yf_r)];
%[Ap,Anb_vab]=esq_interp(Fv_n,v_N,v_P,grad_v_p,rp_e,Ap,Anb_vab,Fu_n);
[Ap,Anb_vab]=esq_interp_upwind(Fv_n,v_N,Ap,Anb_vab,Fu_n);
    

%else
%    V_up=v_P;
%    V_do=v_N;
%    r_z = -1 + 2*(dot(grad_v_p,rp_e))/(V_do-V_up);
%    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
%    Ap=[Ap -Fv_n*(1-0.5*yf_r)];
%    Anb_vab =[Anb_vab -Fv_n*V_do*yf_r*0.5];
    %Contribución de la cara N ,flujo másico en u
%    Ap=[Ap -Fu_n*(1-0.5*yf_r)];
%    Anb_vab =[Anb_vab -Fu_n*V_do*yf_r*0.5];


%end

%Cara E
%Fv_E=0;

%Cara S
Fv_s=0.5*(v_S+v_N)*dA_s*cosS_n*rho;
%rp_s=[0 cV_y(3)-cV_y(2)];%Vector que une los nodos  P y N en vel_v
%{
if Fv_s>0
    V_up=v_S;
    V_do=v_P;
    r_z = -1 + 2*(dot(grad_v_p,rp_s))/(V_do-V_up);
    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
    Ap=[Ap Fv_s*yf_r*0.5];
    Anb_vab =[Anb_vab Fv_s*V_up*(1-0.5*yf_r)];
    %Cara S , contribución de flujo masico de U
%}
%[Ap,Anb_vab]=esq_interp(Fv_s,v_S,v_P,grad_v_p,rp_s,Ap,Anb_vab,Fv_s);
[Ap,Anb_vab]=esq_interp_upwind(Fv_s,v_S,Ap,Anb_vab,Fv_s);


%    Ap=[Ap Fu_s*yf_r*0.5];
%    Anb_vab =[Anb_vab Fu_s*V_up*(1-0.5*yf_r)];
%[Ap,Anb_vab]=esq_interp(Fv_s,v_S,v_P,grad_v_p,rp_s,Ap,Anb_vab,Fu_s);
[Ap,Anb_vab]=esq_interp_upwind(Fv_s,v_S,Ap,Anb_vab,Fu_s);

%{
%else
    V_up=v_P;
    V_do=v_S;
    r_z = -1 + 2*(dot(grad_v_p,rp_e))/(V_do-V_up);
    yf_r= (r_z+ r_z^2)/(r_z^2 + 1);
    Ap=[Ap -Fv_s*(1-0.5*yf_r)];
    Anb_vab =[Anb_vab -Fv_s*V_do*yf_r*0.5];
    
    %Contribución de la cara S ,flujo másico en U
    Ap=[Ap -Fu_s*(1-0.5*yf_r)];
    Anb_vab =[Anb_vab -Fu_s*V_do*yf_r*0.5];

end
%}

%Dp/dy
dp_dy=(p_P-p_S)/(c_y(1)-c_y(2));
delta_V = 0.5*(dA_e +dA_w)*dx;%vol control



%----------------------------- ----------
%--------    Final Equation    ----------
% Final coefficient AP, left side of the equality;
Ap=Ap(~isnan(Ap));
Anb_vab=Anb_vab(~isnan(Anb_vab));
S_dc=S_dc(~isnan(S_dc));


A_P_i = sum(Ap); % Coefficient needed to calculate U_p and pressure correction
d_iJ=-delta_V/(A_P_i*(c_y(1)-c_y(2))); %coeficient pressure corrections  
V_i_j =(1/A_P_i)*(sum(Anb_vab) + sum(S_dc) - dp_dy*delta_V); %New U_P, main output

end
