function [P_i_j,Smb_t] = Pressure_corr(x,y,d_e,d_w,d_s,d_n,U_e,U_w,V_n,V_s,u_N_w,u_N_e,u_S_w,u_S_e,P_n,P_e,P_s,P_w)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Corrección de presiones 
%puntos de la celda 1,2,3,4 Pi,j , entrada x,y 
%x=[0 5 5 0] ejem.
%y=[5 6 2 0] ejem.
%coeficientes gradiente de presiones, entrada d_e,d_w,d_s,d_n
%d_e=1; 
%d_w=1;
%d_s=1;
%d_n=1;

%Velocidades sin corrección *, Entrada 
%U_e=Ui,j
%U_w=Ui,j-1
%V_n=Vi-1,j
%V_s=Vij
%Velocidades adicionales para interpolacion de flujo másico en X, caras N &
%S
%u_N_w =Ui-1,j-1
%u_N_e=2.1;%Ui-1,j
%u_S_w=1.09;%Ui+1,j-1
%u_S_e=2.31;%Ui+1,j

%presiones sin corrección *, entrada
%P_n = Pi-1,j
%P_e = Pi,j+1
%P_s = Pi+1,j
%P_w = Pi,j-1

%flujo masico en X ,caras norte y sur
U_n = 0.25*(u_N_w + u_N_e+ U_w+ U_e);%velocidad promediada en X cara N
U_s =  0.25*(u_S_w + u_S_e+ U_w+ U_e);%velocidad promediada en X cara S
Smb= -U_n*(y(2)-y(1)) +U_s*(y(3)-y(4)); %flujo masico en X ,caras norte y sur

%area lados laterales del V_C
dy_1=y(1)-y(4);
dy_2 =y(2)-y(3);
dx=x(3)-x(4);

AP_p = dy_1*d_w + dy_2*d_e + dx*d_s + dx*d_n;

%Balance de masa total 
Smb_t = -U_e*dy_2 + U_w*dy_1 + dx*(-V_n + V_s) +Smb;

%Presion corregida, salida  
P_i_j = (1/AP_p)*(P_w*dy_1*d_w + P_e*dy_2*d_e + P_s*dx*d_s + P_n*dx*-d_n + Smb_t);

end
