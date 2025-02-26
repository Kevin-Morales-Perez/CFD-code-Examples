%Graficos
clear all
close all
clc
tic

load("mesh_bumpchannel2.mat");
load("bumpinchannel_fields2.mat") 

[n_y,n_x] =size(X);%numero de celdas en X y Y
n_x=n_x-1;
n_y=n_y-1;

dx=X(1,2)-X(1,1);

X1=zeros(n_y,n_x);
Y1=zeros(n_y,n_x);

for i=1:n_y
        for j=1:n_x
            x=[X(n_y+2-i,j),X(n_y+1-i,j),X(n_y+1-i,j+1),X(n_y+2-i,j+1)];
            y=[Y(n_y+2-i,j),Y(n_y+1-i,j),Y(n_y+1-i,j+1),Y(n_y+2-i,j+1)];
            m1=(y(3)-y(1))/dx;
            m2=(y(4)-y(2))/dx;
            c2=y(2)-y(1);
            X1(n_y+1-i,j)=x(1) + c2/(1- m2/m1);
            Y1(n_y+1-i,j)=y(1)+ c2/(m1-m2);

        end
end

U_plot=zeros(n_y,n_x);
V_plot=zeros(n_y,n_x);
P_plot=zeros(n_y,n_x);
B_plot=zeros(n_y,n_x);

for i =1:n_y
    for j=1:n_x
        U_plot(i,j)=0.5*(U_vel(i+1,j)+U_vel(i+1,j+1));
        V_plot(i,j)=0.5*(V_vel(i,j+1)+V_vel(i+1,j+1));
        P_plot(i,j)=0.25*(P_press(i,j+1)+P_press(i+1,j+2)+P_press(i+2,j+1)+P_press(i+1,j));
        B_plot(i,j)=0.25*(b_p(i,j+1)+b_p(i+1,j+2)+b_p(i+2,j+1)+b_p(i+1,j));
    end
end

V_mag = sqrt(U_plot.^2 + V_plot.^2);%MAgnitud de la velocidad
toc

figure(21);
contourf(X1,Y1,U_plot, 21, 'LineStyle', 'none')
title('Velocidad en X')
xlabel('x')
ylabel('y')
colorbar
colormap('jet')

figure(22);
contourf(X1,Y1,V_plot, 21, 'LineStyle', 'none')
title('Velocidad en Y')
xlabel('x')
ylabel('y')
colorbar
colormap('jet')

figure(23);
title("Campo de Velocidades")
hold on
quiver(X1, Y1, U_plot, V_plot, 5, 'k')

figure(24)
contourf(X1,Y1,P_plot)
title("Distribucion de presiones")
xlabel('x','FontSize',16)
ylabel('y','FontSize',16)
colorbar
colormap('jet')

figure(25);
contourf(X1,Y1,V_mag, 21, 'LineStyle', 'none')
title("Magnitud de la velocidad")
xlabel('x')
ylabel('y')
colorbar
colormap('jet')

figure(26);
contourf(X1,Y1,B_plot, 21, 'LineStyle', 'none')
title("Ecuacion de continuidad")
xlabel('x')
ylabel('y')
colorbar
colormap('jet')



