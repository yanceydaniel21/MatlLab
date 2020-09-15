clear all
close all
clc

T1 = 300; %Kelvin
M1 = 0.6; %M1
gamma1 = 1.4; %for air 
D1 = 6/39.37; %inches to meters
Q1 = 2088*1055.06; %Btu/s to Watts

L_crit1 = 4.25; %used from Fanno Flow

x1 = 0; %entrance
n1 = 10; %steps
h1 = L_crit1/n1;
cp1 = 1005; %J/kg*K
m_dot1 = 15.43; %kg/s
T_o1 =T1*(1+(0.5*(gamma1-1))*M1^2);

x_1 = 0:h1:L_crit1;
M_1 = zeros(1,n1);


y=@(x_1,M_1)   M_1*[1+0.5*(gamma1-1)*M_1^2]/(1-M_1^2)*((1+gamma1*M_1^2)/(2*(T_o1+((Q1/(pi()*D1*L_crit1))*x_1)/(L_crit1*m_dot1*cp1)) )*((Q1/(pi()*D1*L_crit1)))/(L_crit1*m_dot1*cp1));

M_1(1) = M1;


for i = 1:n1
    K1 = y(x_1(i),M_1(i));
    K2 = y(x_1(i) + h1/2 , M_1(i) + K1*h1/2);
    K3 = y(x_1(i) + h1/2 , M_1(i) + K2*h1/2);
    K4 = y(x_1(i) + h1 , M_1(i) + K3*h1);
    M_1(i+1) = M_1(i) + (K1 + 2*(K2 + K3) + K4)*h1/6;
end

plot (x_1/L_crit1 , M_1, '-o')
grid on
hold on

T2 = 300; %Kelvin
M2 = 1.2; %M1
gamma2 = 1.4; %for air 
D2 = 6/39.37; %inches to meters
Q2 = 2088*1055.06; %Btu/s to Watts

L_crit2 = 35; %used from Fanno Flow

x2 = 0; %entrance
n2 = 10; %steps
h2 = L_crit2/n2;
cp2 = 1005; %J/kg*K
m_dot2 = 15.43; %kg/s
T_o2 =T2*(1+(0.5*(gamma2-1))*M2^2);

x_2 = 0:h2:L_crit2;
M_2 = zeros(1,n2);


y=@(x_2,M_2)   M_2*[1+0.5*(gamma2-1)*M_2^2]/(1-M_2^2)*((1+gamma2*M_2^2)/(2*(T_o2+((Q2/(pi()*D2*L_crit2))*x_2)/(L_crit2*m_dot2*cp2)) )*((Q2/(pi()*D2*L_crit2)))/(L_crit2*m_dot2*cp2));

M_2(1) = M2;


for i = 1:n2
    K1 = y(x_2(i),M_2(i));
    K2 = y(x_2(i) + h2/2 , M_2(i) + K1*h2/2);
    K3 = y(x_2(i) + h2/2 , M_2(i) + K2*h2/2);
    K4 = y(x_2(i) + h2 , M_2(i) + K3*h2);
    M_2(i+1) = M_2(i) + (K1 + 2*(K2 + K3) + K4)*h2/6;
end

plot (x_2/L_crit2 , M_2, '-o')
hold on

T3 = 300; %Kelvin
M3 = 1.8; %M1
gamma3 = 1.4; %for air 
D3 = 6/39.37; %inches to meters
Q3 = 2088*1055.06; %Btu/s to Watts

L_crit3 = 3.1; %used from Fanno Flow

x3 = 0; %entrance
n3 = 10; %steps
h3 = L_crit3/n3;
cp3 = 1005; %J/kg*K
m_dot3 = 15.43; %kg/s
T_o3 =T3*(1+(0.5*(gamma3-1))*M3^2);

x_3 = 0:h3:L_crit3;
M_3 = zeros(1,n3);


y=@(x_3,M_3)   M_3*[1+0.5*(gamma3-1)*M_3^2]/(1-M_3^2)*((1+gamma3*M_3^2)/(2*(T_o3+((Q3/(pi()*D3*L_crit3))*x_3)/(L_crit3*m_dot3*cp3)) )*((Q3/(pi()*D3*L_crit3)))/(L_crit3*m_dot3*cp3));

M_3(1) = M3;


for i = 1:n3
    K1 = y(x_3(i),M_3(i));
    K2 = y(x_3(i) + h3/2 , M_3(i) + K1*h3/2);
    K3 = y(x_3(i) + h3/2 , M_3(i) + K2*h3/2);
    K4 = y(x_3(i) + h3 , M_3(i) + K3*h3);
    M_3(i+1) = M_3(i) + (K1 + 2*(K2 + K3) + K4)*h3/6;
end

plot (x_3/L_crit3 , M_3, '-o')
hold on

xlim([0 1]);
ylim([0 2]);

legend('M=0.6' , 'M=1.2', 'M=1.8')

title ('Rayleigh Flow')
xlabel ('x/L')
ylabel ('Mach, M')