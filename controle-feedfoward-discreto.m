pkg load control
%matriz identidade
I = [1 0;0 1];
%matrizes A, B e C de G1
A = [0,1;0,-1];
B = [0;1];
C = [1 0];
%Primeiro tempo de amostragem: 200ms
T = 0.2;
%polos em tempo discreto
z1 = 0.207 + i*0.3988;
z2 = 0.207 - i*0.3988;
%matrizes A e B em tempo discreto
Ab = [1,0;0,1] + A*T + A^2*(T^2/2) + A^3*(T^3/6)
Bb = ([1,0;0,1]*T + A*(T^2/2) + A^2*(T^3/6) + A^3*(T^4/24))*B
%calculo de K1 e K2 em tempo discreto
Kb = acker(Ab,Bb,[z1 z2])
%Feedfoward em tempo discreyo
F1 = -inv(C*inv(Ab - Bb*Kb -I)*Bb)
%matrizes Ab e Bb em malha fechada
Abf = Ab - Bb*Kb
Bbf = Bb*F1
%função em espaço de estados em malha fechada
T1 = ss(Abf, Bbf, C, 0, T)
%função de transferncia de T1
tf(T1)
figure
%resposta ao degrau
step(T1)

%o mesmo procedicmento se repete a seguir, porém com tempo de amostragem 
%igual a 50ms
T = 0.05;
z1 = 0.7884 + i*0.2207;
z2 = 0.7884 - i*0.2207;
Ab = [1,0;0,1] + A*T + A^2*(T^2/2) + A^3*(T^3/6)
Bb = ([1,0;0,1]*T + A*(T^2/2) + A^2*(T^3/6) + A^3*(T^4/24))*B
Kb = acker(Ab,Bb,[z1 z2])
F2 = -inv(C*inv(Ab - Bb*Kb -I)*Bb)
Abf = Ab - Bb*Kb
Bbf = Bb*F2
T2 = ss(Abf, Bbf, C,0, T)
tf(T2)
figure
step(T2)







