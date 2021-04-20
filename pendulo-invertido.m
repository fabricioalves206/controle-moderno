pkg load control; 
clear all; 
close all; 
clc; 

global M m l g I f; 

M = 1.00; 
m = 0.30; 
l = 0.25; 
g = 9.81; 

I = m*(2*l)^2/12; 

as = [ 0 0 1 0; 0 0 0 1; 0 g*(m*l)^2/(I*(M+m)+M*m*l^2) 0 0; 0 g*m*l*(M+m)/(I*(M+m)+M*m*l^2) 0 0 ]; 
bs = [ 0; 0; (I+m*l^2)/(I*(M+m)+M*m*l^2); m*l/(I*(M+m)+M*m*l^2) ]; 
cs = [ 1 0 0 0 ]; 

gs = ss(as, bs, cs); 
tf(gs) 

sist = -5.96572;
s1 = -2+i*2.7288; 
s2 = -2-i*2.7288; 
s3 = -5.4249; 
s4 = -5.4249; 

wsist = abs(sist)
  
wa=2*wsist

fa= wa/(2*pi)
ta= 1/fa

ks = acker(as, bs, [s1 s2 s3 s4]) 
fs = - inv(cs*inv(as-bs*ks)*bs) 

% Condicao inicial 
rs = 0; 
x0 = pi - 20/180*pi; 

% Simulacao nao-linear com a retroacao de estados 
ti = 0.001; 
xs = [ 0 x0 0 0; 0 x0 0 0 ]; 

for i = 1:12001 
  if i == 6001 
    rs = 1; 
  end 
  f = fs*rs - ks(1)*xs(2,1) - ks(2)*( xs(2,2) - pi ) - ks(3)*xs(2,3) - ks(4)*xs(2,4); 
  xs = lsode("naolinear",xs(2,:), [0 ti]); 
  ts(i) = (i-1)*ti; 
  us(i) = f; 
  y1(i) = xs(2,1); 
  y2(i) = xs(2,2); 
end 

polos = pole(gs)
zeros = zero(gs)


ta=0.1
z1 = exp(s1*ta)
z2 = exp(s2*ta)
z3 = exp(s3*ta)
z4 = exp(s4*ta)

gz = c2d(gs, ta)

[az, bz, cz, dz, tz] = ssdata(gz);

mc = ctrb(az, bz)
det(mc)

kz = acker(az, bz, [z1 z2 z3 z4])
fz = -inv(cz*inv(az-bz*kz-eye(4))*bz)

% Condicao inicial 
rz = 0; 
x0 = pi - 20/180*pi; 

% Simulacao nao-linear com a retroacao de estados 
xz = [ 0 x0 0 0; 0 x0 0 0 ]; 

for i = 1:12001 
  if i == 6001 
    rz = 1; 
  end 
  if(rem(((i-1)*ti),ta) == 0)
    f = fz*rz - kz(1)*xz(2,1) - kz(2)*(xz(2,2)- pi) - kz(3)*xz(2,3) - kz(4)*xz(2,4); 
  end
  xz = lsode("naolinear", xz(2,:), [0 ti]); 
  tz(i) = (i-1)*ti; 
  uz(i) = f; 
  y1(i) = xz(2,1); 
  y2(i) = xz(2,2); 
end 

figure; 
plot(tz, uz); 
hold on
plot(ts, us);

grid; 
xlabel('Tempo [s]'); 
ylabel('Forca [N]'); 

figure; 
plot(tz, 100*y1); 
grid; 
xlabel('Tempo [s]'); 
ylabel('Posicao do carro [cm]'); 

figure; 
plot(tz, 180/pi*y2); 
grid; 
xlabel('Tempo [s]'); 
ylabel('Angulo do pendulo [graus]');