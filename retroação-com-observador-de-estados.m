pkg load control
A= [ 0 1; 1 0];
B=[0;1];
C=[2 1];
gs=ss(A,B,C);
tf(gs);
pole(gs);
s1 = -4-i*5.46;
s2 = -4+i*5.46;
s3 = -2;
s4=-1;
K=acker(A,B,[s1 s2])
L = acker(A' ,C' ,[s3 s4])'
eig(A-B*K)
eig(A-L*C)
zero(ss(A, L,K))
Af =[A-B*K B*K; zeros(2) A-L*C];
Bf = [0;0;L];
Cf = [C 0 0 ];
ts= ss(Af,Bf,Cf);
tf(ts)
zero(ts)
pole(ts)
cs=ss(A-B*K-L*C,L,K)
gs=tf(gs)
cs=tf(cs)
tt=cs*gs/(1+cs*gs)
step(tt)