close;clear
T=1:0.01:3.5;
a = length(T);
ET = zeros(a,4);
MT = zeros(a,4);
CT = zeros(a,4);
Xi = zeros(a,4);
N=100;% the code generate lattice of size 2N*2N
cr_length =zeros(a,4);
w = zeros(4,N);
tic
% we use parfor instead of for loop to utlize the parallel computing
% feature and use all the cores in the CPU 
parfor i = 1:4
    % initial condtion here is not completly random 75% of the spins are
    % either up or down to avoid metstable states
    s = ((i<3)*(rand(2*N)<0.25)+ (i>2)*(rand(2*N)>0.25) )* -2 + 1;
    for j = 1:a
        [s,E,M,C,x]=ising_sim_B(s,N,T(j),0);
        ET(j,i) = E;% energy
        MT(j,i) = M;% magnatization
        CT(j,i) = C;% heat capacity
        Xi(j,i) = x;
        % calculate corrlation length from correlation function
        cor = corrlation(s,N);
        m = max(max(cor));
        r = find(diff(cor <= m*0.368));%correlation length
        if isempty(r)||(m <(10^-3)) r = 0; end %for small values result 
        cr_length(j,i)= r(1) ;                  % not accurate
        [T(j) i] % printed  just to show code is running
    end
end
toc
subplot(221)
plot(T,mean(ET,2))
subplot(222)
plot(T,mean(CT,2))
subplot(223)
plot(T,MT,'o')
subplot(224)
plot(T,mean(Xi,2),'o')
figure
plot(T,mean(cr_length,2))
%% ising model with magnatic field
N =100;
B = -1.5:0.02:1.5;
T = 1.5:0.5:4.5;
a = length(B);
b = length(T);
MT = zeros(a,b);
parfor j =1:a
    for i =1:b
        s = (rand(2*N)<0.5)*-2 +1;
        [~,~,M,~,~] = ising_sim_B(s,N,T(i),B(j));
        MT(j,i) = M;
        [T(i) B(j)]
    end
end
figure
plot(B,MT)
legend("T =1.5 ","T =2 ","T =2.5 ","T =3 ","T =3.5 ","T =4 ")
%% show the lattice
T=[1 1.5 2 2.27 2.5 3 4];
a = length(T);
N=100;
for i = 1:2
    % initial condtion here is not completly random 75% of the spins are
    % either up or down to avoid metstable states
    s = ((i<2)*(rand(2*N)<0.25)+ (i>1)*(rand(2*N)>0.25) )* -2 + 1;
    for j = 1:a
        [s,~,~,~,~]=ising_sim_B(s,N,T(j),0);
        figure; imagesc(s)
        title("temperature is  " + T(j))
        [T(j) i] % printed  just to show code is running
    end
end
%% block spin transformation
N = 243;
s =(rand(2*N)<0.5)* -2 + 1;
[s,~,~,~,~]=ising_sim_B(s,N,2.5,0);
s1 = block_s_tr(s,2*N);
s2 = block_s_tr(s1,2*N/3);
s3 = block_s_tr(s2,2*N/9);
subplot(221);title("first")
imshow(s,'InitialMagnification',3000)
subplot(222);title("second")
imshow(s1,'InitialMagnification',3000)
subplot(223);title("third")
imshow(s2,'InitialMagnification',3000)
subplot(224);title("third")
imshow(s3,'InitialMagnification',3000)
%%  ising model simulation function 
function [s,E, M, C,x] = ising_sim_B(system,N,T,B)
    l = 2000;% number of iteration
    s =system;
    M = zeros(l,1);
    E = zeros(l,1);
    chosed(:,:,1) = repmat([1 0;0 0],N);
    chosed(:,:,2) = repmat([0 1;0 0],N);
    chosed(:,:,3) = repmat([0 0;1 0],N);
    chosed(:,:,4) = repmat([0 0;0 1],N);
    for iter = 1:l
        DeltaE = 2 * (s.* (circshift(s,[ 0 1])+ circshift(s,[ 0 -1])+...
                circshift(s, [ 1 0])+ circshift(s, [-1 0])+B)) ;
        i = randi([1,4]);
        transitions = ((rand(2*N) < exp(-DeltaE/T)) & chosed(:,:,i))* -2 +1;
        s = s.* transitions;
        M(iter) = sum(sum(s));
        E(iter) = -sum(sum(DeltaE))/2;
    end
    C = var(E(200:end)).*(2*T*N)^-2; %we do not count the first 200
    E = mean(E(200:end))/(4*N^2);  % iteration because the system might not
    x = var(M(200:end))./(4*T*N^2);% be able to reach equiliprium
    M = mean(M(200:end))/(4*N^2);
end
%% block spin transformation
function s = block_s_tr(system,N)
    s =zeros(N/3);
    for i = 1:N/3
        for j = 1:N/3
            s(i,j) = (sum(sum(system(3*i-2:3*i,3*j-2:3*j)))<0)*-2 + 1;
        end
    end
end

%% correlation
function cr = corrlation(s,N)
    cr = zeros(1,N);
    for i = 1:N
        A = 0.25*s.* (circshift(s,[ 0 i])+ circshift(s,[ 0 -i])+...
                circshift(s, [ i 0])+ circshift(s, [-i 0]));
            cr(i) = mean(mean(A))-mean(mean(s)).^2;
    end
end