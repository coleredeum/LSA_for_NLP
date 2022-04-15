%% Import Dataset

M = load('A_med.mat');
dict = load('dict_med.mat');
Q = load('Q_med.mat');

M = struct2array(M);
dict = struct2array(dict);
Q = struct2array(Q);
%

%% Plot of Big Document Matrix 

figure
mesh(M)

%

%% Plot of Dict Document Matrix 

figure
mesh(dict)

%

%% Plot of  Document Matrix 

figure
mesh(Q)

%

%% Cutted U, S, V 

function [Uk, Sk, Vkt] = appR(A, k)
    [U,S,V] = svd(A);
    Uk = U(:, 1:k);
    Sk = S(1:k, 1:k);
    Vkt = V(:, 1:k)';
end
%
%% Get approximation vector for one document at a time 
function dk = VectorApp(Vkt, i)
    dk = Vkt(:,i);
end
%
%% Get new query vector

function qk = NewQue(qt, Uk, Sk)
    Skinv = inv(Sk);
    qk = qk*Uk*Skinv;
end
%
%% Compute Angle for one document 

function cosSIm = CompAngle(a, b)
    cosSim = dot(a,b)/(norm(a)*norm(b));
end
%