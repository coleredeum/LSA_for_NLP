%% Import Dataset

M = load('A_med.mat');
dict = load('dict_med.mat');
Q = load('Q_med.mat');

M = struct2array(M);
dict = struct2array(dict);
Q = struct2array(Q);
%

%% Plot of Big Document Matrix 

%figure
%mesh(M)

%

%% Plot of  Document Matrix 

%figure
%mesh(Q)

%

%% For 50 

[cos50, ind50] = overall(M, 50, Q);

%

%%

[cos500, ind500] = overall(M, 500, Q);

%

%%

[cos1033, ind1033] = overall(M, 1033, Q);

%

%% Plotting the singular values

[U,S,V] = svds(M, 1033);
loglog(diag(S));

%
%% Get approximation vector for one document at a time 
function dk = VectorApp(Vkt, i)
    Vkt = transpose(Vkt);
    dk = Vkt(:,i);
end
%
%% Get new query vector

function qk = NewQue(qt, Uk, Sk)
    Skinv = inv(Sk);
    qk = transpose(qt)*Uk*Skinv;
end
%
%% Compute Angle for one document 

function cosSim = CompAngle(a, b)
    cosSim = dot(a,b)/(norm(a)*norm(b));
end
%

%%
function [all_cosines, all_indices50] = overall(M, k, Q)
    tic 
    %[Uk, Sk, Vkt] = appR(M, 50);
    [U,S,V] = svds(M, k);

    [r c] = size(M);
    count = 0;
    all_indices50 = [];
    all_cosines = [];
    for j = 1:30
        count = count+1;
        q50 = NewQue(Q(:,j), U, S);
        results = zeros(1033,1);
        for i=1:c
            dk = VectorApp(V, i);
            cosSim = CompAngle(q50, dk);
            results(i, 1) = cosSim;
        end

        [resultsS, indices] = sort(results, 'descend');
        all_cosines = [all_cosines resultsS];
        all_indices50 = [all_indices50 indices];
    end
    toc 
end


%