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

%%

Allcos = zeros(30,1033);
tic
for k=1:1:30
    for l = 1:1:1033
        Allcos(k,l) = CompAngle(Q(:,k),M(:,l));
    end
end
toc
Allcos;

%%

disp(size(Allcos(1,:)'));
Allcos(Allcos == 0) = NaN;
save('Allcos', 'Allcos');
nnz(~isnan(Allcos(1,:)))
nnz(~isnan(Allcos(2,:)))
nnz(~isnan(Allcos(3,:)))
nnz(~isnan(Allcos(4,:)))
nnz(~isnan(Allcos(5,:)))
Relations = [];
for i=1:1033
    Relations = [Relations nnz(~isnan(Allcos(1:5,i)))];
end

figure()
hold on
plot(Allcos(1,:)', '-');
plot(Allcos(2,:), '-')
plot(Allcos(3,:), '-')
plot(Allcos(4,:), '-')
plot(Allcos(5,:), '-')
xlabel('Documents')
ylabel('Cosine Values')
legend('Query 1','Query 2','Query 3','Query 4','Query 5')
hold off

%

%%

a1 = sort(Allcos(1,:), 'descend');
a2 = sort(Allcos(2,:), 'descend');
a3 = sort(Allcos(3,:), 'descend');
a4 = sort(Allcos(4,:), 'descend');
a5 = sort(Allcos(5,:), 'descend');
figure()
hold on
plot(Allcos(1,:), '-');
plot(Allcos(2,:), '-');
plot(Allcos(3,:), '-');
plot(Allcos(4,:), '-');
plot(Allcos(5,:), '-');
%xlim([0 400])
xlabel('Documents')
ylabel('Cosine Values')
legend('Query 1','Query 2','Query 3','Query 4','Query 5')
hold off


%
%% Plotting the relation documents to queries 1-5:

figure()
hold on
plot(Relations, 'o');
xlabel('Documents')
ylabel('Corresponding queries')
hold off

%

%% For 50 

[cos50, ind50] = overall(M, 50, Q);
save('cos50', 'cos50');
save('ind50', 'ind50');

%
%%
cos50 = transpose(cos50);
%
%%

b1 = sort(cos50(1,:), 'descend');
b2 = sort(cos50(2,:), 'descend');
b3 = sort(cos50(3,:), 'descend');
b4 = sort(cos50(4,:), 'descend');
b5 = sort(cos50(5,:), 'descend');
figure()
hold on
plot(cos50(:,1), '-');
plot(cos50(:,2), '-');
plot(cos50(:,3), '-');
plot(cos50(:,4), '-');
plot(cos50(:,5), '-');
xlabel('Documents')
ylabel('Cosine Values')
legend('Query 1','Query 2','Query 3','Query 4','Query 5')
hold off


%
%%

[cos500, ind500] = overall(M, 500, Q);

save('cos500', 'cos500');
save('ind500', 'ind500');
%
%% Ordering of cosines

a1 = sort(cos500(:,1), 'descend');
a2 = sort(cos500(:,2), 'descend');
a3 = sort(cos500(:,3), 'descend');
a4 = sort(cos500(:,4), 'descend');
a5 = sort(cos500(:,5), 'descend');
figure()
hold on
plot(cos500(:,1), '-');
plot(cos500(:,2), '-');
plot(cos500(:,3), '-');
plot(cos500(:,4), '-');
plot(os500(:,5), '-');
xlabel('Documents')
ylabel('Cosine Values')
legend('Query 1','Query 2','Query 3','Query 4','Query 5')
hold off

%
%%

[cos1033, ind1033] = overall(M, 1033, Q);
save('cos1033', 'cos1033');
save('ind1033', 'ind1033');

%

%% Plotting the singular values
[U1,S1,V1] = svds(M, 50);
[U,S,V] = svds(M, 500);
figure()
hold on
loglog(diag(S1), 'x');
loglog(diag(S));
xlabel('i-th singular value');
ylabel('Singular value');
legend('50 largest SV','500 largest SV')
grid on  
hold off
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
    %[Uk, Sk, Vkt] = appR(M, 50);
    [U,S,V] = svds(M, k);

    [r c] = size(M);
    count = 0;
    all_indices50 = [];
    all_cosines = [];
    new_Q = NewQue(Q, U, S);
    tic 
    for j = 1:30
        count = count+1;
        %q50 = NewQue(Q(:,j), U, S);
        %disp(size(new_Q));
        q50 = new_Q(j, :);
        results = zeros(1033,1);
        for i=1:c
            dk = VectorApp(V, i);
            cosSim = CompAngle(q50, dk);
            results(i, 1) = cosSim;
        end
        all_cosines = [all_cosines results];
    end
    toc 
    [p, k]= size(all_cosines);
    all_indices = [];
    all_cos_final = [];
    for i=1:k
        [a1, b1] = sort(all_cosines, 'descend');
        all_cos_final = [all_cos_final a1];
        all_indices = [all_indices b1];
    end 
end


%