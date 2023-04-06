function y = rbf(param, N, M, K, inputs)
% N: number of Basis
% M: number of inputs
% k: number of outputs

%% assign paramters
const = param(1:K); % constant for output
c = nan; %center of gaussians
b = nan; %radius of gaussians
w = nan; %weights

count = K+1;
for i=1:N
    for j=1:M
        c = [c, param(count) ];
        b = [b, param(count + 1) ];
        count = count+2;
    end
end

for k = 1:K
    w = [w, param(count) ];
    count = count+1;
end

c = c(2:end);
b = b(2:end);
w = w(2:end);

%% calculate gaussian output
count_c = 1;
count_b = 1;
phi = nan(N,1);
for j = 1:N
    bf = 0;

    for i = 1:M
        num = ( inputs(i) - c(count_c) )^2;
        den = b(count_b)^2;
        count_c = count_c + 1;
        count_b = count_b + 1;

        if den < 10^(-6)
            den = 10^(-6);
        end

        bf = bf + num / den ;
    end
    phi(j) = exp(-bf);
end

%% calculate RBF output
y = nan(1, K);
for k = 1:K
    o = const(k);
    for i = 1:N
        o = o + w(k)*phi(i);
    end
    if o < 0
        o  = 0;
    elseif o > 1
        o = 1;
    end 
    y(k) = o;
end

