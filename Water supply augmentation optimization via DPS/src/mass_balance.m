function [s1,r1] = mass_balance( s, u, q, storage )

% Output:
%      s1 - final storage. 
%      r1 - release over the 24 hours. 
%
% Input:
%       s - initial storage. 
%       u - release decision 
%       q - inflow 


HH = 2; % integration substep
delta = 1/HH;
s_ = nan(HH+1,1);
r_ = nan(HH+1,1);

evap =  0 ;
s_(1) = s;

for i=1:HH
    r_(i+1) = min( s_(i) , max( s_(i) - storage , u ) );
    s_(i+1) = s_(i) + delta*( q - r_(i+1) ) - evap ;
end

s1 = s_(HH);
r1 = mean(r_(2:end));
