function [x,g,E,r,epsilon] = ...
    dcHMCRunOnce(f,nTau, x,g,E,r,epsilon) 
r_slowness = 0.99;
target_r   = 0.651;
epsilon_inc= 1.1;

if isempty(g)
    [E,g] = f(x);                   % set initial gradient and Energy
end
    
p = randn(size(x));                 % initial momentum is Normal(0,1)
H = p * p' / 2 + E;                 % evaluate H(x,p)
    
xnew = x;
gnew = g;
    
for tau = 1:nTau                    % make nTau 'leapfrog' steps
    p = p - epsilon * gnew / 2;     % make half-step in p
    xnew = xnew + epsilon * p;      % make step in x 
    [~,gnew] = f(xnew);             % find new gradient    
    p = p - epsilon * gnew / 2;     % make half-step in p
end
    
[Enew,~] = f(xnew);                 % find new value of H
Hnew = p*p'/2 + Enew;
dH = Hnew - H;                      % decide wheter to accept
    
if dH<0 || rand < exp(-dH)
    accept = true;
    x = xnew;
    g = gnew;
    E = Enew;
else
    accept = false;
end

% adjust step size epsilon to acceptance rate
r = r * r_slowness + accept * (1-r_slowness);
if r>target_r
    epsilon = epsilon * epsilon_inc;
else
    epsilon = epsilon / epsilon_inc;
end
epsilon = epsilon * (1+(0.4*rand-0.2)); % add some random variation
end
