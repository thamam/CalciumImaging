function evt = ogatapoisson(rate,rmax,tmax,pmax)
% OGATAPOISSON   generate Poisson process by Ogata's thinning algorithm
%    T = OGATAPOISSON(RATE,RMAX,TMAX,PMAX) generates a set of event times
%    T drawn from an inhomogeneous Poisson proces with rate RATE. RATE is a
%    function handle that returns the rate of the process at a specified
%    time RATE(t). RMAX should be an upper bound such that RATE(t) <= RMAX.
%    The generation quits when either the time reaches TMAX or the number
%    of points generated reaches PMAX.
 
if isinf(tmax) && isinf(pmax)
    error('No stopping criteria - both TMAX and PMAX are infinite');
end
if isfinite(pmax)
    evt = zeros(pmax,1); % can pre-allocate
else
    evt = zeros(100,1); % may need more, but start with something
end
 
t = 0;
i = 0;
while t < tmax && i < pmax
    t = t-log(rand)/rmax; % advance to next candidate time
    ratefrac = rate(t)/rmax;
    if ratefrac > 1 % cheap enough to check
        fprintf('T is %f, Rate is %f, max rate is %f, frac is %f\n', t, rate(t), rmax, ratefrac)
        error('RATE evaluated to a value greater than RMAX')
    elseif ratefrac < 0
        error('RATE evaluated to a negative value')
    elseif rand < ratefrac && t < tmax % do we keep it?
        i = i+1;
        evt(i) = t;
    end
end
evt = evt(1:i); % trim excess entries, if necessary
end