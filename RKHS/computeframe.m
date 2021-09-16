function [tminInd, mintime] = computeframe(supeps, delta, k )
%computeframe
% Find the minimal distance where k(|tmin|)<supeps
% % cind = mintime in delta units ( time indexes units)
% %
%%

N = 100;
for i=1:10
    timebrackets = (0:delta:N*delta); %prepare delta-spaced grid
    keval = k(0,timebrackets);        % evaluate k(0,t) for all t on the grid
    tminInd = find(keval<supeps,1,'first'); % find first i where  k(0,t_i)<supeps
    if ~isempty(tminInd) % if found such t_i terminate and return tmin
        mintime = tminInd*delta;       
        return;
    else %if not, extend search range
        N=N*2;
    end
end

 if isempty(tminInd)
     error('can not find k<eps in range\n');
 end


end

