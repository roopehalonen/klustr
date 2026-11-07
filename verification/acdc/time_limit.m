function varargout = time_limit(t, y)
% Wall-time limit of one steady state, used by ode15s as an event function (as on the server).
%   time_limit('set', seconds)   start the clock for a new condition
%   time_limit('expired')        true when the limit has been reached
% After the limit the event value changes sign at every call, so every ode15s run stops at its first step.
persistent clock limit n
if ischar(t)
    switch t
        case 'set'
            clock = tic; limit = y; n = 0;
        case 'expired'
            varargout = {toc(clock) > limit};
    end
    return
end
if toc(clock) <= limit
    value = 1;
else
    n = n + 1; value = (-1)^n;
end
varargout = {value, 1, 0};
end
