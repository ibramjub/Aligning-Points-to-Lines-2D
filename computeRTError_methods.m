function error = computeRTError_methods(P,L,B,R,t,z,method,th)
    % Given a candidate alignment (R,t), compute its error with the
    % specified method:
    % method = 1: normal sum of euclidean distances^z
    % method = 2: M-estimator cost function with distance^z and threshold = th
    % The M-estimator cost for each point is min{computeRTError(...,z),threshold}
    % Li contains the direction vector of li
    
    N = length(P);
    error = 0;
    
    for i=1:N
        bi = B(i);
        vi_bot = null(L(:,i)');
        pi = P(:,i);
        pi_aligned = R*pi+t;
        
        erri = norm(vi_bot'*pi_aligned-bi)^z;
        if (method == 2)
            erri = min(erri,th);
        end
        
        error = error + erri;
    end
end