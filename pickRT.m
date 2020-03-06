function [R,t] = pickRT(P,L,B,R_all,t_all,method,z,th)
    % Given a set of many alignments, pick the best alignment for the
    % chosen cost function:
    % method = 1: normal cost function with norm^z
    % method = 2: M-estimator cost function with norm^z and threshold = th
    % The M-estimator cost for each point is min{computeRTError(...,z),threshold}
    
    if (method == 2)
       assert(exist('th','var') == 1);
    elseif (~exist('th','var'))
       th = 1; 
    end
        
    Err = Inf;
    for i=1:2:length(R_all)
       curr_R = R_all(i:i+1,:); 
       curr_t = t_all(i:i+1); 
       curr_Err = computeRTError_methods(P,L,B,curr_R,curr_t,z,method,th);
       if (curr_Err < Err)
           Err = curr_Err;
           Idx = i;
       end
    end
    R = R_all(Idx:Idx+1,:);
    t = t_all(Idx:Idx+1);
end