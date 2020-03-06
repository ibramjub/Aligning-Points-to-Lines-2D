function wrapper()
    % A wrapper function for the main AlignPointsToLines.m function.
    % Generate a random set of 2D points and lines, and find their best
    % alignment.
    
    % Number of points
    N = 10;
    
    % Flag indicating if the points and alignments are generated with or
    % without error. If OptIs0 equals 1 then the points and lines can be 
    % perfectly aligned (with zero total error).
    OptIs0 = 1;
    
    % The set of points
    P_original = rand(2,N)*10;
    
    % The set of lines. The direction vectors of the lines are stored in L
    % and their distances from the origin are stored in B.
    L = rand(2,N); L = L./repmat(sqrt(sum(L.^2,1)),2,1);
    B = rand(1,N)*10; 
    
    if (OptIs0 == 1)
        for i=1:N
           vi_bot = null(L(:,i)');
           pi = vi_bot*B(i)+ randi(N,1)*L(:,i);
           P_original(:,i) = pi;
        end
    end
    
    
    % Add a random rotation R and translation t to the points
    t = rand(2,1)*10;
    v = rand(2,1);
    v = v/norm(v);
    v_bot = null(v');
    R = [v,v_bot];
    if (det(R) < 0)
        R(:,end) = -1*R(:,end);
    end

    % Apply rotation and translation to the points
    P = R*(P_original - repmat(mean(P_original,2),1,N)) - repmat(t,1,N);
    
    % Compute the set of candidate alignments (rotation matrices and
    % translation vectors). One of those alignments is the desired
    % approximate alignment.
    [R_candidates,t_candidates] = AlignPointsToLines(P,L,B);
    
    % Pick the alignment that minimizes the desired cost function:
    % ErrorMethod = 1: normal sum of euclidean distances^z
    % ErrorMethod = 2: M-estimator cost function with distance^z and threshold = th
    z = 2;
    ErrorMethod = 1;
    MEstimatorTh = 50;  % Only used when ErrorMethod = 2
    [R_approx,t_approx] = pickRT(P,L,B,R_candidates,t_candidates,ErrorMethod,z,MEstimatorTh);
    
    % Compute this alignment's error.
    Err = computeRTError_methods(P,L,B,R_approx,t_approx,z,ErrorMethod,MEstimatorTh);
end