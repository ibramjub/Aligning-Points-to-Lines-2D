function [R,t,flag] = calcPose(p1,p2,p3,q1,q2,q3)
    % Find the 2D transformation (if exists) that aligns (p1,p2,p3) with 
    % (q1,q2,q3).
    
    n = 3;
    small_number = 0.001;
    flag = 0;
    
    P = [p1,p2,p3];
    Q = [q1,q2,q3];
    
    P_centered = P-repmat(mean(P,2),1,n);
    Q_centered = Q-repmat(mean(Q,2),1,n);
    
    PQ = P_centered*Q_centered';
    
    [U,D,V] = svd(PQ);
    
    d = det(V*U');
    R = V*[1,0;0,d]*U';
    R = V*U';
    
    diff = R*P_centered-Q_centered;
    err = sum(sqrt(sum(diff.^2,1)));
    
    t = mean(P,2)-mean(Q,2);
    t = -R*mean(P,2)+mean(Q,2);
    
    Q_test = R*P+t;
    diff = Q-Q_test;
    final_err = sum(sqrt(sum(diff.^2,1)));
    
    if (err < small_number)
        flag = 1;
    end 
end