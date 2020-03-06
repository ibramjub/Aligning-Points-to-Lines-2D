function [R_approx,t_approx] = ComputeAlignment(p1,p2,p3,v1,v2,v3,b1,b2,b3)
    % Given 3 points and 3 lines (L(:,i) is the direction vector of li, 
    % B(i) is its distance from the origin), calculate a set of O(1)
    % alignments as in ComputeAlignment(p1,p2,p3,l1,l2,l3)
    % ** dist(pi,li) = ||li_bot'pi - bi||
    % 1) calculate (R,t) instantly, not afterwards
    % 2) test code compared to original one
    % 3) change def of p_li
    
    % This version: Fix to the sign of bi + perp vec
    
    small_number = 0.01;
    plotFlag = 0;
    
    P_trig = [];
    Q_trig = [];
    Z_trig = [];
    R_approx = [];
    t_approx = [];
    
    % Calculate all pairwise auclidean distances
    r12 = norm(p1-p2);
    r13 = norm(p1-p3);
    r23 = norm(p2-p3);
    
    % Perp vectors (fixed by the sign of bi
    v1_bot = sign(b1)*null(v1');
    v2_bot = sign(b2)*null(v2');
    v3_bot = sign(b3)*null(v3');
    b1 = abs(b1);
    b2 = abs(b2);
    b3 = abs(b3);
    
    % The slopes of l1 and l2
    slope_1 = v1(2)/v1(1);
    slope_2 = v2(2)/v2(1);

    % p1 \in l1 and p2 \in l2
%     p_l1 = v1_bot*b1;
%     p_l2 = v2_bot*b2;
%     p_l3 = v3_bot*b3;
    p_l1 = v1_bot*b1+v1;
    p_l2 = v2_bot*b2+v2;
    p_l3 = v3_bot*b3+v3;

    % Intersection points of l1 and l2 with the y-axis
    intersect_1y = p_l1(2)-slope_1*p_l1(1);
    intersect_2y = p_l2(2)-slope_2*p_l2(1);

    if (abs(abs(dot(v1,v2))-1) > small_number)
        % Intersection of l1 and l2
        intersect_12 = lineLineIntersection(v1,p_l1,v2,p_l2);

        % Calculate the matrix Zi for point i as in the Zx theorem
        [Z1,Z2,Pmat,Qmat] = calcZ(v1,v2,r12,r13,r23);

        % Pick the right matrix (Z1 or Z2)
        Zmat = checkZ([p1,p2,p3],Pmat,Qmat,Z1,Z2);

        % Computing the vector ai
        a = v3_bot'*Zmat;

        % Approximating |ax-b3| when assuming that l1 and l2
        % intersect at the ORIGIN! so b3 (the dist from l3 to the original origin)
        % should be changed accordingly
        p_l3_new = p_l3-intersect_12;
        b3_new = p_l3_new'*v3_bot; %%%%%%%% MAYBE NEGAATIVE, so change the direction vector a1 / a2 accordingly
        sign_b3_new = sign(b3_new);
        b3_new = abs(b3_new);
        if (abs(b3_new) < small_number)
           sign_b3_new = 1; 
        end

        % Find unit vectors that minimize |a1'x-b3_new| and |a2'x-b3_new| (if
        % l1 and l2 aren't parallel)
        x = solve_ax_b(sign_b3_new*a',b3_new);
    end
                  
    % Case(I)
    % If l1 and l2 aren't parallel, do the usual ellipse case
    if (abs(abs(dot(v1,v2))-1) > small_number)
        for currVec=1:size(x,2)
            p_approx = Pmat*x(:,currVec)+intersect_12;
            q_approx = Qmat*x(:,currVec)+intersect_12;
            z_approx = Zmat*x(:,currVec)+intersect_12;

            [R,t,flag] = calcPose(p1,p2,p3,p_approx,q_approx,z_approx);
            if (flag)
                R_approx = [R_approx;R];
                t_approx = [t_approx;t];
            end

            % Assert that: p1 \in l1, p2 \in l2
            assert(norm(v1_bot'*(p_approx)-b1) < small_number);
            assert(norm(v2_bot'*(q_approx)-b2) < small_number);
        end      
    
        % Case(II)
        % If l1 and l2 aren't parallel, compute all possible triplets of
        % points such that (and p1 \in l2, p2 \in l2 and (p2-p1) is perp to l2)
        [P_options,Q_options,Z_options] = Trig(r12,r13,r23,v1,v2,b1,b2);

        for o=1:length(Z_options)
            % Find intersection of l3 and third point when moving on a line parallel
            % to l1 (and p1 \in l2, p2 \in l2 and (p2-p1) is perp to l2)
           curr_intersect = lineLineIntersection(v3,p_l3,v1,Z_options(:,o));

           % If z cannot intersect it's line then dont move at all
           isBadNumber = sum(isnan(curr_intersect))+sum(isnan(curr_intersect));
           if (isBadNumber == 0)
               trans = Z_options(:,o)-curr_intersect;
           else
               trans = [0;0];
           end

           p_approx = P_options(:,o) - trans;
           q_approx = Q_options(:,o) - trans;
           z_approx = Z_options(:,o) - trans;
           [R,t,flag] = calcPose(p1,p2,p3,p_approx,q_approx,z_approx);
            if (flag)
                R_approx = [R_approx;R];
                t_approx = [t_approx;t];
            end
        end
        
    % Case(III)
    % If l1 and l2 are parallel and p2 has two options to reach l2
    elseif (abs(b2-b1) < r12)
        % Set p to be some point on l1, and find the two possible locations of
        % point q (by finding the intersection between l2 and the
        % circle with radius r12.
        [qx,qy] = linecirc(slope_2,intersect_2y,p_l1(1),p_l1(2),r12);
        if (isinf(slope_2))
            qx = [p_l2(1);p_l2(1)];
            qy = [p_l1(2)+sqrt(r12^2-(b2-p_l1(1))^2);p_l1(2)-sqrt(r12^2-(b2-p_l1(1))^2)];
        end

        % For each of the 2 possible locations for q on l2
        for currq = 1:length(qx)
           [zx,zy] = circcirc(p_l1(1),p_l1(2),r13,qx(currq),qy(currq),r23);
           if (sum(isnan(zx)+isnan(zy)) > 0)    % If answer contains NaNs (3 points are on a line)
               q = [qx(currq);qy(currq)];

               pq_vec = q-p_l1;
               pq_slope = pq_vec(2)/pq_vec(1);
               pq_intersect_y = p_l1(2)-pq_slope*p_l1(1);

               [zx,zy] = linecirc(pq_slope,pq_intersect_y,p_l1(1),p_l1(2),r13);

               if (abs(norm([zx(1);zy(1)]-q)-r23) < small_number)
                   zx = zx(1);
                   zy = zy(1);
               else
                   zx = zx(2);
                   zy = zy(2);
               end 
           end
           
           % For each of the 2 possible locations of z for the current p and q
           for o=1:length(zx)
               curr_intersect = lineLineIntersection(v3,p_l3,v1,[zx(o);zy(o)]);

               % If z cannot intersect it's line then dont move at all
               isBadNumber = sum(isnan(curr_intersect))+sum(isnan(curr_intersect));
               if (isBadNumber == 0)
                   trans = [zx(o);zy(o)]-curr_intersect;
               else
                   trans = [0;0];
               end

               p_approx = p_l1 - trans;
               q_approx = [qx(currq);qy(currq)] - trans;
               z_approx = [zx(o);zy(o)] - trans;
               
               [R,t,flag] = calcPose(p1,p2,p3,p_approx,q_approx,z_approx);
                if (flag)
                    R_approx = [R_approx;R];
                    t_approx = [t_approx;t];
                end
            end
        end 
        
    % Case(IV)
    % If l1 and l2 are parallel and q has only 1 critical point (1 closest position to l2)
    else
        % Find the position of q which is perp to v1 (and v2) and is the closest to l2
        q1 = p_l1+r12*v1_bot;
        q2 = p_l1-r12*v1_bot;

        if (norm(v2'*q1-b2) < norm(v2'*q2-b2))
            q = q1;
        else
            q = q2;
        end
        
        % Find Z's positions
        [zx,zy] = circcirc(p_l1(1),p_l1(2),r13,q(1),q(2),r23);

        if (sum(isnan(zx)+isnan(zy)) > 0)    % If answer contains NaNs (3 points are on a line)
%            disp('here!'); 
           pq_vec = q-p_l1;
           pq_slope = pq_vec(2)/pq_vec(1);
           pq_intersect_y = p_l1(2)-pq_slope*p_l1(1);
           [zx,zy] = linecirc(pq_slope,pq_intersect_y,p_l1(1),p_l1(2),r13);
           
           if (abs(norm([zx(1);zy(1)]-q)-r23) < small_number)
               zx = zx(1);
               zy = zy(1);
           else
               zx = zx(2);
               zy = zy(2);
           end
        end
        
        for o=1:length(zx)
           curr_intersect = lineLineIntersection(v3,p_l3,v1,[zx(o);zy(o)]);

           % If z cannot intersect it's line then dont move at all
           isBadNumber = sum(isnan(curr_intersect))+sum(isnan(curr_intersect));
           if (isBadNumber == 0)
               trans = [zx(o);zy(o)]-curr_intersect;
           else
               trans = [0;0];
           end

           p_approx = p_l1 - trans;
           q_approx = q - trans;
           z_approx = [zx(o);zy(o)] - trans;
           [R,t,flag] = calcPose(p1,p2,p3,p_approx,q_approx,z_approx);
            if (flag)
                R_approx = [R_approx;R];
                t_approx = [t_approx;t];
            end
        end
    end

    if (plotFlag == 1)
        figure;
        numAlignments = size(R_approx,1)/2;
        for i=1:numAlignments
           curr_fig = subplot(1,numAlignments,i); hold on
           plotResults2([p1,p2,p3],[v1,v2,v3],[b1,b2,b3],R_approx(2*i-1:2*i,:),t_approx(2*i-1:2*i),curr_fig);
        end
    end
end