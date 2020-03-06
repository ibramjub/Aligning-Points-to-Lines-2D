function [P_trig,Q_trig,Z_trig] = Trig(r1,r2,r3,vj,vk,bj,bk)
    % Given r1,r2,r3>0 and two lines defined by two unit vectors (the line directions) and two
    % scalars, find all possible options for a point z which satisfies:
    % dist(z,p)=r2,dist(z,q)=r3,dist(p,q)=r1,p \in lj q \in lk,
    % dot((p-q),vk) = 0;    
    small_number = 0.01;
    plotData = 0;
    P_trig = [];
    Q_trig = [];
    Z_trig = [];
    
    vj_perp = null(vj');
    vk_perp = null(vk');
    
    cos_theta = dot(vj,vk);
    sin_theta = sqrt(1-cos_theta^2);
    
    % Distance between the intersection point and p (the intersection,p,q
    % form a right triangle)
    a = r1/sin_theta;  
    b = a*cos_theta;
    assert(r1^2+b^2-a^2 < small_number);
    
    % The intersection point of lj and lk
    inter_jk = lineLineIntersection(vj,bj*vj_perp,vk,bk*vk_perp);
    
    % The intersection of lj,lk and the y-axis
    inter_jy = lineLineIntersection(vj,bj*vj_perp,[0;1],[0;0]);
    inter_ky = lineLineIntersection(vk,bk*vk_perp,[0;1],[0;0]);
    
    % The slopes of lj and lk
    slope_k = vk(2)/vk(1);
    slope_j = vj(2)/vj(1);
        
    % Find the two positions for point p and two possible positions for q
%     [px,py] = linecirc(slope_j,inter_jy(2),inter_jk(1),inter_jk(2),a);
%     [qx,qy] = linecirc(slope_k,inter_ky(2),inter_jk(1),inter_jk(2),b);
    p1 = inter_jk+a*vj;
    p2 = inter_jk-a*vj;
    px = [p1(1),p2(1)];
    py = [p1(2),p2(2)];
    q1 = inter_jk+b*vk;
    q2 = inter_jk-b*vk;
    qx = [q1(1),q2(1)];
    qy = [q1(2),q2(2)];
    
    % Find right matching between two options of p and two options of q
    if (norm([qx(1);qy(1)]-[px(1);py(1)])-r1 > small_number)
       qx([1,2]) = qx([2,1]); 
       qy([1,2]) = qy([2,1]); 
    end

    for i=1:length(px)
        if (r3 ~= 0)
            [zx,zy] = circcirc(px(i),py(i),r2,qx(i),qy(i),r3);
        else  % If third point is exactly the second point
            zx = qx(i);
            zy = qy(i);
        end
        for o=1:length(zx)
            P_trig = [P_trig,[px(i);py(i)]];
            Q_trig = [Q_trig,[qx(i);qy(i)]];
            Z_trig = [Z_trig,[zx(o);zy(o)]];
        end
    end

    % Plots
    if (plotData)
        l1 = [inter_jk';inter_jk'+vj'];
        l2 = [inter_jk';inter_jk'+vk'];
        origin = [0;0];
        plot(l1(:,1),l1(:,2),'color','b'); hold on
        plot(l2(:,1),l2(:,2),'color','b'); hold on
        plot(origin(1),origin(2),'+','color','b'); hold on;
    end
    
    for i=1:size(Z_trig,2)
       assert(abs(vj_perp'*P_trig(:,i)-bj) < small_number);
       assert(abs(vk_perp'*Q_trig(:,i)-bk) < small_number);
       assert(abs(norm(P_trig(:,i)-Q_trig(:,i))-r1) < small_number);
       assert(abs(norm(Z_trig(:,i)-P_trig(:,i))-r2) < small_number);
       assert(abs(norm(Z_trig(:,i)-Q_trig(:,i))-r3) < small_number);
       assert(abs(dot(P_trig(:,i)-Q_trig(:,i),vk)) < small_number);
       if (plotData == 1)
           plot(P_trig(1,i),P_trig(2,i),'+','color','red'); hold on
           plot(Q_trig(1,i),Q_trig(2,i),'+','color','yellow'); hold on
           plot(Z_trig(1,i),Z_trig(2,i),'+','color','green'); hold on
       end
    end
end