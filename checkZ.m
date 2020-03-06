function Zmat = checkZ(P,Pmat,Qmat,Z1,Z2)
    % Check if the right ellipse matrix is Z1 or Z2 by:
    % 1) place p1 at the origin and p2 at the positive side of the Y-axis
    % 2) place Pmat*x at the origin and Qmat*x at the positive side of the Y-axis
    % 3) pick the matrix Zi (i=1,2) such that ||Zi*x-p3|| is the smallest
    
    small_number = 0.0001;
    
    % random unit vector
    x = rand(2,1);
    x = x/norm(x);
    
    % Align p1_new with the origin and p2_new with the positive side of the Y-axis
    p2_new = Qmat*x-Pmat*x;
    p3_new = [null(p2_new')';p2_new']*(Z1*x-Pmat*x);
    
    % Align p1 with the origin and p2 with the positive side of the Y-axis
    p2 = P(:,2)-P(:,1);
    p3 = [null(p2')';p2']*(P(:,3)-P(:,1));
    
    if (norm(p3-p3_new) <= small_number)
        Zmat = Z1;
    else
        Zmat = Z2;
    end
    
    p3_final = [null(p2_new')';p2_new']*(Zmat*x-Pmat*x);
    assert(norm(p3-p3_final) <= small_number);
end