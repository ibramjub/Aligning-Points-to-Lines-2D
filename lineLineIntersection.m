function intersection = lineLineIntersection(v1,o1,v2,o2)
    % Given two lines defined by a unit vector and point, find the
    % intersection points between the two lines
    
    % p1 is a another point on l1, p2 is another point on l2
    p1 = o1+v1;
    p2 = o2+v2;
    
    x1 = o1(1);
    x2 = p1(1);
    x3 = o2(1);
    x4 = p2(1);
    
    y1 = o1(2);
    y2 = p1(2);
    y3 = o2(2);
    y4 = p2(2);
    
    denominator = det([det([x1 1; x2 1]), det([y1 1; y2 1]); det([x3 1; x4 1]), det([y3 1; y4 1])]);
    x_inter = det([det([x1 y1; x2 y2]), det([x1 1; x2 1]); det([x3 y3; x4 y4]), det([x3 1; x4 1])])/denominator;
    y_inter = det([det([x1 y1; x2 y2]), det([y1 1; y2 1]); det([x3 y3; x4 y4]), det([y3 1; y4 1])])/denominator;
    
    intersection = [x_inter; y_inter];
end