function [Z,Zs,Z2,Z3] = calcZ(v1,v2,r1,r2,r3)
    % Given 2 unit vector and r1,r2,r3 \in \REAL that satisfy the triangle
    % inequality and r1 > 0, calculate the matrix Z as in lemma Zx
    
    v1_bot = null(v1');
    if (dot(v2,v1_bot) < 0)
       v1_bot = -1*v1_bot; 
    end
    d1 = (r1^2+r2^2-r3^2)/(2*r1);
    d2 = sqrt(r2^2-d1^2);
    s = dot(v1,v2)/sqrt(1-dot(v1,v2)^2);
    Z1 = [0,-1;1,0];
    Z2 = r1*[s*v1, v1];
    Z3 = r1*[s*v1+v1_bot,[0;0]];
    Z = Z2+(d1/r1)*(Z3-Z2)+(d2/r1)*Z1*(Z3-Z2);
    Zs = Z2+(d1/r1)*(Z3-Z2)-(d2/r1)*Z1*(Z3-Z2);
end