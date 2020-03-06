%{
"""*****************************************************************************************
MIT License
Copyright (c) 2020 Ibrahim Jubran, Dan Feldman
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*****************************************************************************************"""
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% This code is based on Algorithm 2 from the paper:
% Aligning Points to Lines: Provable Approximations
% By: Ibrahim Jubran and Dan Feldman
% Please cite the paper when using the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R_approx,t_approx] = AlignPointsToLines(P,L,B)
    % Align the set of points given in P (2XN) to the set of lines, where
    % the i'th line is defined by its direction vector L(:,i) and its
    % distance from the origin B(i).
    
    
    % Num of points
    n = length(P);
    R_approx = [];
    t_approx = [];
    
    counter = 1;

    jkl_indices = permn(1:n,3);
    for t = 1:length(jkl_indices)
        j = jkl_indices(t,1);
        k = jkl_indices(t,2);
        l = jkl_indices(t,3);
        if (length(unique([j,k,l]))~=3)
            continue
        end
        
        counter = counter + 1;
        try
            [curr_R,curr_t] = ComputeAlignment(P(:,j),P(:,k),P(:,l),L(:,j),L(:,k),L(:,l),B(j),B(k),B(l));
            R_approx = [R_approx;curr_R];
            t_approx = [t_approx;curr_t];
        catch

        end
    end
end