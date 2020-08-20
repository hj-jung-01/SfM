function [Q] = rotm2quat(R)
    % Preallocation of variables.
    Q = zeros(1, 4);
    q0 = 0;
    q1 = 0;
    q2 = 0;
    q3 = 0;
    
    % Convert 3x3 rotation matrix to quaternion.
    % 1. Find largest squared element.
    q0Square = 0.25 * (1 + R(1, 1) + R(2, 2) + R(3, 3));
    q1Square = 0.25 * (1 + R(1, 1) - R(2, 2) - R(3, 3));
    q2Square = 0.25 * (1 - R(1, 1) + R(2, 2) - R(3, 3));
    q3Square = 0.25 * (1 - R(1, 1) - R(2, 2) + R(3, 3));
    
    % Solve for quaternion elements.
    if ((q0Square > q1Square) && (q0Square > q2Square) && (q0Square > q3Square))
        q0 = sqrt(q0Square);
        q1 = 0.25 * (R(3, 2) - R(2, 3)) / q0;
        q2 = 0.25 * (R(1, 3) - R(3, 1)) / q0;
        q3 = 0.25 * (R(2, 1) - R(1, 2)) / q0;
    elseif ((q1Square > q0Square) && (q1Square > q2Square) && (q1Square > q3Square))
        q1 = sqrt(q1Square);
        q0 = 0.25 * (R(3, 2) - R(2, 3)) / q1;
        q2 = 0.25 * (R(1, 2) + R(2, 1)) / q1;
        q3 = 0.25 * (R(1, 3) + R(3, 1)) / q1;
    elseif ((q2Square > q0Square) && (q2Square > q1Square) && (q2Square > q3Square))
        q2 = sqrt(q2Square);
        q0 = 0.25 * (R(1, 3) - R(3, 1)) / q2;
        q1 = 0.25 * (R(1, 2) + R(2, 1)) / q2;
        q3 = 0.25 * (R(2, 3) + R(3, 2)) / q2;
    else
        q3 = sqrt(q3Square);
        q0 = 0.25 * (R(2, 1) - R(1, 2)) / q3;
        q1 = 0.25 * (R(1, 3) + R(3, 1)) / q3;
        q2 = 0.25 * (R(2, 3) + R(3, 2)) / q3;
    end

    Q = [q0; q1; q2; q3];
    Q = Q';
end