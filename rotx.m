function Rx = rotx(theta)
% (* rotation matrices; to rotate a 3d point p_i;
% right hand rule thumb pointint to x direction we left \
% multiply the rotation matrix to it: p_I= RIi *p_i  *)

Rx = [ 1 0 0
    0 cos(theta) -sin(theta)
    0 sin(theta) cos(theta)];

end
