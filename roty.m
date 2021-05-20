function Ry = roty(theta) 
% (* rotation matrices; to rotate a 3d point p_i; 
% right hand rule thumb points to +y direction we left \
% multiply the rotation matrix to it: p_I= RIi *p_i  *)

c = cos(theta);
s = sin(theta); 

Ry = [ c 0 s 
    0 1 0 
    -s 0 c]; 

end
