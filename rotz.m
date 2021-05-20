function Rz = rotz(theta) 
% (* rotation matrices; to rotate a 3d point p_i; 
% right hand rule thumb pointint to +z direction 
% we left multiply the rotation matrix to it: p_I= RIi *p_i  *)

c = cos(theta);
s = sin(theta); 

Rz = [ c -s 0 
    s c 0 
     0 0 1]; 

end
