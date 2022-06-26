function angle=angfor2vector(vector1,vector2)
angle=atan2(norm(cross(vector1,vector2)), dot(vector1,vector2));
angle=angle/2/pi*360;
end