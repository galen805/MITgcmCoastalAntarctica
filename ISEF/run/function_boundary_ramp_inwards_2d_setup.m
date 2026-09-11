function output_2d_matrix = function_boundary_ramp_inwards_2d_setup(bathy_matrix,boundary_rest_value,center_value,const_rest_length,ramp_to_center_length)

% From the boundary to the inside, inward direction
% first, we have const values along west/east/north boundaries
% then, there is a ramp inwards to the center value

% Example
% bathy_matrix = ones(360,480);
% boundary_rest_value = 270;
% center_value = 200;
% const_rest_length = 10;
% ramp_to_center_length = 10;


output_2d_matrix = center_value*ones(size(bathy_matrix));

ramp_vector = linspace(boundary_rest_value, center_value, ramp_to_center_length);

for k = 1:const_rest_length
    output_2d_matrix(1:const_rest_length,:) = boundary_rest_value;
    output_2d_matrix((end-const_rest_length):end,:) = boundary_rest_value;
    output_2d_matrix(:,(end-const_rest_length):end) = boundary_rest_value;
end

for k = 1:ramp_to_center_length
    output_2d_matrix(const_rest_length+k,1:(end-const_rest_length-k+1)) = ramp_vector(k);
    output_2d_matrix((end-const_rest_length-k+1),1:(end-const_rest_length-k+1)) = ramp_vector(k);
    output_2d_matrix((const_rest_length+k):(end-const_rest_length-k+1),(end-const_rest_length-k+1)) = ramp_vector(k);
end