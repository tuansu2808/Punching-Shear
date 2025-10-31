function Hole = FCN2_hole_processing(number_of_hole_input,...
                                    x_hole_input, y_hole_input, ...
                                    Check_radius)

number_of_hole = number_of_hole_input;
% theta_fraction = 0.0001;

Hole = struct;
for i1 = 1:number_of_hole
    Hole(i1).x_hole = x_hole_input(i1,:);
    Hole(i1).y_hole = y_hole_input(i1,:);
end

for i1 = 1:number_of_hole

    x_hole = Hole(i1).x_hole;
    y_hole = Hole(i1).y_hole;

    
    x_far_1 = x_hole(1)*Check_radius*1.2;
    y_far_1 = y_hole(1)*Check_radius*1.2;
    x_far_2 = x_hole(2)*Check_radius*1.2;
    y_far_2 = y_hole(2)*Check_radius*1.2;


    shape = polyshape([0,x_far_1,x_far_2],[0,y_far_1,y_far_2]);

    Hole(i1).shape = shape;
end