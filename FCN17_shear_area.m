function downstand_sum = FCN17_shear_area(downstand_sum,Perimeter_num)

for i0 = 1:Perimeter_num
    ele_tem = downstand_sum(i0).All_element;
    num_ele = length(ele_tem(:,3));
    area_section = sum(ele_tem(:,3).*(ones(num_ele,1)-ele_tem(:,4)).*ele_tem(:,5));
    downstand_sum(i0).area = area_section;
end