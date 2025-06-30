figure
hold on;
axis equal;
for i=1:size(Data_linear,2)
    d = Data_linear{i};
    plot(d(1,:), d(2,:));
    
end

figure
hold on;
axis equal;
for i=1:size(Data_linear,2)
    d = Data_linear{i};
    plot(d(3,:), d(4,:));
    
end
