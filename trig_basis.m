function trig_Mat = trig_basis(OM)
%Returns sin, cos, tan of the original matrix
trig_Mat = [sin(OM) cos(OM) tan(OM)];
% Alternative sin-cos-tan of each column order
%     trig_Mat = [];
%     for i = 1:size(OM,2) %loop through columns
%         trig_Mat = [trig_Mat sin(OM(:,i)) cos(OM(:,i)) tan(OM(:,i))];
%     end
end