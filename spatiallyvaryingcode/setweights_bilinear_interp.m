
%%

%version = beadcalib;
%if simulation; version = 1; end

[Ny, Nx, Py, Px, Nz] = size(psfsr);

weightbilin = zeros(Ny, Nx, Py, Px);
x1 = 1:Nx; x2 = 1:Ny;
[X1, X2] = meshgrid(x1, x2);

% if version==1
%     rangei = 1:(Py-1);
%     rangej = 1:(Px-1);
% elseif version == 2
    rangei = 1:(Py-1);
    rangej = 2:(Px);
% end


for i = rangei
    fprintf('%d / %d \n', i, Py-1);
    for j = rangej
        
        % find nearest 4 points
            % if beads = v1
%         if version == 1
%         x11 = xpos(i, j); y11 = ypos(i, j);
%         x12 = xpos(i, j+1); y12 = ypos(i, j+1);
%         x21 = xpos(i+1, j); y21 = ypos(i+1, j);
%         x22 = xpos(i+1, j+1); y22 = ypos(i+1, j+1);
%         elseif version == 2
            % if beads = v2
        x11 = xpos(i, j); y11 = ypos(i, j);
        x12 = xpos(i, j-1); y12 = ypos(i, j-1);
        x21 = xpos(i+1, j); y21 = ypos(i+1, j);
        x22 = xpos(i+1, j-1); y22 = ypos(i+1, j-1);
%         end
        
        targetarea = zeros(Ny, Nx);
        targetarea(X2 <= (y12-y11)/(x12-x11)*(X1 - x11) + y11 & ...
            X2 > (y21-y22)/(x21-x22)*(X1 - x22) + y22 & ...
            X1 <= (x21-x11)/(y21-y11)*(X2 - y11) + x11 & ...
            X1 > (x12-x22)/(y12-y22)*(X2 - y22) + x22 ...
            ) = 1;
        
        Ay = X2; if round(x21-x11) == 0; Ax = ones(size(X1))*x11; else Ax = (x21-x11)./(y21-y11).*(X2-y11)+ x11; end
        Cy = X2; if round(x12-x22) == 0; Cx = ones(size(X1))*x22; else Cx = (x12-x22)./(y12-y22).*(X2-y22)+ x22; end
        
        A11dist = sqrt((Ax-x11).^2 + (Ay-y11).^2);
        A21dist = sqrt((Ax-x21).^2 + (Ay-y21).^2);
        A11 = A21dist./(A11dist+A21dist);
        A21 = A11dist./(A11dist+A21dist);
        
        
        C12dist = sqrt((Cx-x12).^2 + (Cy-y12).^2);
        C22dist = sqrt((Cx-x22).^2 + (Cy-y22).^2);
        C22 = C12dist./(C22dist+C12dist);
        C12 = C22dist./(C22dist+C12dist);
        
        Apdist = sqrt((Ax-X1).^2+(Ay-X2).^2);
        Cpdist = sqrt((Cx-X1).^2+(Cy-X2).^2);
        
        p11 = A11.*(Cpdist./(Apdist+Cpdist));
        p21 = A21.*(Cpdist./(Apdist+Cpdist));
        p12 = C12.*(Apdist./(Apdist+Cpdist));
        p22 = C22.*(Apdist./(Apdist+Cpdist));
        
%         if version == 1
%         weightbilin(:, :, i, j) = p11.*targetarea + weightbilin(:, :, i, j);
%         weightbilin(:, :, i, j+1) = p12.*targetarea + weightbilin(:, :, i, j+1);
%         weightbilin(:, :, i+1, j) = p21.*targetarea + weightbilin(:, :, i+1, j);
%         weightbilin(:, :, i+1, j+1) = p22.*targetarea + weightbilin(:, :, i+1, j+1);
%         elseif version==2
        weightbilin(:, :, i, j) = p11.*targetarea + weightbilin(:, :, i, j);
        weightbilin(:, :, i, j-1) = p12.*targetarea + weightbilin(:, :, i, j-1);
        weightbilin(:, :, i+1, j) = p21.*targetarea + weightbilin(:, :, i+1, j);
        weightbilin(:, :, i+1, j-1) = p22.*targetarea + weightbilin(:, :, i+1, j-1);
%         end
            
    end
end
%%
weightsds = zeros(Ny, Nx, Py*Px);
for i = 1:Py
    for j = 1:Px
        weightsds(:, :, j+Px*(i-1)) = weightbilin(:,:,i,j);
        
    end
end

pcds = pc_all;
disp('Done creating bilinear weights')



