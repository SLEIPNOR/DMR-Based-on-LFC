function MatCrop = ZeroCrop(IM_mat)

% This function crops images after applying scaled and rotation recovering
% It should be invoked as:
%
%   IC = ZeroCrop(IM_mat)
%
% where 
%   MatCrop, is image group after cropping an with no zero point 
%           in image.
%           
%   IM_mat, is the input image group, which should be a tensor x*y*n.
%
%   Method, is first finding the row and col whose pixels are all zero,
%           and then finding the optimal cropping way to generate cropped
%           images.
%  Zhenzhe Han
%  12/2022


% setting parameters
s = size(IM_mat);
a = s(1);
b = s(2);
len = s(3);
size_row_min = 1;
size_row_max = a ;
size_col_min = 1 ;
size_col_max = b ;
% All zeros pixels rows and cols searching

for img_No = 1:len

    for row = 1:a

        zero_num = sum(IM_mat(row,:,img_No)==0);

        if zero_num ~=b

            size_row_min = max(row,size_row_min);
            break
    
        end

    end

end 



for img_No = 1:len

    for row = a:-1:1

        zero_num = sum(IM_mat(row,:,img_No)==0);

        if zero_num ~= b

            size_row_max = min(row,size_row_max);
            break
    
        end

    end

end 



for img_No = 1:len

    for col = 1:b

        zero_num = sum(IM_mat(:,col,img_No)==0);

        if zero_num ~= a

            size_col_min = max(col,size_col_min);
            break
    
        end

    end

end 



for img_No = 1:len

    for col = b:-1:1

        zero_num = sum(IM_mat(:,col,img_No)==0);

        if zero_num ~= a

            size_col_max = min(col,size_col_max);
            break
    
        end

    end

end 

% first cropping applied
MatCrop = IM_mat(size_row_min:size_row_max,size_col_min:size_col_max,:);

% finding the optimal way for second cropping
row=[];
col=[];

for img_No = 1:len

    [row_temp,col_temp] = find(MatCrop(:,:,img_No)==0);
    row = [row;row_temp];
    col = [col;col_temp];
    

end


del_rows=[];
del_cols=[];
while (~isempty(row))
    max_row = max(mode(row));
    pos_row = find(row==max_row);
    max_col = max(mode(col));
    pos_col = find(col==max_col);

    if length(pos_row) >= length(pos_col)

        del_rows = [del_rows;max_row];
        row(pos_row,:) = [];
        col(pos_row,:) = [];
        
    else 
        del_cols = [del_cols;max_col];       
        col(pos_col,:) = [];
        row(pos_col,:) = [];
    end
end

MatCrop(del_rows,:,:)=[];
MatCrop(:,del_cols,:)=[];


end