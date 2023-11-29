function patches = imagePatchSplit(image, npatch)
    % 输入参数：
    % image: 输入图像
    % npatch: 切割成 n * n 块
    
    % 获取图像尺寸
    [rows, cols] = size(image);
    
    % 计算每个块的尺寸
    patchSize = floor([rows, cols] / npatch);
    
    % 确保每个块的尺寸是整数
    patchSize = patchSize + rem([rows, cols], npatch * patchSize) / npatch;
    
    % 预分配单元格数组
    patches = cell(npatch, npatch);
    
    % 切割图像并存储到单元格数组中
    for i = 1:npatch
        for j = 1:npatch
            rowStart = 1 + (i - 1) * patchSize(1);
            rowEnd = min(i * patchSize(1), rows);
            
            colStart = 1 + (j - 1) * patchSize(2);
            colEnd = min(j * patchSize(2), cols);
            
            patches{i, j} = image(rowStart:rowEnd, colStart:colEnd);
        end
    end
end
