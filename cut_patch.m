
function abs_sum_result = cut_patch(A_o,window_size)

% 假设你的矩阵为A
% A_o = randi([1, 10], 344, 346); % 这里使用了一个随机生成的矩阵，你可以替换为你的实际矩阵

% 定义窗口大小
% window_size = 3;

% 利用im2col函数生成滑动窗口的列
A = padarray(A_o, [(window_size-1)/2, (window_size-1)/2], 0, 'both');
columns = im2col(A, [window_size, window_size], 'sliding');

% 计算列的总数，即生成的小矩阵数量
num_columns = size(columns, 2);

% 将列还原为小矩阵
small_matrices = reshape(columns, [window_size, window_size, num_columns]);


abs_sum_result = mean(abs(small_matrices), [1 2]);

abs_sum_result = reshape(abs_sum_result,size(A_o,1),size(A_o,2));

