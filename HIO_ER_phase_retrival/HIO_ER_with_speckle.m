%% 调用generate_speckle产生散斑图
speckle_intensity = generate_speckle('N', 2048, 'W', 2, 'S', 2, 'D',2048);
% 'N':二维数组大小;'W':圆形平滑滤波器直径;'S':% 相位乘法因子(产生2pi);'D':光阑直径

%注意：光阑大小对傅立叶振幅有影响，会影响恢复结果。光阑太小会失去部分高频信号，恢复结果较差。

% 绘制散斑图像
figure;
set(gcf, 'Position', [100, 100, 550, 550]);
imagesc(speckle_intensity);
colormap([zeros(256, 1), linspace(0, 1, 256)', zeros(256, 1)]);
caxis([0, 0.5]);
axis off;
% colorbar;



% ========================
% 添加 400 px 的标尺
scale_px = 400;                          % 标尺长度：400像素
x0 = 200;                                % 起点横坐标
y0 = size(speckle_intensity,1) - 200;    % 起点纵坐标（靠近下边）
% 绘制白色线段
line([x0, x0 + scale_px], [y0, y0], 'Color', 'w', 'LineWidth', 4);
% 添加注释文字
text(x0 + scale_px/2, y0 - 120, '400 pixels', ...
    'Color', 'w', 'FontSize', 30, 'HorizontalAlignment', 'center');


%% 生成数字板/双缝/三缝图案
% 生成矩形双缝/三缝
% 竖直双缝（左右排列）名称：'double_slit_vertical'
% 水平双缝（上下排列）名称：'double_slit_horizontal'
% 竖直三缝（左右排列）名称：'triple_slit_vertical'
% 水平三缝（上下排列）名称：'triple_slit_horizontal'
% params.width = 30;  % 单个竖孔的宽度
% params.height = 200;  % 竖孔的高度
% params.gap = 30;  % 竖孔间距
% binaryImg = createMask('double_slit_vertical', [512, 512], params);

% 生成数字板
params.text='7'; % 选择数字
params.fontSize=80; % 设置数字大小
binaryImg = createMask('digit', [512, 512], params);

% 画出binaryImg 
figure;imshow(binaryImg); %一个512x512的逻辑数组(0,1)



% ============ 添加像素标尺 ============
scale_px = 50;  % 标尺长度：100 像素
x0 = 50;         % 标尺起始横坐标
y0 = size(binaryImg,1) - 50;  % 标尺纵坐标（放在下边缘）
% 画标尺线段
line([x0 x0 + scale_px], [y0 y0], 'Color','w', 'LineWidth', 3);
% 添加文字说明
text(x0 + scale_px/2, y0 - 30, '50 pixels', ...
     'Color','w', 'HorizontalAlignment','center', 'FontSize', 10);

%% 计算卷积得到散斑模糊图像

%第一种方式计算卷积：conv2函数
% blurred_speckle = conv2(speckle_intensity, double(binaryImg), 'same');

%第二种方式计算卷积：FT法
% 计算傅里叶变换
I_fft = fft2(speckle_intensity);
binaryImg_fft = fft2(binaryImg, 2048, 2048); % 直接填充到2048x2048大小
% 频域相乘（卷积定理）
blurred_speckle = abs(fftshift(ifft2(I_fft .* binaryImg_fft)));
%归一化结果
blurred_speckle = blurred_speckle / max(blurred_speckle(:));
blurred_speckle = blurred_speckle - mean(blurred_speckle(:));


% 绘制卷积后的散斑图像
figure;
imagesc(blurred_speckle);
colormap([zeros(256, 1), linspace(0, 1, 256)', zeros(256, 1)]);
axis off;



% ========================
% 添加 400 px 的标尺
scale_px = 400;                          % 标尺长度：400像素
x0 = 200;                                % 起点横坐标
y0 = size(speckle_intensity,1) - 200;    % 起点纵坐标（靠近下边）
% 绘制白色线段
line([x0, x0 + scale_px], [y0, y0], 'Color', 'w', 'LineWidth', 4);
% 添加注释文字
text(x0 + scale_px/2, y0 - 120, '400 pixels', ...
    'Color', 'w', 'FontSize', 30, 'HorizontalAlignment', 'center');


%% 自相关的FT开根号/频谱 （图的自相关做FT，得到功率谱，开根号得到频域振幅（频谱））
autocorr_binaryImg=fftshift(real(ifft2(abs(fft2(blurred_speckle)).^2)));%先求图的自相关
% 取自相关中心区域
center_crop = autocorr_binaryImg(1024-255:1024+256, 1024-255:1024+256); % crop出 


% 然后做FT再开根号
power_spectrum_binaryImg=abs(fftshift(fft2(center_crop)));%根据自相关做FT求功率谱
fourier_amplitute_binaryImg=sqrt(power_spectrum_binaryImg);%开根号得到频谱/频域振幅


%画出傅立叶振幅
figure;
% imagesc(fourier_amplitute_binaryImg);
imagesc(log(1+fourier_amplitute_binaryImg));%对数增强视觉效果
colorbar;


%% 高斯滤波
H = fspecial('gaussian', size(fourier_amplitute_binaryImg),0.5);  % 30 可调
filtered_mag = imfilter(fourier_amplitute_binaryImg, H, 'replicate');


%画出高斯滤波后的傅立叶振幅
figure;
% imagesc(filtered_mag);
imagesc(log(1+filtered_mag));%对数增强视觉效果
axis off;


%% HIO-ER相位恢复

% 初始化相位恢复的输入频域振幅，可以是图像自相关的FT开根号/图像直接FT
 mag_spec = filtered_mag; 
% mag_spec = fft_binaryImg; % 数字板的FT是频域振幅，不用开根号

% 迭代次数
iter_num = 90;  

% 初始化随机相位
g = rand(512, 512);

% 窗口参数
img_size = 512;  % 原始图像尺寸
win_size =500; % 窗口尺寸

% 生成 win_size x win_size 矩阵，初始化为 0
win = zeros(img_size, img_size);  
% 设定窗口区域
win(img_size/2 - win_size/2 : img_size/2 + win_size/2, ...
    img_size/2 - win_size/2 : img_size/2 + win_size/2) = 1; 

for i = 1:iter_num
    %=================HIO======================== 效果较好
    g2 = projM(g, mag_spec); % 频域投影操作
    g2 = g2 .* win; % 应用窗函数（支持域约束）
    g = (g2 >= 0) .* g2 + (g2 < 0) .* (g - 0.9 .* g2); % HIO 更新公式

    %=================ER========================
    % g = projM(projSup(g, win), mag_spec);

    % 显示当前恢复的图像
    imshow(abs(g(1:win_size, 1:win_size)), 'InitialMagnification', 200);
   
    
      % ============ 添加像素标尺 ============
       scale_px = 50;  % 标尺长度：100 像素
       x0 = 50;         % 标尺起始横坐标
       y0 = size(binaryImg,1) - 50;  % 标尺纵坐标（放在下边缘）
       % 画标尺线段
        line([x0 x0 + scale_px], [y0 y0], 'Color','w', 'LineWidth', 3);
        % 添加文字说明
         text(x0 + scale_px/2, y0 - 30, '50 pixels', ...
         'Color','w', 'HorizontalAlignment','center', 'FontSize', 10);

    title(strcat('迭代步数 ', num2str(i)));
    pause(0.01);  % 短暂暂停以观察恢复过程
end


%%频域投影函数
function g2 = projM(g, mag_spec)
    G = fftshift(fft2(g));   % 计算傅里叶变换
    G2 = mag_spec .* (G ./ (abs(G) + eps)); % 仅调整幅度，保持相位
    g2 = ifft2(ifftshift(G2));  % 逆变换回空间域
end


