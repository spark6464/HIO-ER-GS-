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
params.fontSize=200; % 设置数字大小
binaryImg = createMask('digit', [512, 512], params);

% 画出binaryImg
figure;imshow(binaryImg); %一个512x512的逻辑数组(0,1)

% ============ 添加像素标尺 ============

scale_px = 100;  % 标尺长度：100 像素
x0 = 50;         % 标尺起始横坐标
y0 = size(binaryImg,1) - 50;  % 标尺纵坐标（放在下边缘）

% 画标尺线段
line([x0 x0 + scale_px], [y0 y0], 'Color','w', 'LineWidth', 3);

% 添加文字说明
text(x0 + scale_px/2, y0 - 30, '100 pixels', ...
     'Color','w', 'HorizontalAlignment','center', 'FontSize', 10);

%% 图的傅立叶变换（FT），频域振幅（频谱）
fft_binaryImg=abs(fftshift(fft2(binaryImg)));
figure;
% imagesc(fft_binaryImg);
imagesc(log(1+fft_binaryImg));%对数增强视觉效果
axis off;
% colorbar;

%% 图的自相关 （图做FT，再算模值的平方，即功率谱，再逆FT）（Wiener-Khinchin:自相关与功率谱互为FT）
autocorr_binaryImg=fftshift(real(ifft2(abs(fft2(binaryImg)).^2)));%FT法
figure;
imagesc(autocorr_binaryImg );
%imagesc(log(1+autocorr_binaryImg));%对数增强视觉效果
colormap([zeros(256, 1), linspace(0, 1, 256)', zeros(256, 1)]);
axis off;
colorbar;

% ============ 添加像素标尺 ============

scale_px = 100;  % 标尺长度：100 像素
x0 = 50;         % 标尺起始横坐标
y0 = size(binaryImg,1) - 60;  % 标尺纵坐标（放在下边缘）

% 画标尺线段
line([x0 x0 + scale_px], [y0 y0], 'Color','w', 'LineWidth', 3);

% 添加文字说明
text(x0 + scale_px/2, y0 - 30, '100 pixels', ...
     'Color','w', 'HorizontalAlignment','center', 'FontSize', 20);

%% 自相关的FT/功率谱 （图的自相关做FT，得到功率谱）
autocorr_binaryImg=fftshift(real(ifft2(abs(fft2(binaryImg)).^2)));%先求图的自相关
power_spectrum_binaryImg=abs(fftshift(fft2(autocorr_binaryImg)));%根据自相关做FT求功率谱
figure;
% imagesc(power_spectrum_binaryImg);
imagesc(log(1+power_spectrum_binaryImg));%对数增强视觉效果
axis off;
colorbar;

%% 自相关的FT开根号/频谱 （图的自相关做FT，得到功率谱，开根号得到频域振幅（频谱））
autocorr_binaryImg=fftshift(real(ifft2(abs(fft2(binaryImg)).^2)));%先求图的自相关
power_spectrum_binaryImg=abs(fftshift(fft2(autocorr_binaryImg)));%根据自相关做FT求功率谱
fourier_amplitute_binaryImg=sqrt(power_spectrum_binaryImg);%开根号得到频谱/频域振幅
figure;
% imagesc(fourier_amplitute_binaryImg);
imagesc(log(1+fourier_amplitute_binaryImg));%对数增强视觉效果
axis off;
colorbar;

%% HIO-ER相位恢复
% 
% % 初始化相位恢复的输入频域振幅，可以是图像自相关的FT开根号/图像直接FT
 mag_spec = fourier_amplitute_binaryImg; 
% mag_spec = fft_binaryImg; % 数字板的FT是频域振幅，不用开根号

% 迭代次数
iter_num = 90;  

% 初始化随机相位
g = rand(512, 512);

% 窗口参数
img_size = 512;  % 原始图像尺寸
win_size = 500; % 窗口尺寸

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
    imshow(g(1:win_size, 1:win_size), 'InitialMagnification', 200);
    title(strcat('迭代步数 ', num2str(i)));
    pause(0.01);  % 短暂暂停以观察恢复过程
end

%%频域投影函数
function g2 = projM(g, mag_spec)
    G = fftshift(fft2(g));   % 计算傅里叶变换
    G2 = mag_spec .* (G ./ abs(G)); % 仅调整幅度，保持相位
    g2 = ifft2(ifftshift(G2));  % 逆变换回空间域
end

