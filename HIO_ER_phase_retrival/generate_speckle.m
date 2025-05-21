function speckle_intensity = generate_speckle(varargin)
    % 生成散斑图像的函数
    % 示例调用：
    % speckle_intensity = generate_speckle('N', 2048, 'W', 20, 'S', 2, 'D', 128);
    
    % 解析输入参数
    p = inputParser;
    addParameter(p, 'N', 2048); % 二维数组大小
    addParameter(p, 'W', 20);   % 圆形平滑滤波器直径
    addParameter(p, 'S', 2);    % 相位乘法因子
    addParameter(p, 'D', 128);  % 光阑直径
    parse(p, varargin{:});
    
    N = p.Results.N;
    W = p.Results.W;
    S = p.Results.S;
    D = p.Results.D;
    
    %% 创建卷积核（低通滤波器）
    kernel = zeros(N, N);
    [x, y] = meshgrid(1:N, 1:N);
    radius = sqrt((x - N/2).^2 + (y - N/2).^2);
    kernel(radius < W/2) = 1 / sqrt(pi * (W/2)^2);

    %% 生成无关相位数组（随机相位场）
    uncorrelated_phase = S * pi * randn(N, N);
    
    %% 进行傅里叶变换卷积，得到相关相位场
    correlated_phase = ifft2(fft2(uncorrelated_phase) .* fft2(kernel));
    
    %% 计算扩散器传输的场（复指数形式）
    diffuser = exp(1i * correlated_phase);
    
    %% 创建圆形光阑（限制光束传播）
    stop = zeros(N, N);
    stop(radius < D/2) = 1;
    
    %% 计算入射光场（散射光的傅里叶变换）
    incident_field = fftshift(fft2(diffuser));
    transmitted_field = incident_field .* stop;
    
    %% 计算图像平面上的光强分布
    speckle_field = 1/N * ifft2(transmitted_field);
    speckle_intensity = abs(speckle_field).^2;
    speckle_intensity = speckle_intensity / max(speckle_intensity(:));
   
end
