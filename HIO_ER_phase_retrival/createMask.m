%% 生成各种二值掩膜函数：数字板，条纹图
function binaryImg = createMask(type, img_size, params)
    binaryImg = zeros(img_size);
    switch type
        %% **数字板**
        case 'digit'
            digitText = params.text;
            fontSize = params.fontSize;
            position = [img_size(2)/2, img_size(1)/2];
           RGB = insertText(binaryImg, position, digitText, 'FontSize', fontSize, ...
                 'BoxOpacity', 0, 'TextColor', 'white', 'AnchorPoint', 'Center', ...
                 'Font', 'Arial');  % 更换为支持更多字符的字体
            grayImg = rgb2gray(RGB);
            binaryImg = imbinarize(grayImg);

        %% **竖直双缝（左右排列）**
        case 'double_slit_vertical'
            w = params.width;  % 缝宽
            h = params.height; % 缝高
            gap = params.gap;  % 间隔
            binaryImg(img_size(1)/2-h/2:img_size(1)/2+h/2, img_size(2)/2-gap/2-w:img_size(2)/2-gap/2) = 1;
            binaryImg(img_size(1)/2-h/2:img_size(1)/2+h/2, img_size(2)/2+gap/2:img_size(2)/2+gap/2+w) = 1;

        %% **水平双缝（上下排列）**
        case 'double_slit_horizontal'
            w = params.height;  % 缝宽
            h = params.width; % 缝高
            gap = params.gap;  % 间隔
            binaryImg(img_size(1)/2-gap/2-h:img_size(1)/2-gap/2, img_size(2)/2-w/2:img_size(2)/2+w/2) = 1;
            binaryImg(img_size(1)/2+gap/2:img_size(1)/2+gap/2+h, img_size(2)/2-w/2:img_size(2)/2+w/2) = 1;

        %% **竖直三缝（左右排列）**
        case 'triple_slit_vertical'
            w = params.width;
            h = params.height;
            gap = params.gap;
            binaryImg(img_size(1)/2-h/2:img_size(1)/2+h/2, img_size(2)/2-2*gap/2-w:img_size(2)/2-2*gap/2) = 1;
            binaryImg(img_size(1)/2-h/2:img_size(1)/2+h/2, img_size(2)/2-w/2:img_size(2)/2+w/2) = 1;
            binaryImg(img_size(1)/2-h/2:img_size(1)/2+h/2, img_size(2)/2+2*gap/2:img_size(2)/2+2*gap/2+w) = 1;

        %% **水平三缝（上下排列）**
        case 'triple_slit_horizontal'
            w = params.height;
            h = params.width;
            gap = params.gap;
            binaryImg(img_size(1)/2-2*gap/2-h:img_size(1)/2-2*gap/2, img_size(2)/2-w/2:img_size(2)/2+w/2) = 1;
            binaryImg(img_size(1)/2-h/2:img_size(1)/2+h/2, img_size(2)/2-w/2:img_size(2)/2+w/2) = 1;
            binaryImg(img_size(1)/2+2*gap/2:img_size(1)/2+2*gap/2+h, img_size(2)/2-w/2:img_size(2)/2+w/2) = 1;
    end
end
