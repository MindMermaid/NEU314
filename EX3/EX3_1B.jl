##=====Circular Rotation Function

"""
circulate_pixels_R(arg1, arg2) -- circulates red channel of an image
Reads an image file, shifts specified rows (p) of red channel in image and returns shifted image

Args:
   arg1 (str): filename
   arge2 (int): p -- number of pixels by which to shift p rows of the red channel

Returns:
   out (image): returns shifted image
"""

function circulate_pixels_R(filename, p)
    
       img = imread(filename);
    
    img_top = img[1:p+1,:,1];
    img_bottom = img[p+1:end,:,1];
    
    img[1:end-p,:,1] = img_bottom;
    img[end-p:end,:,1] = img_top;  
    
    return imshow(img)
end

img = imread("el-capitan.png")


fig = figure(figsize=(10,4))
suptitle("Circulate R-channel");
    subplot(121)
        imshow(img)
        title("Original")
        axis("off")
    subplot(122)
        circulate_pixels_R("el-capitan.png",180)
        title("Shifted R-channel")
        axis("off")
