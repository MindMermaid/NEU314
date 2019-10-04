# Preparation
img = imread("el-capitan.png");
figure(1)
    imshow(img);

"""
image_RGB(arg1) -- outputs RGB values of pixels in image file
Reads an image file, extracts RGB values from the data of each pixel and
stores each value of R, G & B inside a separate array

Args:
   arg1 (str): filename

Returns:
   out (int): returns 3 arrays containing the R, G & B values of every image pixel
"""

function image_RGB(filename)

    img = imread(filename);
    R = img[:,:,1];
    G = img[:,:,2];
    B = img[:,:,3];

     return [R, G, B]
end

##===== Swap RGB values to create new images!
image1 = image_RGB("el-capitan.png");
#BGR
image2 = cat(image1[3], image1[2], image1[1]; dims=3);
#BRG
image3 = cat(image1[3], image1[1], image1[2]; dims=3);
#GBR
image4 = cat(image1[2], image1[3], image1[1]; dims=3);
#GRB
image5 = cat(image1[2], image1[1], image1[1]; dims=3);
#RBG
image6 = cat(image1[1], image1[3], image1[2]; dims=3);
##===== Create a figure to display swapped RGB values
fig = figure(figsize=(10,4))
    suptitle("RGB Value Swaps")
        subplot(231)
            imshow(img);
            title("RGB (original)")
            axis("off")
        subplot(232)
            imshow(image2);
            title("BGR")
            axis("off")
        subplot(233)
            imshow(image3);
            title("BRG")
            axis("off")
        subplot(234)
            imshow(image4);
            title("GBR")
            axis("off")
        subplot(235)
            imshow(image5);
            title("GRB")
            axis("off")
        subplot(236)
            imshow(image6);
            title("RBG")
            axis("off")
##=====
