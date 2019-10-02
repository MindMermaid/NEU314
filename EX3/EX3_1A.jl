# Orientation
cd("NEU314\\EX3");
include("standard_start.jl");

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

function image_RGB(filename::AbstractString)

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

##===== Create a figure to display swapped RGB values
fig = figure("Flip RGB",figsize=(10,4))
    suptitle("RGB Value Swaps")
        subplot(121)
            imshow(img);
            title("RGB (original)")
            axis("off")
        subplot(122)
            imshow(image2);
            title("BGR")
            axis("off")
##=====
