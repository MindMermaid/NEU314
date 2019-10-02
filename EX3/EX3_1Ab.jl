# Orientation
cd("NEU314\\EX3")
include("standard_start.jl");

# Preparation
img = imread("el-capitan.png");
imshow(img);
size(img);

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

    return R, G, B
end
