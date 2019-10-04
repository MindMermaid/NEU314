##===== Preparation
img = imread("eye.png");
figure(1)
imshow(img);
##=====

##=====Circular Rotation Function

"""
pixel_circulate(arg1, arg2) -- circulates specified number of rows of pixels in an image
Reads an image file, shifts specified rows of an image and returns shifted image

Args:
   arg1 (str): filename
   arge2 (int): p -- number of pixels by which to shift p rows

Returns:
   out (image): returns 3 arrays containing the R, G & B values of every image pixel
"""

function pixel_circulate(filename, p)
