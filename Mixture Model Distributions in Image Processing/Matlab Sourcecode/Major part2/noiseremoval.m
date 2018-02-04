img=imread('4.jpg');
G = fspecial('gaussian',[5 5],2);
Ig = imfilter(img,G,'same');
figure(1),imshow(Ig);