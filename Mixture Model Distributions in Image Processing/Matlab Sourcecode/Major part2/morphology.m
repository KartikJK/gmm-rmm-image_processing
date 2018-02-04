
img=imread('4.jpg');
n1=img;
x=im2bw(n1);
figure(2),imshow(x);
se = strel('disk',3);
BW= imerode(x,se);
figure(3),imshow(BW);
h=im2double(n1);
i=im2double(BW);
v=imadd(h,i);
im2(:,:,1)=v;
im2(:,:,2)=v;
im2(:,:,3)=v;
for i=1:256
for j=1:256
if v(i,j)>1.5:1.8
%im2(i,j,1)=255;im2(i,j,3)=0;im2(i,j,2)=0;else im2(i,j,1)=v(i,j);
im2(i,j,2)=255;im2(i,j,3)=0;im2(i,j,1)=0;else im2(i,j,2)=v(i,j);
%im2(i,j,3)=255;im2(i,j,2)=0;im2(i,j,1)=0;else im2(i,j,3)=v(i,j);
%im2(i,j,3)=255;im2(i,j,2)=0;im2(i,j,1)=255;
%im2(i,j,2)=255;im2(i,j,3)=0;im2(i,j,1)=255;
end
end
end
figure(4),imshow(im2);
Img=BW;
%b=abs(BW);
r=mean(mean(BW));
volume=r