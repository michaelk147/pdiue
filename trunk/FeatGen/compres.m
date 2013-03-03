img = imread('/home/michael/pedestrian_datasets/INRIAPerson/96X160H96/Train/pos/crop001001a.png');
%img = imread('/home/michael/pedestrian_datasets/caltech/data-INRIA/learning/train_neg_hard/set02V000_00000329a_hard_1.png');
stride = 8;
nw = 64;
nh = 128;
oBin = 9;
cells = 4;

% skip extra space
hspace = (size(img,1)-nh)/2 /stride;
wspace = (size(img,2)-nw)/2 /stride;

hogspww = nw/stride-2;
hogspwh = nh/stride-2;

fullmF = hog(double(img));
figure(1);
im(img);
fullmFV = hogDraw(fullmF);

figure(2);
im(fullmFV);

mF = fullmF(1+hspace:1+hspace+hogspwh - 1,1+wspace:1+wspace+hogspww-1,:);
mFV = hogDraw(mF);
figure(3);
im(mFV);


FR = reshape(F,[hogspwh hogspww cells*oBin]);
FRV = hogDraw(FR);
figure(4);
im(FRV);