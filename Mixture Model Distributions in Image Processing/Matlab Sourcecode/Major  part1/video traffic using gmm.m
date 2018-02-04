clear all

v = VideoReader('C:\Users\kartik\Desktop\major\visiontraffic.avi');
% -----------------------  frame size variables -----------------------
%fr = step(v); % read the next video frame
%fr = v(1).cdata; % read in 1st frame as background frame
fr = read(v,1);
fr_bw = rgb2gray(fr);     % convert background to greyscale
fr_size = size(fr);             
width = fr_size(2);
height = fr_size(1);
fg = zeros(height, width);
bg_bw = zeros(height, width);

% --------------------- gmm variables -----------------------------------

C = 3;                                  % number of gaussian components (typically 3-5)
M = 3;                                  % number of background components
D = 2.5;                                % positive deviation threshold
alpha = 0.005;                           % learning rate (between 0 and 1) (from paper 0.01)
thresh = 0.25;                          % foreground threshold (0.25 or 0.75 in paper)
sd_init = 30/255;                            % initial standard deviation (for new components) var = 36 in paper
w = zeros(height,width,C);              % initialize weights array
mean = zeros(height,width,C);           % pixel means
sd = zeros(height,width,C);             % pixel standard deviations
u_diff = zeros(height,width,C);         % difference of each pixel from mean
p = alpha/(1/C);                        % initial p variable (used to update mean and sd)
rank = zeros(1,C);                      % rank of components (w/sd)


% --------------------- initialize component means and weights -----------

pixel_depth = 8;                        % 8-bit resolution
pixel_range = 2^pixel_depth -1;         % pixel range (# of possible values)

for i=1:height
    for j=1:width
        for k=1:C
            
            mean(i,j,k) = rand*pixel_range;     % means random (0-255)
            w(i,j,k) = 1/C;                     % weights uniformly dist
            sd(i,j,k) = sd_init;                % initialize to sd_init
            
        end
    end
end

%--------------------- process frames -----------------------------------

for n = 1:150
    %fr = step(v);
    %fr = v(n).cdata;       % read in frame
    %fr = read(v,[1 Inf]);
    fr = read(v,n);
    fr_bw = rgb2gray(fr);       % convert frame to grayscale
    figure(1),subplot(3,1,1),imshow(fr)
    % calculate difference of pixel values from mean
    for m=1:C
        u_diff(:,:,m) = abs(double(fr_bw) - double(mean(:,:,m)));
    end
     
    % update gaussian components for each pixel
    for i=1:height
        for j=1:width
            
            match = 0;
            for k=1:C                       
                if (abs(u_diff(i,j,k)) <= D*sd(i,j,k))       % pixel matches component
                    
                    match = 1;                          % variable to signal component match
                    
                    % update weights, mean, sd, p
                    w(i,j,k) = (1-alpha)*w(i,j,k) + alpha;
                    p = alpha/w(i,j,k);                  
                    mean(i,j,k) = (1-p)*mean(i,j,k) + p*double(fr_bw(i,j));
                    sd(i,j,k) =   sqrt((1-p)*(sd(i,j,k)^2) + p*((double(fr_bw(i,j)) - mean(i,j,k)))^2);
                else                                    % pixel doesn't match component
                    w(i,j,k) = (1-alpha)*w(i,j,k);      % weight slighly decreases
                    
                end
            end
            
            w(i,j,:) = w(i,j,:)./sum(w(i,j,:));
            
            bg_bw(i,j)=0;
            for k=1:C
                bg_bw(i,j) = bg_bw(i,j)+ mean(i,j,k)*w(i,j,k);
            end
            
            % if no components match, create new component
            if (match == 0)
                [min_w, min_w_index] = min(w(i,j,:));  
                mean(i,j,min_w_index) = double(fr_bw(i,j));
                sd(i,j,min_w_index) = sd_init;
            end

            rank = w(i,j,:)./sd(i,j,:);             % calculate component rank
            rank_ind = [1:1:C];
            
            % sort rank values
            for k=2:C               
                for m=1:(k-1)
                    
                    if (rank(:,:,k) > rank(:,:,m))                     
                        % swap max values
                        rank_temp = rank(:,:,m);  
                        rank(:,:,m) = rank(:,:,k);
                        rank(:,:,k) = rank_temp;
                        
                        % swap max index values
                        rank_ind_temp = rank_ind(m);  
                        rank_ind(m) = rank_ind(k);
                        rank_ind(k) = rank_ind_temp;    

                    end
                end
            end
            
            % calculate foreground
            match = 0;
            k=1;
            
            fg(i,j) = 0;
            while ((match == 0)&&(k<=M))

                if (w(i,j,rank_ind(k)) >= thresh)
                    if (abs(u_diff(i,j,rank_ind(k))) <= D*sd(i,j,rank_ind(k)))
                        fg(i,j) = 0;
                        match = 1;
                    else
                        fg(i,j) = fr_bw(i,j);     
                    end
                end
                k = k+1;
            end
        end
    end
    
    %figure(1),subplot(3,1,1),imshow(fr)
    subplot(3,1,2),imshow(uint8(bg_bw))
    subplot(3,1,3),imshow(uint8(fg)) 

    se = strel('square', 3);
filteredForeground = imopen(uint8(fg), se);
figure(2); imshow(filteredForeground); title('Clean Foreground');

blobAnalysis = vision.BlobAnalysis('BoundingBoxOutputPort', true, ...
    'AreaOutputPort', false, 'CentroidOutputPort', false, ...
    'MinimumBlobArea', 150);
bbox = step(blobAnalysis, logical(filteredForeground));

result = insertShape(fr, 'Rectangle', bbox, 'Color', 'green');

numCars = size(bbox, 1);
result = insertText(result, [10 10], numCars, 'BoxOpacity', 1, ...
    'FontSize', 14);
figure (3); imshow(result); title('Detected Cars');
    
    %Mov1(n)  = im2frame(uint8(fg),gray);           % put frames into movie
    %Mov2(n)  = im2frame(uint8(bg_bw),gray);           % put frames into movie
    
end
      
%movie2avi(Mov1,'mixture_of_gaussians_output','fps',3);           % save movie as avi 
%movie2avi(Mov2,'mixture_of_gaussians_background','fps',30);           % save movie as avi 

 
