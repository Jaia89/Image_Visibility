% ****************************************************************************
% VISIBILITY PATCHES 
% 
% This code can be redistributed and/or modified under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
%  
% This program is distributed by the author in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%  
% 
% Please cite:
% 
% [1] Visibility graphs for image processing
%     Jacopo Iacovacci and Lucas Lacasa
%     Transactions on Pattern Analysis and Machine Intelligence  (2019)
% 
%
% (c) Jacopo Iacovacci (mriacovacci@hotmail.it)
% *************************************************************************




function Z=visibilityPatches(I,stride,criterion)

% ****************************************************************************
% This code takes in input an image and extract from each channel the frequency 
% distribution associated to visibility patches of size 3x3 pixels. 
%
% Input
% I : matrix of the input image; for bw images I is a double/uint8 matrix of
% size NxN; for rgb images I is a double/uint8 matrix of size NxNx3.
%     
% stride : integer value representing the distance in unit of pixels at which patches 
%          are extracted; stride=1 extracts a patch around each pixel of the image 
%          (except borders).   
%
% criterion : algorithm used to extract patches; choose 'horizontal' or 'natural'.    
% 
% Output
% Z : matrix of size Mx256 containing the frequency values of the 256
% visibility patches in each of the M image channels.
%*****************************************************************************


motif_size=3;
pow2=[1,2,4,8,16,32,64,128];

n_rows=size(I,1);
n_cols=size(I,2);

stride=round(stride);

if (n_rows~=n_cols)
    error('Input image must be square')
end

if (stride>=n_rows-motif_size+1 || stride>=n_cols-motif_size+1)
    error('Stride length exceeds image size')
end

if (((~strcmp(criterion,'horizontal'))+(~strcmp(criterion,'natural')))>1)
    error('Criterion string must be <natural> or <horizontal>')
end


Z=[];
I=double(I);


if (strcmp(criterion,'horizontal'))
    
    for m=1:size(I,3)
        
        Freq=zeros(1,256);
        count=0;
        
        for i=1:stride:n_rows-motif_size+1
            for j=1:stride:n_cols-motif_size+1
                
                count=count+1;
                string=zeros(1,8);
                M=I(i:i+motif_size-1,j:j+motif_size-1,m);
                
                if(M(1,1)>M(1,2)&&M(1,3)>M(1,2))
                    string(1)=1;
                end
                
                if(M(1,3)>M(2,3)&&M(3,3)>M(2,3))
                    string(2)=1;
                end
                
                if(M(3,1)>M(3,2)&&M(3,3)>M(3,2))
                    string(3)=1;
                end
                
                if(M(1,1)>M(2,1)&&M(3,1)>M(2,1))
                    string(4)=1;
                end
                
                if(M(1,2)>M(2,2)&&M(3,2)>M(2,2))
                    string(5)=1;
                end
                
                if(M(2,1)>M(2,2)&&M(2,3)>M(2,2))
                    string(6)=1;
                end
                
                if(M(1,1)>M(2,2)&&M(3,3)>M(2,2))
                    string(7)=1;
                end
                
                if(M(1,3)>M(2,2)&&M(3,1)>M(2,2))
                    string(8)=1;
                end
                
                label=sum(string.*pow2)+1;
                Freq(label)=Freq(label)+1;
                
            end
        end
        Freq= Freq./count;
        Z=[Z;Freq];
        
    end
end

if (strcmp(criterion,'natural'))
    
    for m=1:size(I,3)
        
        Freq=zeros(1,256);
        count=0;
        
        for i=1:stride:n_rows-motif_size+1
            for j=1:stride:n_cols-motif_size+1
                
                count=count+1;
                string=zeros(1,8);
                M=I(i:i+motif_size-1,j:j+motif_size-1,m);
                
                if(M(1,3)>2*M(1,2)-M(1,1))
                    string(1)=1;
                end
                
                if(M(3,3)>2*M(2,3)-M(1,3))
                    string(2)=1;
                end
                
                if(M(3,1)>2*M(3,2)-M(3,3))
                    string(3)=1;
                end
                
                if(M(1,1)>2*M(2,1)-M(3,1))
                    string(4)=1;
                end
                
                if(M(3,2)>2*M(2,2)-M(1,2))
                    string(5)=1;
                end
                
                if(M(2,3)>2*M(2,2)-M(2,1))
                    string(6)=1;
                end
                
                if(M(3,3)>2*M(2,2)-M(1,1))
                    string(7)=1;
                end
                
                if(M(3,1)>2*M(2,2)-M(1,3))
                    string(8)=1;
                end
                
                label=sum(string.*pow2)+1;
                Freq(label)=Freq(label)+1;
                
            end
        end
        
        Freq= Freq./count;
        Z=[Z;Freq];
        
    end
end


end
