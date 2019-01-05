% ****************************************************************************
% IMAGE VISIBILITY GRAPH 
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
% [1] Visibility graphs of random scalar fields and spatial data
%     Lucas Lacasa and Jacopo Iacovacci
%     Physical Review E 96(1) (2017)
% 
% [2] Visibility graphs for image processing
%     Jacopo Iacovacci and Lucas Lacasa
%     Transactions on Pattern Analysis and Machine Intelligence  (2019)
% 
%
% (c) Jacopo Iacovacci (mriacovacci@hotmail.it)
% *************************************************************************



function Edge_list=imageVisibilityGraph(I,criterion,lattice)

% ****************************************************************************
% This code takes in input a square matrix and extract the corresponding 
% image visibility graph.
%
% Input
% I : input matrix (double/uint8 matrix of size NxN);
%
% criterion : algorithm used to extract patches; choose 'horizontal' or 'natural'.   
%
% lattice : to include in the output graph the lattice (set TRUE) or not (set FALSE). 
%
% Output
% Edge_list : image visibility graph in the form of an edge list.
%*****************************************************************************

n_rows=size(I,1);
n_columns=size(I,2);

if (n_rows~=n_columns)
    error('Input image must be square')
end

if (((~strcmp(criterion,'horizontal'))+(~strcmp(criterion,'natural')))>1)
    error('Criterion string must be <natural> or <horizontal>')
end

if (~islogical(lattice))
    error('lattice must be logical: TRUE to save lattice structure, FALSE otherwise')
end


I=double(I);
edge_list_index=0;

if (strcmp(criterion,'horizontal'))

for i=1:n_rows
    for j=1:n_columns        
        %%%%%%  lattice
        if(lattice) 
        if (j<n_columns)            
            edge_list_index=edge_list_index+1;            
            ei(edge_list_index)=n_columns*(i-1)+j;
            ej(edge_list_index)=n_columns*(i-1)+j+1;
        end
        if(i<n_rows && j<n_columns)
            edge_list_index=edge_list_index+1;
            ei(edge_list_index)=n_columns*(i-1)+j;
            ej(edge_list_index)=n_columns*(i+1-1)+j+1;
        end
        if(i<n_rows)           
            edge_list_index=edge_list_index+1;            
            ei(edge_list_index)=n_columns*(i-1)+j;
            ej(edge_list_index)=n_columns*(i+1-1)+j;
        end
        if(i<n_rows && j>1)
            edge_list_index=edge_list_index+1;            
            ei(edge_list_index)=n_columns*(i-1)+j;
            ej(edge_list_index)=n_columns*(i+1-1)+j-1;
        end               
        end
        %%%%%%
       if(j<n_columns-1)
            k=j+1;
            for c=j+2:n_columns
                cond=1;
                for l=k:c-1
                    if((I(i,l)>=I(i,j))||(I(i,l)>=I(i,c)))
                        cond=0;
                        k=l;
                        break
                    end
                end
                if(I(i,l)>=I(i,j))
                    break
                end
                if(cond==1)
                    edge_list_index=edge_list_index+1;
                    ei(edge_list_index)=n_columns*(i-1)+j;
                    ej(edge_list_index)=n_columns*(i-1)+c;
                end
            end
        end
        
        if(j<n_columns-1 && i<n_rows-1)
            kj=j+1;
            ki=i+1;
            diag_lenght=min(n_rows-i,n_columns-j);
            r=ki;
            for c = j+2:j+diag_lenght
                r=r+1;
                cond=1;
                li=ki-1;
                for lj=kj:c-1
                    li=li+1;
                    if((I(li,lj)>=I(i,j))||(I(li,lj)>=I(r,c)))
                        cond=0;
                        ki=li;
                        kj=lj;
                        break
                    end
                end
                if(I(li,lj)>=I(i,j))
                    break
                end
                if(cond==1)
                    edge_list_index=edge_list_index+1;
                    ei(edge_list_index)=n_columns*(i-1)+j;
                    ej(edge_list_index)=n_columns*(r-1)+c;
                end
            end
        end
        
        if(i<n_rows-1)
            k=i+1;
            for r=i+2:n_rows
                cond=1;
                for l=k:r-1
                    if((I(l,j)>=I(i,j))||(I(l,j)>=I(r,j)))
                        cond=0;
                        k=l;
                        break
                    end
                end
                if(I(l,j)>=I(i,j))
                    break
                end
                
                if(cond==1)
                    edge_list_index=edge_list_index+1;
                    ei(edge_list_index)=n_columns*(i-1)+j;
                    ej(edge_list_index)=n_columns*(r-1)+j;
                end
            end
        end
        
        if(j>2 && i<n_rows-1)
            kj=j-1;
            ki=i+1;
            diag_lenght=min(n_rows-i,j-1);
            c=kj;
            for r = i+2:i+diag_lenght
                c=c-1;
                cond=1;
                lj=kj+1;
                for li=ki:r-1
                    lj=lj-1;
                    if((I(li,lj)>=I(i,j))||(I(li,lj)>=I(r,c)))
                        cond=0;
                        ki=li;
                        kj=lj;
                        break
                    end
                end
                
                if(I(li,lj)>=I(i,j))
                    break
                end
                if(cond==1)
                    edge_list_index=edge_list_index+1;
                    ei(edge_list_index)=n_columns*(i-1)+j;
                    ej(edge_list_index)=n_columns*(r-1)+c;
                end
            end
        end
   end
end

end


if (strcmp(criterion,'natural'))
    
for i=1:n_rows    
    for j=1:n_columns
         
        %%%%%% lattice
        if(~lattice) 
        if (j<n_columns)            
            edge_list_index=edge_list_index+1;            
            ei(edge_list_index)=n_columns*(i-1)+j;
            ej(edge_list_index)=n_columns*(i-1)+j+1;
        end
        if(i<n_rows && j<n_columns)
            edge_list_index=edge_list_index+1;
            ei(edge_list_index)=n_columns*(i-1)+j;
            ej(edge_list_index)=n_columns*(i+1-1)+j+1;
        end
        if(i<n_rows)           
            edge_list_index=edge_list_index+1;            
            ei(edge_list_index)=n_columns*(i-1)+j;
            ej(edge_list_index)=n_columns*(i+1-1)+j;
        end
        if(i<n_rows && j>1)
            edge_list_index=edge_list_index+1;            
            ei(edge_list_index)=n_columns*(i-1)+j;
            ej(edge_list_index)=n_columns*(i+1-1)+j-1;
        end               
        end
        %%%%%%          
        if(j<n_columns-1)
            for c=j+2:n_columns
                cond=1;
                for l=j+1:c-1                    
                    if(I(i,l)>=I(i,j)+(I(i,c)-I(i,j))*(l-j)/(c-j))
                        cond=0;
                        break;      
                    end                   
                end                               
                
                if(cond==1)                    
                    edge_list_index=edge_list_index+1;                    
                    ei(edge_list_index)=n_columns*(i-1)+j;
                    ej(edge_list_index)=n_columns*(i-1)+c;                                      
                end
            end            
        end
        
        if(j<n_columns-1 && i<n_rows-1)            
            diag_lenght=min(n_rows-i,n_columns-j);
            r=i+2;
            for c = j+2:j+diag_lenght                
                cond=1;
                li=i+1;
                for lj=j+1:c-1
                    if(I(li,lj)>=I(i,j)+(I(r,c)-I(i,j))*(lj-j)/(c-j))
                        cond=0;
                        break
                    end
                    li=li+1;                   
                end                
                if(cond==1)
                    edge_list_index=edge_list_index+1;                    
                    ei(edge_list_index)=n_columns*(i-1)+j;
                    ej(edge_list_index)=n_columns*(r-1)+c;
                end                
                 r=r+1;
            end            
        end
        
        if(i<n_rows-1)
            for r=i+2:n_rows
                cond=1;
                for l=i+1:r-1
                    if(I(l,j)>=I(i,j)+(I(r,j)-I(i,j))*(l-i)/(r-i))
                        cond=0;
                        break
                    end                    
                end                                                
                if(cond==1)                    
                    edge_list_index=edge_list_index+1;
                    ei(edge_list_index)=n_columns*(i-1)+j;
                    ej(edge_list_index)=n_columns*(r-1)+j;
                end                
            end            
        end
        
        if(j>2 && i<n_rows-1)
            diag_lenght=min(n_rows-i,j-1);
            c=j-2;
            for r = i+2:i+diag_lenght
                cond=1;
                lj=j-1;
                for li=i+1:r-1
                    if(I(li,lj)>=I(i,j)+(I(r,c)-I(i,j))*(lj-j)/(c-j))    
                        cond=0;
                        break
                    end
                    lj=lj-1;
                end
                if(cond==1)
                    edge_list_index=edge_list_index+1;
                    ei(edge_list_index)=n_columns*(i-1)+j;
                    ej(edge_list_index)=n_columns*(r-1)+c;
                end
                 c=c-1;  
            end
        end
    end
end
end


Edge_list=[ei',ej'];


end


    
    