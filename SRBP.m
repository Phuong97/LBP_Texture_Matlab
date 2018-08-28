% thuc hien tren anh gray
%SRBP thuc hien quet theo ban kinh R = ban kinh
% tra ve 5 array chua cac gia tri ham mean,med,var,min,max
% gop 5 array = 1 hist tra ve

%Matrix luu cac gia tri diem anh qua moi R = bán kính
%MatrixBinary luu gia tri nhi phan 1 0 1:duong 0:am

%gom 3 buoc duoc ap dung trong 1 o voi R = ban kính
% sau khi duyet cac phan tu cua o ta thuc hien cac buoc 

%row 59 : buoc 1: STATISTICAL CHARACTERIZATION
% row 85 : buoc 2: RADIAN DIFFERENCE
%row : buoc 3: STATISTICAL RADIAL BINARY PATTERNS

%function mapping = SRBP(R)
%function [Smean, Smed,Svar,Smin,Smax] = SRBP(R)
function hist1 = SRBP(image,R)
image = imread(image);
image = rgb2gray(image);  %chuyen thanh gray
image = imadjust(image,stretchlim(image),[]);
[r c] = size(image); 
sMean=[]; 
sMed=[];
sVar = [];
sMin=[];
sMax=[];

for i= R+1 : r-R
    for j=R+1 : c-R
     
        for k =1:R %duyet cac vong tu 1 den R
           N = k*8; % so diem tren 1 duong chon (1,8) (2,16 ) (3,24)...
           goc = 360 / N; % so do goc
        
           values = zeros(); % 
                for z= N :-1: 1 %chay tung goc
                    goc1 = goc*z;
                    
                    if(k*cosd(goc1)>0 && k*sind(goc1)<0)    
                        x = ceil(k*cosd(goc1));
                        y = ceil(-k*sind(goc1));
                    elseif(k*cosd(goc1)<0 && k*sind(goc1)>0)
                        x = floor(k*cosd(goc1));
                        y = floor(-k*sind(goc1));
                    elseif(k*cosd(goc1)<0 && k*sind(goc1)<0)
                         x = floor(k*cosd(goc1));
                         y = ceil(-k*sind(goc1));
                    else
                         x = ceil(k*cosd(goc1));
                         y = floor(-k*sind(goc1));
                    end
                     P = get_pixel(image,i+x,j+y);
                     values = [values,P];    % cho vao 1 array
                 end
           
             values = double(values);
             
             L = createMatrix(k,values);
             
             
             %update matrix sau moi lan chay xong 1 vong 
            [r1 c1] = size(L);       
                for i1 = 1 : r1
                    for j1 =1 : c1
                        if(L(i1,j1)~=0)
                          Matrix1(i1,j1) =  L(i1,j1);
                        end
                    end
                
                end
          
            %
        end % end R
       
         %Tinh Mean Max ...
        %chuyen cac gia tri ve do dai chuan la 8 theo
        %BRINT(8=8;16->8,24->8...)
    
       
        [r2 c2] = size(Matrix1);
       
             
        for ii = 1 : r2
             index = 0;
             count = 1;
             value1=[];
            for jj = 1 : c2  %Tinh BRINT 
                if(count <= ii && Matrix1(ii,jj)~=0)
                    index = index + Matrix1(ii,jj) + count-1;                       
                    count=count+1;
                    if(count > ii)
                        index = index ./ii;
                        value1= [value1,index];
                        count =1;
                        index = 0;                     
                    end   
                        
                end             
              
            end
                 value1 = double(value1);
                 L1 = createMatrix1(ii,value1);
                 
                 [r3 c3] = size(L1);       
                 for ii1 = 1 : r3
                    for jj1 =1 : c3
                        if(L1(ii1,jj1)~=0)
                          Matrix(ii1,jj1) =  L1(ii1,jj1);
                        end
                    end
                
                  end
              
        end %end row ii
      
     
        
         %MatrixBinary = zeros(5,R);
         % tinh cac gia tri mean,med... cua cac R = 1,2,3..
           
           
%STATISTICAL CHARACTERIZATION
          [r5 c5] = size(Matrix);
          for i2 = 1 : r5
            mean1 = mean(Matrix(i2,:));
            MatrixBinary(1,i2) = mean1;

            
            med1 = median(Matrix(i2,:));
            MatrixBinary(2,i2) = med1;
                       
            var1= var(Matrix(i2,:));
            MatrixBinary(3,i2) = var1;
            
            min1 = min(Matrix(i2,:));
            MatrixBinary(4,i2) = min1;      
          
            max1 =max(Matrix(i2,:));
            MatrixBinary(5,i2) = max1;      
          end
           
        
            % Het 1 o bat dau tinh toan ra gia tri tai R o do ma chuyen sang dong ke tiep
            % chuyen ve nhi phan va nhan de tinh ra gia tri
         

%RADIAN DIFFERENCE
 [r4 c4] = size(MatrixBinary);
          for i5=1:r4
              for j5=c4:-1:1
              if(j5==1)
                  MatrixBinary(i5,j5) = MatrixBinary(i5,j5)- double(get_pixel(image,R+1,R+1));
              else
                  MatrixBinary(i5,j5) = MatrixBinary(i5,j5) - MatrixBinary(i5,j5-1);
              end
              end
          end
           
%STATISTICAL RADIAL BINARY PATTERNS
          MatrixBinary = binaryMatrix(MatrixBinary);
          LL = multiMatrix(MatrixBinary);
          sMean = [sMean,LL(1)];
          sMed  = [sMed,LL(2)];
          sVar  = [sVar,LL(3)];
          sMin  = [sMin,LL(4)];
          sMax  = [sMax,LL(5)];
     
    end %end column
    hist1 = histogram([sMean sMed sVar sMin sMax]);
    
end % end Row      

    
end % end hàm
 
% Lay gia tri tai vi tri x,y trong 1 anh
function L =get_pixel(image,idx,idy)
try
    L = image(idx,idy);
catch
    L = 0;
end
end


% tao matran moi theo 1 dong so 0 o dau
function L = createMatrix(R,values)
    for i =2: length(values)
     L(R,i-1) = values(i);
    end
end

%k co so 0 o dau
function L = createMatrix1(R,values)
    for i =1: length(values)
     L(R,i) = values(i);
    end
end

% tao matran moi theo 1 cot
function L = createMatrixColumn(R,values)
    for i = 1: length(values)
     L(i,R-1) = values(i);
    end
end

% chuyen ve 1 va 0
function L = binaryMatrix(matrix)
[r c] = size(matrix);
    for i = 1 :r 
        for j =1 :c
            if(matrix(i,j)>0)
                matrix(i,j)=1;
            else
                matrix(i,j)=0;
            end
        end
    end
    L = matrix;
end

%nhan 2 matran
function L = multiMatrix(matrix)
    [r c] = size(matrix);
    for i = 0 : c-1
    matrix1(i+1,1) = 2.^i;
    end
   L =matrix*matrix1; 
end
