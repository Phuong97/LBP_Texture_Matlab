
%R=3, P=3*8=24, d=1, mapping type = riu2

function H1 = A_RMLBP(I,R,m,d,mapping)
I = imread(I);
I =rgb2gray(I);
I = imadjust(I,stretchlim(I),[]);

[rs , cs] = size(I);

if(R<=2)
    neighbors = R * 8;
else
    neighbors = 24;
end
goc = 360/neighbors;

spoints = zeros(neighbors,2);


bins = mapping.num;

result =[];
m1 = floor(m/2);
Rn=0;

for R1 =R-m1 : R+m1 % c�c v�ng R
% Angle step
Rn= Rn+1;

originx= R1+1;
originy= R1+1;

dx = rs-originx+1;
dy = cs-originy+1;

value1=[];   

for i = originx : 6
    for j =originy : 6
        
    %duy?t ?i?m 
            for k =  1: neighbors
                    if(R1*cosd(k*goc)>0 && R1*sind(k*goc)<0)    
                       spoints(k,1) = ceil(R1*cosd((k)*goc))+i;
                       spoints(k,2) = ceil(-R1*sind((k)*goc))+j;
                    elseif(R1*cosd(k*goc)<0 && R1*sind(k*goc)>0)
                       spoints(k,1) = floor(R1*cosd((k)*goc))+i;
                       spoints(k,2) = floor(-R1*sind((k)*goc))+j;
                    elseif(R1*cosd(k*goc)<0 && R1*sind(k*goc)<0)
                       spoints(k,1) = floor(R1*cosd((k)*goc))+i;
                       spoints(k,2) = ceil(-R1*sind((k)*goc))+j;
                    else
                       spoints(k,1) = ceil(R1*cosd((k)*goc))+i;
                       spoints(k,2) = floor(-R1*sind((k)*goc))+j;
                    end
            end %end neighbors
            
             values = get_pixel(I,spoints);
             value = meanS(R1,m,d,values);
             
       center = double(I(i,j));
       value = double(value);
       value = value(:) - center;
       value1 = [value1 ,sum(value)];
       result = createMatrix(Rn,value1);
    end 
    
end % end 1R
   
        %UPDATE             
                    for j1 =1 : size(result,2)
                        if(result(Rn,j1)~=0)
                          result1(Rn,j1) =  result(Rn,j1);
                        end
                    end
                
end %end nR
result1 = abs(result1);

% Mapping riu2
for i =1: size(result1,1)
    for j =1:size(result1,2)
        result1(i,j) = mapping.table(result1(i,j)+1);
    end
end

% %Hist
H1 = I;
%H1 = hist([result1(1,:)],0:(bins-1));
 %result1(2,:) result1(3,:)
end % END H�M

function L =get_pixel(I , spoints)
   L=[];
    spoints = spoints(:);
    j = size(spoints,1)/2;
    for i = 1 : j
        L1 =  I(spoints(i),spoints(j+i));
        L = [L,L1];
    end
end

function L = meanS(R,m,d,values)
      L=[];      
      values = double(values);
    col = size(values,2);
    for j = 1 : col    
        index = 0;
            for z =  -(m-1)/2 : (m-1)/2
                index = index + ((values(1,j)*(R+z*d))/m);
            end
            L =[L,index];
    end
   
end

function L = createMatrix(R,values)
    for i =1: length(values)
     L(R,i) = values(i);
    end
end