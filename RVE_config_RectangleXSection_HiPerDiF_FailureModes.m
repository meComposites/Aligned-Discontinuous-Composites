% function [RVE,nc] = RVE_config_RectangleXSection_HiPerDiF (n_rows,n_columns,nc,Vc,x_section_arrangement)
function [RVE,nc] = RVE_config_RectangleXSection_HiPerDiF_FailureModes(RVE_loop_counter,n_rows,n_columns,nc,Vc,x_section_arrangement,c_width,n_columns_c,HiPerDiF_fuzzyness,Fuzzy_Factor,rng_switch)
% function [RVE,nc] = RVE_config_RectangleXSection ()
%% RVE_config
% FJ - generates an RVE cross-section configuration with a specified
% pattern of glass and carbon fibres
% n_rows=16;
% n_columns=372;
% nc=1300;
% Vc=0.4;
% x_section_arrangement=5;
switch rng_switch
    case 0
        rng(RVE_loop_counter);
    case 1
        rng('shuffle','simdTwister')       
    otherwise
        rng(rng_switch);
end

HiPerDiF_arrangement=4;
HiPerDiF_rows_resolution=floor(n_rows/4);
HiPerDiF_columns_resolution=floor(n_columns/9);
HiPerDiF_Vc_A=1.0;
HiPerDiF_nc_A=floor(HiPerDiF_Vc_A*HiPerDiF_rows_resolution*HiPerDiF_columns_resolution);
HiPerDiF_Vc_B=0;
HiPerDiF_nc_B=floor(HiPerDiF_Vc_B*HiPerDiF_rows_resolution*HiPerDiF_columns_resolution);

RVE=zeros(n_rows,n_columns);


if x_section_arrangement==1
    
    % checkerboard pattern
    % ### FJ - na modified to support n_rows x n_columns rectangle RVE ###
    x_check_c=min(floor(sqrt(nc/5)), ceil(min(n_rows,n_columns)/3));
    
    row1=[ones(x_check_c),zeros(x_check_c,n_columns-2*x_check_c),ones(x_check_c)];
    row2_height=floor((n_rows-3*x_check_c)/2);
    if row2_height <0
        row2_height=0;
    end
    row2=zeros(row2_height,n_columns);
    row3_Lzeros_width=ceil((n_columns-x_check_c)/2);
    row3_Rzeros_width=floor((n_columns-x_check_c)/2);
    row3=[zeros(x_check_c,row3_Lzeros_width),ones(x_check_c),zeros(x_check_c,row3_Rzeros_width)];
    row4_height=ceil((n_rows-3*x_check_c)/2);
    if row4_height <0
        row4_height=0;
    end
    row4=zeros(row4_height,n_columns);
    row5=[ones(x_check_c),zeros(x_check_c,(n_columns-2*x_check_c)),ones(x_check_c)];
    
    RVE=[row1;row2;row3;row4;row5];
    
    rem_c=nc-sum(sum(RVE));
    check_1_counter=0;
    check_2_counter=0;
    check_3_counter=0;
    check_4_counter=0;
    check_5_counter=0;
    check_6_counter=0;
    
    while rem_c>=1
        ident_check=mod(rem_c,6);
        
        if ident_check==1
            check_1_counter=check_1_counter+1;
            column_counter=ceil(check_1_counter/x_check_c);
            row_counter=mod(check_1_counter,x_check_c)+1;
            row_locn=row_counter;
            column_locn=x_check_c+column_counter;
            
            RVE(int8(row_locn),int8(column_locn))=1;
            rem_c=nc-sum(sum(RVE));
        end
        
        if ident_check==2
            check_2_counter=check_2_counter+1;
            column_counter=ceil(check_2_counter/x_check_c);
            row_counter=mod(check_2_counter,x_check_c)+1;
            row_locn=row_counter;
            column_locn=n_columns-(x_check_c-1)-column_counter;
            
            RVE(int8(row_locn),int8(column_locn))=1;
            rem_c=nc-sum(sum(RVE));
        end
        
        if ident_check==3
            check_3_counter=check_3_counter+1;
            column_counter=ceil(check_3_counter/x_check_c);
            row_counter=mod(check_3_counter,x_check_c)+1;
            row_locn=x_check_c+row2_height+row_counter;
            column_locn=row3_Lzeros_width+1-column_counter;
            
            RVE(int8(row_locn),int8(column_locn))=1;
            rem_c=nc-sum(sum(RVE));
        end
        
        if ident_check==4
            check_4_counter=check_4_counter+1;
            column_counter=ceil(check_4_counter/x_check_c);
            row_counter=mod(check_4_counter,x_check_c)+1;
            row_locn=x_check_c+row2_height+row_counter;
            column_locn=n_columns-row3_Rzeros_width+column_counter;
            
            RVE(int8(row_locn),int8(column_locn))=1;
            rem_c=nc-sum(sum(RVE));
        end
        
        if ident_check==5
            check_5_counter=check_5_counter+1;
            column_counter=ceil(check_5_counter/x_check_c);
            row_counter=mod(check_5_counter,x_check_c)+1;
            row_locn=n_rows-x_check_c+row_counter;
            column_locn=x_check_c+column_counter;
            
            RVE(int8(row_locn),int8(column_locn))=1;
            rem_c=nc-sum(sum(RVE));
        end
        
        if ident_check==0
            check_6_counter=check_6_counter+1;
            column_counter=ceil(check_6_counter/x_check_c);
            row_counter=mod(check_6_counter,x_check_c)+1;
            row_locn=n_rows-x_check_c+row_counter;
            column_locn=n_columns-x_check_c+1-column_counter;
            
            RVE(int8(row_locn),int8(column_locn))=1;
            rem_c=nc-sum(sum(RVE));
        end
    end
    
elseif x_section_arrangement==2
    % intermingled sandwich pattern
    fraction_intermingled=0.3;
    RVE_sandwich_n_rows=floor((1-fraction_intermingled)*nc/(n_columns*2));
    Vf=0.35;
    dc=0.010;
    dg=0.007;
    
    % FJ - generated to calculate dimensions of continuum elements
    Vc_intermingled=(1-(2*RVE_sandwich_n_rows*n_columns/nc))*Vc;
    
    nc_sandwich=RVE_sandwich_n_rows*n_columns;
    % FJ - calculate the actual dimensions of each of the zones,keeping Vf
    % constant and accounting for he differences in fibre sizes. Note that
    % the width dimension must remain constant despite the differences in
    % Vf and fibre diameters, therefore the intermingled height is adjusted
    % to facilitate this.
    RVE_sandwich_area=(nc_sandwich*pi()*(dc/2)^2)/Vf;
    RVE_width=sqrt((n_columns/RVE_sandwich_n_rows)*RVE_sandwich_area);
    RVE_sandwich_height=RVE_sandwich_area/RVE_width;
    RVE_sandwich_area_check=RVE_width*RVE_sandwich_height;
    RVE_intermingled_area=(((nc-2*RVE_sandwich_n_rows*n_columns)*pi()*(dc/2)^2)+((n_rows*n_columns-nc)*pi()*(dg/2)^2))/Vf;
    RVE_intermingled_height=RVE_intermingled_area/RVE_width;
    RVE_height=RVE_intermingled_height+2*RVE_sandwich_height;
    RVE_area=(nc*pi()*(dc/2)^2+(n_rows*n_columns-nc)*pi()*(dg/2)^2)/Vf;
    RVE_area_check=RVE_width*RVE_height;
    
    RVE_sandwich=ones(RVE_sandwich_n_rows,n_columns);
    
    RVE=[RVE_sandwich;zeros((n_rows-2*RVE_sandwich_n_rows),n_columns);RVE_sandwich];
    
    while sum(sum(RVE))~=nc
        RVE((RVE_sandwich_n_rows+ceil((n_rows-2*RVE_sandwich_n_rows)*rand())),ceil(n_columns*rand()))=1;
    end
    
elseif x_section_arrangement==3
    % FJ - Create a matrix with a random placement of glass and carbon fibres
    % FJ - start by defining a matrix with all glass fibres
    RVE=zeros(n_rows,n_columns); %  Glass = 0  -  Carbon = 1
    % FJ - replace a fibre with a carbon fibre until the number of carbon
    % fibres matches nc
    while sum(sum(RVE))~=nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    % FJ - still not sure what this means
    % RVE=RVEg;
    
elseif x_section_arrangement==4
    % FJ - create an RVE with one set of fibres only
    
    if Vc==0
        RVE=zeros(n_rows,n_columns);
    elseif Vc==1
        RVE=ones(n_rows,n_columns);
    else
        error('Vc can only be 1 or zero for this RVE x-section arrangement');
    end
    
elseif x_section_arrangement==5
    % make your own HiPerDiF arrangement
    % Assign HiPerDiF arrangement, based on limited resolution of HiPerDiF
    % cells
    
    % "Dispersed - 1": diagonal checkerboard
    if HiPerDiF_arrangement==1
        HiPerDiF_RVE=['A' 'B' 'B' 'A' 'B' 'B' 'A' 'B' 'B';'B' 'B' 'A' 'B' 'B' 'A' 'B' 'B' 'A';'B' 'A' 'B' 'B' 'A' 'B' 'B' 'A' 'B';'A' 'B' 'B' 'A' 'B' 'B' 'A' 'B' 'B'];
    end
    % "Dispersed - 2": vertical stripes
    if HiPerDiF_arrangement==2
        HiPerDiF_RVE=['B' 'A' 'B' 'B' 'A' 'B' 'B' 'A' 'B';'B' 'A' 'B' 'B' 'A' 'B' 'B' 'A' 'B';'B' 'A' 'B' 'B' 'A' 'B' 'B' 'A' 'B';'B' 'A' 'B' 'B' 'A' 'B' 'B' 'A' 'B'];
    end
    % "Dispersed - 3": enclosed boxes
    if HiPerDiF_arrangement==3
        HiPerDiF_RVE=['B' 'B' 'B' 'B' 'B' 'B' 'B' 'B' 'B';'B' 'A' 'A' 'A' 'B' 'A' 'A' 'A' 'B';'B' 'A' 'A' 'A' 'B' 'A' 'A' 'A' 'B';'B' 'B' 'B' 'B' 'B' 'B' 'B' 'B' 'B'];
    end
    % "Dispersed - 4": blocked
    if HiPerDiF_arrangement==4
        HiPerDiF_RVE=['B' 'B' 'B' 'A' 'A' 'A' 'B' 'B' 'B';'B' 'B' 'B' 'A' 'A' 'A' 'B' 'B' 'B';'B' 'B' 'B' 'A' 'A' 'A' 'B' 'B' 'B';'B' 'B' 'B' 'A' 'A' 'A' 'B' 'B' 'B'];
    end
    % Count number of HiPerDiF cells in arangemnts
    [rows_HiPerDiF_cells,columns_HiPerDiF_cells]=size(HiPerDiF_RVE);
    
    RVE=zeros((rows_HiPerDiF_cells*HiPerDiF_rows_resolution),(columns_HiPerDiF_cells*HiPerDiF_columns_resolution));
    
    for i=1:rows_HiPerDiF_cells
        for j=1:columns_HiPerDiF_cells
            if HiPerDiF_RVE(i,j)=='A'
                HiPerDiF_count_nc_A=0;
                while HiPerDiF_count_nc_A~=HiPerDiF_nc_A
                    RVE_row_counter=((i-1)*HiPerDiF_rows_resolution+1)+(floor(HiPerDiF_rows_resolution*rand()));
                    RVE_column_counter=(((j-1)*HiPerDiF_columns_resolution+1)+(floor(HiPerDiF_columns_resolution*rand())));
                    RVE(RVE_row_counter,RVE_column_counter)=1;
                    HiPerDiF_count_nc_A=sum(sum(RVE(((i-1)*HiPerDiF_rows_resolution+1):(i*HiPerDiF_rows_resolution),((j-1)*HiPerDiF_columns_resolution+1):(j*HiPerDiF_columns_resolution))));
                end
            elseif HiPerDiF_RVE(i,j)=='B'
                HiPerDiF_count_nc_B=0;
                while HiPerDiF_count_nc_B~=HiPerDiF_nc_B
                    RVE_row_counter=((i-1)*HiPerDiF_rows_resolution+1)+(floor(HiPerDiF_rows_resolution*rand()));
                    RVE_column_counter=(((j-1)*HiPerDiF_columns_resolution+1)+(floor(HiPerDiF_columns_resolution*rand())));
                    RVE(RVE_row_counter,RVE_column_counter)=1;
                    HiPerDiF_count_nc_B=sum(sum(RVE(((i-1)*HiPerDiF_rows_resolution+1):(i*HiPerDiF_rows_resolution),((j-1)*HiPerDiF_columns_resolution+1):(j*HiPerDiF_columns_resolution))));
                end
            else
                error('Incorrect HiPerDiF RVE configuration');
            end
        end
    end
    
elseif x_section_arrangement==6
    %     RVE(1:n_rows,1:112)=0;
    %     RVE(1:n_rows,113:167)=1;
    %     RVE(1:n_rows,168:n_columns)=0;
    
    % for Vc of 0.2
    RVE(1:n_rows,1:124)=0;
    RVE(1:n_rows,125:154)=1;
    RVE(1:n_rows,155:n_columns)=0;
    
    % Dispersed 3
elseif x_section_arrangement==7
    %         % test
    %         RVE(1,:)=0;
    %         RVE(2:3,1:73)=0;
    %         RVE(2:3,74:103)=1;
    %         RVE(2:3,104:175)=0;
    %         RVE(2:3,176:206)=1;
    %         RVE(2:3,207:279)=0;
    %         RVE(4,:)=0;
    
    g_height=floor(HiPerDiF_rows_resolution);
    c_height=floor(2*HiPerDiF_rows_resolution+1);
    
    c_width=round((nc/(c_height))/2);
    g_width=floor((n_columns-2*c_width)/3);

    
    RVE(1:g_height,1:end)=0;
    RVE(g_height+1:g_height+c_height,1:g_width)=0;
    RVE(g_height+1:g_height+c_height,g_width+1:(g_width+c_width))=1;
    RVE(g_height+1:g_height+c_height,(g_width+c_width+1):(2*g_width+c_width))=0;
    RVE(g_height+1:g_height+c_height,(2*g_width+c_width+1):(2*g_width+2*c_width))=1;
    RVE(g_height+1:g_height+c_height,(2*g_width+2*c_width+1):n_columns)=0;
    RVE(g_height+c_height+1:end,1:end)=0;
    
    % Dispersed 2
elseif x_section_arrangement==8
    c_width=ceil(nc/((n_rows)*3));
    g_width=floor((n_columns-3*c_width)/6);
    
    RVE(1:n_rows,(g_width+1):(g_width+c_width))=1;
    RVE(1:n_rows,(3*g_width+c_width+1):(3*g_width+2*c_width))=1;
    RVE(1:n_rows,(5*g_width+2*c_width+1):(5*g_width+3*c_width))=1;
        
    % blocked
elseif x_section_arrangement==9
    c_width=ceil(nc/((n_rows)));
    g_width=floor((n_columns-c_width)/2);
    RVE(1:n_rows,1:g_width)=0;
    RVE(1:n_rows,(g_width+1):(g_width+c_width))=1;
    RVE(1:n_rows,(g_width+c_width+1):end)=0;
    
    
    % Dispersed 1 - chequered
elseif x_section_arrangement==10
    n_blocks_width=9;
    n_cblocks_width=3;
    n_gblocks_width=3;
    c_width=round(nc/((n_rows*n_cblocks_width)));
    g_width=round((n_rows*n_columns-nc)/((n_rows*n_gblocks_width)));
       
    RVE(((0)*HiPerDiF_rows_resolution+1):(1*HiPerDiF_rows_resolution),1:c_width)=1;
    RVE(((0)*HiPerDiF_rows_resolution+1):(1*HiPerDiF_rows_resolution),(c_width+g_width+1):(2*c_width+g_width))=1;
    RVE(((0)*HiPerDiF_rows_resolution+1):(1*HiPerDiF_rows_resolution),(2*c_width+2*g_width+1):(3*c_width+2*g_width))=1;
    
    RVE(((1)*HiPerDiF_rows_resolution+1):(2*HiPerDiF_rows_resolution+1),(2/3*(g_width+c_width)+1):(c_width+2/3*(g_width+c_width)))=1;
    RVE(((1)*HiPerDiF_rows_resolution+1):(2*HiPerDiF_rows_resolution+1),(c_width+g_width+2/3*(g_width+c_width)+1):(2*c_width+g_width+2/3*(g_width+c_width)))=1;
    RVE(((1)*HiPerDiF_rows_resolution+1):(2*HiPerDiF_rows_resolution+1),(2*c_width+2*g_width+2/3*(g_width+c_width)+1):(3*c_width+2*g_width+2/3*(g_width+c_width)))=1;
    
    RVE(((2)*HiPerDiF_rows_resolution+2):(3*HiPerDiF_rows_resolution+2),(1/3*(g_width+c_width)+1):(1/3*(g_width+c_width)+c_width))=1;
    RVE(((2)*HiPerDiF_rows_resolution+2):(3*HiPerDiF_rows_resolution+2),(c_width+g_width+1/3*(g_width+c_width)+1):(2*c_width+g_width+1/3*(g_width+c_width)))=1;
    RVE(((2)*HiPerDiF_rows_resolution+2):(3*HiPerDiF_rows_resolution+2),(2*c_width+2*g_width+1/3*(g_width+c_width)+1):(3*c_width+2*g_width+1/3*(g_width+c_width)))=1;
        
    RVE(((3)*HiPerDiF_rows_resolution+3):(4*HiPerDiF_rows_resolution+2),1:c_width)=1;
    RVE(((3)*HiPerDiF_rows_resolution+3):(4*HiPerDiF_rows_resolution+2),(c_width+g_width+1):(2*c_width+g_width))=1;
    RVE(((3)*HiPerDiF_rows_resolution+3):(4*HiPerDiF_rows_resolution+2),(2*c_width+2*g_width+1):(3*c_width+2*g_width))=1;  
    
        
    % Horizontal stripes
elseif x_section_arrangement==11
    c_height=floor(nc/n_columns);
    remaining_c=nc-n_columns*c_height;
    g_height=floor((n_rows-c_height)/2);
    RVE(g_height+1:g_height+c_height,:)=1;
    RVE(g_height+c_height+1,1:remaining_c)=1;
    
    % FJ - larger size vertical stripes
    
elseif x_section_arrangement==12
    %     n_columns_c=3;
    %     c_width=floor(nc/(n_rows*n_columns_c));
    g_width=ceil(((n_rows*n_columns-nc)/(n_columns_c+1))/n_rows);
    column_counter=g_width;
    for i=1:n_columns_c
        RVE(:,column_counter+1:column_counter+c_width)=1;
        column_counter=column_counter+c_width+g_width;
    end
    
    remaining_c=nc-sum(sum(RVE(:,:)));
    remaining_c_columns=floor(remaining_c/n_rows);
    %
    
    
    %     if remaining_c>0
    %
    %
    %         if remaining_c_columns>0
    %             columns_to_add=floor(remaining_c_columns/n_columns_c);
    %             column_counter2=g_width+c_width+1;
    %             for i=1:n_columns_c
    %                 RVE(:,column_counter2:column_counter2+columns_to_add-1)=1;
    %                 column_counter2=column_counter2+g_width+c_width+columns_to_add;
    %             end
    %             RVE(1:(remaining_c-n_rows*remaining_c_columns),column_counter2)=1;
    %         else
    %             RVE(1:remaining_c,column_counter2-c_width-g_width+1)=1;
    %         end
    %     end
    
elseif x_section_arrangement ==13
    block_columns=n_columns/c_width;
    block_counter=1;
    for i=1:block_columns
        
        if mod(i,2)==0
            RVE(1:c_width,((block_counter-1)*c_width+1):(block_counter)*c_width)=1;
            block_counter=block_counter+1;
        elseif mod(i,2)==1
            RVE(c_width+1:2*c_width,((block_counter-1)*c_width+1):(block_counter)*c_width)=1;
            block_counter=block_counter+1;
        end
    end
    
    for j=3:block_columns
        if mod(j,2)==1
            RVE((j-1)*c_width+1:j*c_width,:)=RVE(1:c_width,:);
        elseif mod(j,2)==0
            RVE((j-1)*c_width+1:j*c_width,:)=RVE(c_width+1:2*c_width,:);
        end
    end
    
elseif x_section_arrangement ==14
    gap_counter=round((1+(nc/(n_rows*n_columns-nc)))/(nc/(n_rows*n_columns-nc)));
    initial_gap=round(rand()*4);
    for j=1:n_columns
       RVE(1+initial_gap:gap_counter:end,j)=1;
       initial_gap=mod(initial_gap+2,gap_counter);        
    end
    
elseif x_section_arrangement==15
    RVE(4:4:12,:)=1;

elseif x_section_arrangement==16
    load('RVE_spreadsheet.mat')
    RVE=RVE_spreadsheet;
    
% blocks of 4 carbon fibres    
elseif x_section_arrangement==17
%     load('RVE_4blocks.mat')
%     % for when RVE is 18x387
%     
%     dbstop if error
%     
% %     rows_offset=round(rand()*4);
%     rows_offset=2;
%     
% %     rows_offset=4;
%     
%     Basic_RVE=zeros(n_rows,8);
%     Basic_RVE(1+rows_offset:2+rows_offset,1:2)=1;
%     Basic_RVE(9+rows_offset:10+rows_offset,1:2)=1;
%     Basic_RVE(:,5:6)=Basic_RVE(:,1:2);
%     Basic_RVE(5+rows_offset:6+rows_offset,3:4)=1;
%     Basic_RVE(:,7:8)=Basic_RVE(:,3:4);
%     Basic_RVE(13+rows_offset:14+rows_offset,3:4)=1;
%          
%     RVE=[zeros(n_rows,1) repmat(Basic_RVE,1,floor(n_columns/8)) zeros(n_rows,1)];
%     
%     RVE(1+rows_offset:2+rows_offset,end-1:end)=1;
%     
%     RVE=[RVE zeros(n_rows,1)];
    
    RVE(3:4,7:4:end)=1;
    RVE(3:4,8:4:end)=1;
    RVE(11:12,7:4:end)=1;
    RVE(11:12,8:4:end)=1;
    
    RVE(7:8,3:5:end)=1;
    RVE(7:8,4:5:end)=1;
    RVE(15:16,3:5:end)=1;
    RVE(15:16,4:5:end)=1;
    
    RVE(:,end-4:end)=0;
    RVE(7:8,3:4)=0;
    RVE(15:16,3:4)=0;
    RVE(3:4,end-5:end-4)=0;
    RVE(11:12,end-9:end)=0;
        
elseif x_section_arrangement==18
%      load('RVE_9blocks.mat')
%     % for when RVE is 18x387
%     rows_offset=round(rand()*4);
%     columns_offset=round(rand()*4);
    rows_offset=2;
    columns_offset=2;


    RVE=zeros(n_rows,n_columns);
    RVE(1+rows_offset:3+rows_offset,[5+columns_offset:7:end-10+columns_offset 6+columns_offset:7:end-10+columns_offset 7+columns_offset:7:end-7+columns_offset])=1;
    RVE(6+rows_offset:8+rows_offset,[1+columns_offset:8:end 2+columns_offset:8:end 3+columns_offset:8:end])=1;
    RVE(12+rows_offset:14+rows_offset,[5+columns_offset:8:end 6+columns_offset:8:end 7+columns_offset:8:end])=1;
    
    RVE(:,end)=0;
    
% thin ply    
elseif x_section_arrangement==19
    NumPlies=floor(nc/n_columns);
    PlySpacing=round(n_rows/NumPlies)-1;
    RVE(PlySpacing:PlySpacing:NumPlies*PlySpacing,:)=1;
    
% practical perfectly spaced (checkerboard limitation)
elseif x_section_arrangement==20
    % FJ - Create a matrix with a random placement of glass and carbon fibres
    % FJ - start by defining a matrix with all glass fibres
    RVE=zeros(n_rows,n_columns); %  Glass = 0  -  Carbon = 1
    % FJ - replace a fibre with a carbon fibre until the number of carbon
    % fibres matches nc
    while sum(sum(RVE))~=nc
        fibre_location=[ceil(n_rows*rand()),ceil(n_columns*rand())];
        
        % if proposed fibre location doesn't have a carbon neighbour on any
        % of its sides...
        if RVE(max([1 fibre_location(1)-1]),fibre_location(2))==0 && RVE(min([n_rows fibre_location(1)+1]),fibre_location(2))==0 && RVE(fibre_location(1),max([1 fibre_location(2)-1]))==0 && RVE(fibre_location(1),min([n_columns fibre_location(2)+1]))==0
            RVE(fibre_location(1),fibre_location(2))=1;
        end
    end
    % FJ - still not sure what this means
    % RVE=RVEg;
    
    % practical perfectly spaced (agglomerated fibres limitation)
elseif x_section_arrangement==21
    % FJ - Create a matrix with a random placement of glass and carbon fibres
    % FJ - start by defining a matrix with all glass fibres
    RVE=zeros(n_rows,n_columns); %  Glass = 0  -  Carbon = 1
    % FJ - replace a fibre with a carbon fibre until the number of carbon
    % fibres matches nc
    while sum(sum(RVE))~=nc
        fibre_location=[ceil(n_rows*rand()),ceil(n_columns*rand())];
        
        % if proposed fibre location doesn't have a carbon neighbour on any
        % of its sides...
        if RVE(max([1 fibre_location(1)-1]),fibre_location(2))==0 && RVE(min([n_rows fibre_location(1)+1]),fibre_location(2))==0 && RVE(fibre_location(1),max([1 fibre_location(2)-1]))==0 && RVE(fibre_location(1),min([n_columns fibre_location(2)+1]))==0 && ...
            RVE(max([1 fibre_location(1)-2]),fibre_location(2))==0 && RVE(min([n_rows fibre_location(1)+2]),fibre_location(2))==0 && RVE(fibre_location(1),max([1 fibre_location(2)-2]))==0 && RVE(fibre_location(1),min([n_columns fibre_location(2)+2]))==0 && ...
            RVE(max([1 fibre_location(1)-1]),max([1 fibre_location(2)-1]))==0 && RVE(min([n_rows fibre_location(1)+1]),max([1 fibre_location(2)-1]))==0 && RVE(max([1 fibre_location(1)-1]),min([n_columns fibre_location(2)+1]))==0 && RVE(min([n_rows fibre_location(1)+1]),min([n_columns fibre_location(2)+1]))==0
        
            RVE(fibre_location(1),fibre_location(2))=1;
        end
    end
    % FJ - still not sure what this means
    % RVE=RVEg;
    
% perfectly spaced with border of glass    
elseif x_section_arrangement ==22
    gap_counter=round((1+(nc/(n_rows*n_columns-nc)))/(nc/(n_rows*n_columns-nc)));
    initial_gap=round(rand()*4);
    for j=2:n_columns-1
       RVE(2+initial_gap:gap_counter:end-1,j)=1;
       initial_gap=mod(initial_gap+2,gap_counter);        
    end
    
    % perfectly spaced with double border of glass    
elseif x_section_arrangement ==23
    gap_counter=round((1+(nc/(n_rows*n_columns-nc)))/(nc/(n_rows*n_columns-nc)));
    initial_gap=round(rand()*4);
    for j=3:n_columns-2
       RVE(3+initial_gap:gap_counter:end-2,j)=1;
       initial_gap=mod(initial_gap+2,gap_counter);        
    end
    % chequered configuration, but with block corners touching
elseif x_section_arrangement ==24
    n_blocks_width=9;
    n_cblocks_width=3;
    n_gblocks_width=3;
    c_width=round(nc/((n_rows*n_cblocks_width)));
    g_width=round((n_rows*n_columns-nc)/((n_rows*n_gblocks_width)));

    RVE(((0)*HiPerDiF_rows_resolution+1):(1*HiPerDiF_rows_resolution),1:c_width)=1;
    RVE(((0)*HiPerDiF_rows_resolution+1):(1*HiPerDiF_rows_resolution),(c_width+g_width+1):(2*c_width+g_width))=1;
    RVE(((0)*HiPerDiF_rows_resolution+1):(1*HiPerDiF_rows_resolution),(2*c_width+2*g_width+1):(3*c_width+2*g_width))=1;
    
    RVE(((1)*HiPerDiF_rows_resolution+1):(2*HiPerDiF_rows_resolution+1),((g_width)+1):((g_width+c_width)))=1;
    RVE(((1)*HiPerDiF_rows_resolution+1):(2*HiPerDiF_rows_resolution+1),((2*g_width+c_width)+1):(2*(g_width+c_width)))=1;
    RVE(((1)*HiPerDiF_rows_resolution+1):(2*HiPerDiF_rows_resolution+1),(3*g_width+2*c_width+1):(3*(g_width+c_width)))=1;
    
    RVE(((2)*HiPerDiF_rows_resolution+2):(3*HiPerDiF_rows_resolution+2),((g_width-c_width)+1):((g_width)))=1;
    RVE(((2)*HiPerDiF_rows_resolution+2):(3*HiPerDiF_rows_resolution+2),((2*g_width)+1):(2*g_width+c_width))=1;
    RVE(((2)*HiPerDiF_rows_resolution+2):(3*HiPerDiF_rows_resolution+2),(3*g_width+c_width+1):(3*g_width+2*c_width))=1;
        
    RVE(((3)*HiPerDiF_rows_resolution+3):(4*HiPerDiF_rows_resolution+2),((g_width-2*c_width)+1):((g_width-c_width)))=1;
    RVE(((3)*HiPerDiF_rows_resolution+3):(4*HiPerDiF_rows_resolution+2),((2*g_width-c_width)+1):(2*g_width))=1;
    RVE(((3)*HiPerDiF_rows_resolution+3):(4*HiPerDiF_rows_resolution+2),(3*g_width+1):(3*g_width+c_width))=1;
    
% 16-block cross-section
elseif x_section_arrangement ==25
    RVE(2:5,9:18:end)=1;
    RVE(2:5,10:18:end)=1;
    RVE(2:5,11:18:end)=1;
    RVE(2:5,12:18:end)=1;
    
    RVE(6:9,18:18:end)=1;
    RVE(6:9,19:18:end)=1;
    RVE(6:9,20:18:end)=1;
    RVE(6:9,21:18:end)=1;
    
    RVE(10:13,9:18:end)=1;
    RVE(10:13,10:18:end)=1;
    RVE(10:13,11:18:end)=1;
    RVE(10:13,12:18:end)=1;
    
    RVE(14:17,18:18:end)=1;
    RVE(14:17,19:18:end)=1;
    RVE(14:17,20:18:end)=1;
    RVE(14:17,21:18:end)=1;
    
    RVE(:,end-5:end)=0;
    
% 25-block cross-section
elseif x_section_arrangement ==26
    RVE(2:6,28:20:end)=1;
    RVE(2:6,29:20:end)=1;
    RVE(2:6,30:20:end)=1;
    RVE(2:6,31:20:end)=1;
    RVE(2:6,32:20:end)=1;
    
    RVE(7:11,18:20:end-20)=1;
    RVE(7:11,19:20:end-20)=1;
    RVE(7:11,20:20:end-20)=1;
    RVE(7:11,21:20:end-20)=1;
    RVE(7:11,22:20:end-20)=1;
    
    RVE(12:16,28:20:end)=1;
    RVE(12:16,29:20:end)=1;
    RVE(12:16,30:20:end)=1;
    RVE(12:16,31:20:end)=1;
    RVE(12:16,32:20:end)=1;
    
% all carbon fibre 
elseif x_section_arrangement == 27
    RVE=ones(n_rows,n_columns);
    
% isolated fibre with larger carbon fibres 
elseif x_section_arrangement ==28
    RVE(1:2:end,1:2:end)=1;
    RVE(2:2:end,2:2:end)=1;
    
    % thin-ply with larger carbon fibres 
elseif x_section_arrangement ==29
    NumPlies=round(nc/n_columns);
    PlySpacing=round(n_rows/NumPlies);
    RVE(PlySpacing:PlySpacing:NumPlies*PlySpacing,:)=1;
       
else
    error('invalid RVE x-section arrangement selection');
end

% if size(RVE,1)~=n_rows
%     error('incorrect row dimensions for HiPerDiF arrangement');
% elseif size(RVE,2)~=n_columns
%     error('incorrect row dimensions for HiPerDiF arrangement');
% end
if Vc==0
    RVE=zeros(n_rows,n_columns);
end

% figure('units','normalized','position',[0.1 0.3 0.7 0.1])
% % figure('units','normalized','position',[0.1 0.7 0.7 0.1])
% imagesc(RVE);
% colormap(flipud(colormap))
% axis equal tight

if HiPerDiF_fuzzyness==0
    1;

elseif HiPerDiF_fuzzyness==1
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    while sum(sum(RVE))>nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    
elseif HiPerDiF_fuzzyness==2
    % ensure accurate Vc by adding random carbons
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    % Change a percentage of fibres to glass
    while sum(sum(RVE))>round((1-Fuzzy_Factor)*nc)
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    % change a percentage of fibres to carbon - nc should be preserved
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end

% Accurate fuzzy factor which doesn't replace fibres that have moved
% already
elseif HiPerDiF_fuzzyness==3
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
    % ensure accurate Vc by adding random carbons
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    while sum(sum(RVE))>nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            % choose a random position
            change_position_glass=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end
    
% Accurate, row-specific fuzzy factor    
elseif HiPerDiF_fuzzyness==4
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
    % ensure accurate Vc by adding random carbons
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    while sum(sum(RVE))>nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            % choose a random position
            change_position_glass=[change_position_carbon(1,1),ceil(n_columns*rand())];
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end    
    
% Accurate, row-specific fuzzy factor with some row movement within rows
% around RVE block
elseif HiPerDiF_fuzzyness==5
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
    % ensure accurate Vc by adding random carbons
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    while sum(sum(RVE))>nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            % choose a random position
            change_position_glass=[(change_position_carbon(1,1)~=n_rows)*(floor(((change_position_carbon(1,1)-1)/HiPerDiF_rows_resolution))*HiPerDiF_rows_resolution + ceil(rand()*4))+(change_position_carbon(1,1)==n_rows)*change_position_carbon(1,1),ceil(n_columns*rand())];
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end    
    
    % Accurate, row-specific fuzzy factor with some row movement around
    % original row position
elseif HiPerDiF_fuzzyness==6
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
    % ensure accurate Vc by adding random carbons
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    while sum(sum(RVE))>nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            % choose a random position
            vertical_migration_std_dev=HiPerDiF_rows_resolution;
            change_position_glass=[ceil(abs((change_position_carbon(1,1)-(vertical_migration_std_dev/2)+rand()*vertical_migration_std_dev))),ceil(n_columns*rand())];
            change_position_glass(1,1)=(change_position_glass(1,1)>n_rows)*(n_rows-(change_position_glass(1,1)-n_rows))+(change_position_glass(1,1)<=n_rows)*change_position_glass(1,1);      
            
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end  
    
    % Fuzziness based on limited carbon fibre migration through-thickness
    % and through-width
    elseif HiPerDiF_fuzzyness==7
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
    % ensure accurate Vc by adding random carbons
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    while sum(sum(RVE))>nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            % choose a random position
            vertical_migration_std_dev=HiPerDiF_rows_resolution;
            horizontal_migration_std_dev=HiPerDiF_columns_resolution*5;
            change_position_glass=[ceil(abs((change_position_carbon(1,1)-(vertical_migration_std_dev/2)+rand()*vertical_migration_std_dev))),ceil(abs((change_position_carbon(1,2)-(horizontal_migration_std_dev/2)+rand()*horizontal_migration_std_dev)))];
            change_position_glass(1,1)=(change_position_glass(1,1)>n_rows)*(n_rows-(change_position_glass(1,1)-n_rows))+(change_position_glass(1,1)<=n_rows)*change_position_glass(1,1);      
            change_position_glass(1,2)=(change_position_glass(1,2)>n_columns)*(n_columns-(change_position_glass(1,2)-n_columns))+(change_position_glass(1,2)<=n_columns)*change_position_glass(1,2);      
            
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end 

        % Fuzziness based on normal random number distribution (verical and
        % horizontal)
    elseif HiPerDiF_fuzzyness==8
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
    % ensure accurate Vc by adding random carbons
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    while sum(sum(RVE))>nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            change_position_glass=[0 0];
            while change_position_glass(1,1)>n_rows || change_position_glass(1,1)<1 || change_position_glass(1,2)>n_columns || change_position_glass(1,2)<1
            % choose a random position
            vertical_migration_std_dev=HiPerDiF_rows_resolution/4;
            horizontal_migration_std_dev=HiPerDiF_columns_resolution*1.25;
            change_position_glass=[ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev)))),ceil(abs((change_position_carbon(1,2)+normrnd(0,horizontal_migration_std_dev))))];
%             change_position_glass(1,1)=(change_position_glass(1,1)>n_rows)*(n_rows-(change_position_glass(1,1)-n_rows))+(change_position_glass(1,1)<=n_rows)*change_position_glass(1,1);      
%             change_position_glass(1,2)=(change_position_glass(1,2)>n_columns)*(n_columns-(change_position_glass(1,2)-n_columns))+(change_position_glass(1,2)<=n_columns)*change_position_glass(1,2);     
            end
            
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end
    
    % normal distro with strata
elseif HiPerDiF_fuzzyness==9
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
    % ensure accurate Vc by adding random carbons
    while sum(sum(RVE))<nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
    end
    
    while sum(sum(RVE))>nc
        RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
    end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
%     width_norm_mean=round(normrnd(0,HiPerDiF_columns_resolution*2,1,n_rows));
    width_norm_mean=ceil((2*rand(1,n_rows)-1)*HiPerDiF_columns_resolution*2.5);
    width_norm_mean(2:2:end)=width_norm_mean(1:2:end);
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            change_position_glass=[0 0];
            while_counter=0;
            std_dev_factor=0.75;
            while change_position_glass(1,1)>n_rows || change_position_glass(1,1)<1 || change_position_glass(1,2)>n_columns || change_position_glass(1,2)<1
            % choose a random position
            vertical_migration_std_dev=(HiPerDiF_rows_resolution-2)/4;
            horizontal_migration_std_dev=31*std_dev_factor;
            change_position_glass=[ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev)))),ceil(abs((change_position_carbon(1,2)+normrnd(width_norm_mean(1,change_position_carbon(1,1)),horizontal_migration_std_dev))))];
          
            % make sure you don't get stuck in a loop if the fibres in the
            % corner can't be moved
            while_counter=while_counter+1;
            if while_counter>(n_rows*n_columns)
                std_dev_factor=std_dev_factor*1.1;
%                 width_norm_mean(1,change_position_carbon(1,1))=ceil((2*rand()-1)*HiPerDiF_columns_resolution*2.5);
            end
%             change_position_glass(1,1)=ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev))));
%             
%             change_position_glass(1,2)=ceil(abs((change_position_carbon(1,2)+normrnd(width_norm_mean,horizontal_migration_std_dev))));
            
            end
            
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end 
    
        % normal distro with strata
elseif HiPerDiF_fuzzyness==10
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
%     % ensure accurate Vc by adding random carbons
%     while sum(sum(RVE))<nc
%         RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
%     end
%     
%     while sum(sum(RVE))>nc
%         RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
%     end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
%     width_norm_mean=round(normrnd(0,HiPerDiF_columns_resolution*2,1,n_rows));
    %width_norm_mean=ceil((2*rand(1,n_rows)-1)*HiPerDiF_columns_resolution*2.5);
    width_norm_mean=round(normrnd(0,n_columns/3,[1 n_rows]));
    width_norm_mean(2:2:end)=width_norm_mean(1:2:end);
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            change_position_glass=[0 0];
            while_counter=0;
            std_dev_factor=0.75;
            
            % choose a random position
            vertical_migration_std_dev=(HiPerDiF_rows_resolution-2)/2;
%             vertical_migration_std_dev=0.0001;
            horizontal_migration_std_dev=31*std_dev_factor;
            change_position_glass=[ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev)))),ceil(abs((change_position_carbon(1,2)+normrnd(width_norm_mean(1,change_position_carbon(1,1)),horizontal_migration_std_dev))))];
            
            while change_position_glass(1,1)>n_rows || change_position_glass(1,1)<1 || change_position_glass(1,2)>n_columns || change_position_glass(1,2)<1
            % choose a random position
%             vertical_migration_std_dev=(HiPerDiF_rows_resolution-2)/4;
%             horizontal_migration_std_dev=31*std_dev_factor;
%             change_position_glass=[ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev)))),ceil(abs((change_position_carbon(1,2)+normrnd(width_norm_mean(1,change_position_carbon(1,1)),horizontal_migration_std_dev))))];
            if change_position_glass(1,1)>n_rows
                change_position_glass(1,1)=n_rows-(change_position_glass(1,1)-n_rows);          
            elseif change_position_glass(1,1)<1
                change_position_glass(1,1)=1-(change_position_glass(1,1)-1);              
            elseif change_position_glass(1,2)>n_columns
                change_position_glass(1,2)=n_columns-(change_position_glass(1,2)-n_columns);    
            elseif change_position_glass(1,2)<1
                change_position_glass(1,2)=1-(change_position_glass(1,2)-1);  
            end
            % make sure you don't get stuck in a loop if the fibres in the
            % corner can't be moved
            while_counter=while_counter+1;
            if while_counter>(n_rows*n_columns)
                error('too many iterations to converge fuzziness model')
                std_dev_factor=std_dev_factor*1.1;
%                 width_norm_mean(1,change_position_carbon(1,1))=ceil((2*rand()-1)*HiPerDiF_columns_resolution*2.5);
            end
%             change_position_glass(1,1)=ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev))));
%             
%             change_position_glass(1,2)=ceil(abs((change_position_carbon(1,2)+normrnd(width_norm_mean,horizontal_migration_std_dev))));
            
            end
            
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end 
        % normal distro with strata
elseif HiPerDiF_fuzzyness==11
    if n_rows*n_columns*Fuzzy_Factor > 2*nc
        error(['Fuzzy_Factor too high for model to converge. Try Fuzzy_Factor lower than ' num2str((2*nc/(n_rows*n_columns)),4)])
    end
    
    % ensure accurate Vc by adding random carbons
%     while sum(sum(RVE))<nc
%         RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=1;
%     end
%     
%     while sum(sum(RVE))>nc
%         RVE(ceil(n_rows*rand()),ceil(n_columns*rand()))=0;
%     end
    
    % create a matrix that records whether the values have been changed (0 = no changes, 1 = a change has been made)
    RVE_changes=zeros(n_rows,n_columns);
    RVE_sum_orig=sum(sum(RVE));
%     width_norm_mean=round(normrnd(0,HiPerDiF_columns_resolution*2,1,n_rows));
    width_norm_mean=round((2*rand(1,n_rows)-1)*HiPerDiF_columns_resolution*2*1);
%     width_norm_mean=round(normrnd(0,n_columns/3,[1 n_rows]));
    width_norm_mean(2:2:end)=width_norm_mean(1:2:end);
    % complete this loop until Fuzzy_factor*100 percent changes have been made
    while sum(sum(RVE_changes))<round(n_rows*n_columns*Fuzzy_Factor)
        while RVE_sum_orig<(sum(sum(RVE))+1)
            % choose a random position
            change_position_carbon=[ceil(n_rows*rand()),ceil(n_columns*rand())];
            % if the fibre is a carbon, and it hasn't been changed before, then change it to a glass
            if RVE(change_position_carbon(1,1),change_position_carbon(1,2))==1 && RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))==0
                RVE(change_position_carbon(1,1),change_position_carbon(1,2))=0;
                % indicate a change has been made in this position
                RVE_changes(change_position_carbon(1,1),change_position_carbon(1,2))=1;
                %                 RVE_sum=RVE_sum+1;
                break
            end
        end
        while sum(sum(RVE))<RVE_sum_orig
            change_position_glass=[0 0];
            while_counter=0;
            std_dev_factor=0.75;
            
            % choose a random position
            vertical_migration_std_dev=(HiPerDiF_rows_resolution/4)*1;
%             vertical_migration_std_dev=0.0001;
            horizontal_migration_std_dev=round(HiPerDiF_columns_resolution*0.75)*std_dev_factor*1;
            change_position_glass=[ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev)))),ceil(abs((change_position_carbon(1,2)+normrnd(width_norm_mean(1,change_position_carbon(1,1)),horizontal_migration_std_dev))))];
            
            while change_position_glass(1,1)>n_rows || change_position_glass(1,1)<1 || change_position_glass(1,2)>n_columns || change_position_glass(1,2)<1
            % choose a random position
%             vertical_migration_std_dev=(HiPerDiF_rows_resolution-2)/4;
%             horizontal_migration_std_dev=31*std_dev_factor;
%             change_position_glass=[ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev)))),ceil(abs((change_position_carbon(1,2)+normrnd(width_norm_mean(1,change_position_carbon(1,1)),horizontal_migration_std_dev))))];
            if change_position_glass(1,1)>n_rows
                change_position_glass(1,1)=n_rows-(change_position_glass(1,1)-n_rows);          
            elseif change_position_glass(1,1)<1
                change_position_glass(1,1)=1-(change_position_glass(1,1)-1);              
            elseif change_position_glass(1,2)>n_columns
                change_position_glass(1,2)=n_columns-(change_position_glass(1,2)-n_columns);    
            elseif change_position_glass(1,2)<1
                change_position_glass(1,2)=1-(change_position_glass(1,2)-1);  
            end
            % make sure you don't get stuck in a loop if the fibres in the
            % corner can't be moved
            while_counter=while_counter+1;
            if while_counter>(n_rows*n_columns)
                error('too many iterations to converge fuzziness model')
                std_dev_factor=std_dev_factor*1.1;
%                 width_norm_mean(1,change_position_carbon(1,1))=ceil((2*rand()-1)*HiPerDiF_columns_resolution*2.5);
            end
%             change_position_glass(1,1)=ceil(abs((change_position_carbon(1,1)+normrnd(0,vertical_migration_std_dev))));
%             
%             change_position_glass(1,2)=ceil(abs((change_position_carbon(1,2)+normrnd(width_norm_mean,horizontal_migration_std_dev))));
            
            end
            
            % if the fibre is a glass, and it hasn't been changed before, then change it to a carbon
%             if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 && RVE_changes(change_position_glass(1,1),change_position_glass(1,2))==0
            if RVE(change_position_glass(1,1),change_position_glass(1,2))==0 
                RVE(change_position_glass(1,1),change_position_glass(1,2))=1;
                % inidicate a change has been made in this position
                RVE_changes(change_position_glass(1,1),change_position_glass(1,2))=1;
                %                 RVE_sum=RVE_sum-1;
                break
            end
        end
    end 
    
else 
    error('Incorrect HiPerDiF_fuzziness setting.')
end


nc=sum(sum(RVE(:,:)));

% % figure();
% figure('units','normalized','position',[0.1 0.3 0.7 0.1])
% % figure('units','normalized','position',[0.1 0.1 0.7 0.7])
% imagesc(1:n_columns,1:n_rows,RVE);
% map=[217/256 217/256 217/256;127/256 127/256 127/256];
% % map=[0 0 0;1 1 204/255];
% colormap(flipud(map))
% %set(gca,'position',[0 0 1 1],'units','normalized')
% axis equal tight
% axis off
end