%% Parameters Defination
ExcuteTimes = ExcuteTimes_t;
ParticleNum = ParticleNum_t;             % Number of Particles
MaxIteration = MaxIteration_t;          % Max Iterations
MaxChaosIteration = MaxChaosIteration_t;
dim = dim_t;                     % Dimension
RangeL = RangeL_t;
RangeR = RangeR_t;
XMax = XMax_t;                    % Range of X
%  func_test = 'f1';
save_str_head = ['Figures-',datestr(now,'dd-mm-yyyy'),'/',func_test,'/','I',num2str(MaxIteration),'_D',num2str(dim),'_M',num2str(ParticleNum)];
mkdir(save_str_head);

%---------Chaos-----
% |
% |___    A
% |   \
% |    \
% |     \
% --------------------->
%    |M  |N           |Maxiteration

%---------beta-----
% |
% |_________    
% |         \
% |          \
% |           \
% --------------------->
%    |Offset|           |Maxiteration
A = A_t;
ParticleVari = A;
M = M_t;
N = N_t;
offset = offset_t;

alpha = alpha_t;
sig = sig_t;
excuted =0;
excuted2 =0;


%% Observing Workspace
data_dcqpso=zeros(1,MaxIteration);   
sigma=zeros(1,MaxIteration);
alpha_mat = zeros(1,MaxIteration);
mposition = zeros(MaxIteration,dim);
p = zeros(ParticleNum,dim,MaxIteration);

%% Initiallation

    
    tentmap = rand(ParticleNum,dim,1);
    tentmap2 = rand(MaxIteration,dim,1);
    tentmap2 = mod(2*tentmap2,1);
    
    x=(RangeR- RangeL)*rand(ParticleNum,dim,1) + RangeL;    % Initialize Swarm
    x2=(RangeR- RangeL)*tentmap2 + RangeL;
    
    pbest=x;    % Private Extremum Particles [size*Dim]
    gbest=zeros(1,dim); % Global Extremum Particle [1*Dim]
    
    mbest = zeros(MaxIteration,dim);
    
    % ffx = zeros(ParticleNum, dim+1);
    
    f_x = zeros(1,ParticleNum);
    f_pbest = zeros(1,ParticleNum);
    
    for i=1:ParticleNum
        f_x(i)=eval([func_test, '(x(i,:))']);  % Calculate f(x) for Particles.
        f_pbest(i)=f_x(i);  % Initialize Private Extremum.
    end
    
    f_x2 = zeros(1,MaxChaosIteration);
    
    g=min(find(f_pbest==min(f_pbest(1:ParticleNum))));  % Refresh Global Extremum Particle's Indice
    gbest=pbest(g,:);       % Refresh Global Best Particle Vector.
    f_gbest=f_pbest(g);     % Refresh Global Best Value.
    
    MINIMUM=f_gbest;
    
    data_dcqpso=zeros(1,MaxIteration);
    sigma=zeros(1,MaxIteration);

%% Start
    tic
    for t=1:MaxIteration

        switch beta_select
            case 1                  % standard linear model.
                if t > (M+offset)*MaxIteration
                    beta=-0.5/(MaxIteration*(1-(M+offset)))*t + 0.5/(1-(M+offset)) +0.5;
                else
                    beta = 1;
                end
        end
       
        ffx = [x pbest tentmap f_pbest' f_x'];
        ffx = sortrows(ffx, 3*dim+2);
        x = ffx(:,1:dim);
        pbest = ffx(:,dim+1:2*dim);
        tentmap = ffx(:,2*dim+1:3*dim);
        f_pbest = (ffx(:,3*dim+1))';
        f_x = (ffx(:,3*dim+2))';
        
        f_average=sum(f_x(1:(ParticleNum-ParticleVari)))/(ParticleNum-ParticleVari);
        sigma(1,t)=sqrt(sum(((f_x(1:(ParticleNum-ParticleVari))-f_average)/f_average).^2));

        mbest(t,:)=sum(pbest(1:(ParticleNum - ParticleVari),:))/(ParticleNum - ParticleVari);   % Average Best Position.

        if t < M*MaxIteration
            ParticleVari = A;
        elseif t > N*MaxIteration
            ParticleVari = 0;
        else
            ParticleVari = round((-A/(MaxIteration*(N-M)))*t + A*N/(N-M));
        end
        
        data_dcqpso(1,t)=MINIMUM;
        
        if t > 1
            alpha_mat(t) = (data_dcqpso(t-1) - data_dcqpso(t))/data_dcqpso(t);
        end
        
        if (t > chaos_range(1)*MaxIteration) && ( alpha_mat(t) < alpha) && (t < chaos_range(2)*MaxIteration) && (sigma(1,t) < sig)
            excuted = excuted +1;
                lamda = 1;
                for k=1:MaxChaosIteration
                    tentmap2(t,:) = mod(2*tentmap2(t,:),1) + 0.001*rand(1,dim,1);
                    lamda = lamda * 0.5;
                    x2(k,:) = gbest + lamda * ((RangeR- RangeL) * tentmap2(t,:) + RangeL);
                    if (find(abs(x2(k,:))>=XMax))
                        temp_logic = abs(x2(k,:)) >= XMax;
                        temp_r = (RangeR- RangeL)*rand(1,dim,1) + RangeL;
                        temp = x2(k,:);
                        temp(temp_logic) = temp_r(temp_logic);
                        %                     x2(abs(x2)>=XMax) = (RangeR- RangeL)*rand() + RangeL;
                        x2(k,:) = temp;
                    end
                    f_x2(k) = eval([func_test, '(x2(k,:))']);
                    if f_x2(k)<f_gbest
                        excuted2 = excuted2 + 1;
                        f_gbest = f_x2(k);
                        gbest = x2(k,:);
                        continue;
                    end
%                      lamda = lamda * 0.93;
                end
        end
        
        for i=1:ParticleNum  
            
            if i <= (ParticleNum - ParticleVari)
%                 fi=rand(1,dim);
                fi=betarnd(b1,b2,[1 dim]);
%                 p=(b1*fi.*pbest(i,:)+b2*(1-fi).*gbest);
%                 p=sign(p).*min(abs(p),XMax);
                p(i,:,t)= c1 * fi.*pbest(i,:) + c2 * (1-fi).*gbest;
                if (find(abs(p(i,:,t))>=XMax))
                    temp = p(i,:,t);
                    temp_logical = abs(temp)>=XMax;
                    temp_r = (RangeR- RangeL)*rand(1,dim,1) + RangeL;
                    temp(temp_logical) = temp_r(temp_logical);
                    p(i,:,t) = temp;
                end
                u=rand(1,dim);
                b=beta*abs(mbest(t,:)-x(i,:));
                v=-log(u);
                y=p(i,:,t)+((-1).^ceil(0.5+rand(1,dim))).*b.*v;
                if (find(abs(y)>=XMax))
                    %                 x(i,:)=sign(y).*min(abs(y),XMax);   % Refresh Particle Location and Restriction
                    temp_logic = abs(y)>=XMax;
                    temp_r = (RangeR- RangeL)*rand(1,dim,1) + RangeL;
                    y(temp_logic) = temp_r(temp_logic);
                end
                x(i,:)=y;   % Refresh Particle Location and Restriction
            else
                tentmap(i,:) = mod(2*tentmap(i,:),1) + 0.001 * rand(1,dim,1);
                x(i,:) = mbest(t,:) + 0.1 * ParticleVari * (RangeR- RangeL)*tentmap(i,:) + RangeL;
%                 x(i,:) = sign(x(i,:)).*min(abs(x(i,:)),XMax);
%                 x(i,:) = (RangeR- RangeL)*((tentmap(i,:)+1)/2) + RangeL;
%                 x(i,:) = mbest(t,:) + 0.8 * ((RangeR- RangeL)*((tentmap(i,:)+1)/2) + RangeL);
%                 temp = x(i,:);
%                 temp(abs(temp)>=XMax) = (RangeR- RangeL)*rand() + RangeL;
%                 x(i,:) = temp;
                if (find(abs(x(i,:))>=XMax))
                    temp = x(i,:);
                    temp2 = (RangeR- RangeL) * rand(1,dim,1) + RangeL;
                    temp_logic = abs(temp)>=XMax;
                    temp(temp_logic) = temp2(temp_logic);
                    x(i,:) = temp;
                end
            end
                f_x(i)=eval([func_test, '(x(i,:))']);
                if f_x(i)<f_pbest(i)
                    pbest(i,:)=x(i,:);
                    f_pbest(i)=f_x(i);
                end
                if f_pbest(i)<f_gbest
                    gbest=pbest(i,:);
                    f_gbest=f_pbest(i);
                end            
                MINIMUM=f_gbest;   
        end
%           data_dcqpso(1,t)=MINIMUM;
    end
    Ans_dcqpso(2,z) = toc;
    Ans_dcqpso(1,z) = MINIMUM;
    
    Curve_All_dcqpso(z,:) = data_dcqpso;
    
    % Save data_dcqpso figure
    figure('NumberTitle', 'off', 'Name', 'Minimum and Iterations');
    hold on;
    xlabel('Iteration');
    ylabel('f(x)');
    plot(data_dcqpso);
    hold off;
    save_str = [save_str_head, '/', 'data_', num2str(z), '_DCQPSO', '.fig'];
    saveas(gcf, save_str);
    close(gcf);
    
    % Save sigma figure
    figure('NumberTitle', 'off', 'Name', 'Sigma and Iterations');
    hold on;
    xlabel('Iteration');
    ylabel('sigma');
    plot(sigma);
    hold off;
    save_str = [save_str_head, '/', 'sigma_', num2str(z), '_DCQPSO', '.fig'];
    saveas(gcf, save_str);
    close(gcf);
    
    % Save alpha figure
    figure('NumberTitle', 'off', 'Name', 'Alpha and Iterations');
    hold on;
    xlabel('Iteration');
    ylabel('Alpha');
    plot(alpha_mat);
    hold off;
    save_str = [save_str_head, '/', 'alpha_', num2str(z), '_DCQPSO', '.fig'];
    saveas(gcf, save_str);
    close(gcf);
    
    % Save workspace .mat and Hash Encrypt
    save_str = [save_str_head, '/', 'workspace_data_', num2str(z), '_DCQPSO'];
    save(save_str);
    ToHash = [save_str, '.mat'];
    save_str = [save_str_head, '/', 'workspace_data_', num2str(z), '_DCQPSO', '.hash'];
    fid = fopen(save_str, 'w');
    fprintf(fid,'Hash SHA-256 encrypt saved workspace file,"workspace_data_*_DCQPSO.mat", :\r\n');
    fprintf(fid,'%s \r\n',DataHash(ToHash,Opt));
    fclose(fid);


    
    
 