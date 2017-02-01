

%% Sphere Test

for aa = 1:2
    for bb = 1:3
        for cc = 1:3
            
            clearvars -except aa bb cc
            Opt.Method = 'SHA-256';
            Opt.Format = 'HEX';
            Opt.Input = 'file';
            ExcuteTimes_t = 50;
            func_test = 'f1';
            func_name = 'Sphere';
            ParticleNum_t=20*cc;             % Number of Particles
            MaxIteration_t=500+aa*500;          % Max Iterations
            dim_t=10*bb+10;                     % Dimension
            RangeL_t=-50;
            RangeR_t=50;
            XMax_t=50;
            MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
            MaxVariFactor = 0.1;                            % Max Chaos Vari Ratio
            M_t = 0.05;                                      % Chaos Slides Control Factor
            N_t = 0.3;                                      % Chaos Slides Control Factor
            offset_t = -1 * M_t;                                   % Offset between beta and chaos.
            A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
            alpha_t = 0.1;                               % Convergence Judgement factor
            beta_select = 1;                                % QPSO Beta factor mode.
            chaos_range = [0.5 0.9];                        % Chaos range
            sig_t = 10;                                    % Convergence Judgement factor
            b1 = 1;
            b2 = 7;
            c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
        end
    end
end


%% Rosenbrock Test

for aa = 1:2
    for bb = 1:3
        for cc = 1:3
            
            clearvars -except aa bb cc
            Opt.Method = 'SHA-256';
            Opt.Format = 'HEX';
            Opt.Input = 'file';
            ExcuteTimes_t = 50;
            func_test = 'f2';
            func_name = 'Rosenbrock';
            ParticleNum_t=20*cc;             % Number of Particles
            MaxIteration_t=500+aa*500;          % Max Iterations
            dim_t=10*bb+10;                     % Dimension
            RangeL_t=-20;
            RangeR_t=20;
            XMax_t=20;
            MaxChaosIteration_t = round(1/25 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
            MaxVariFactor = 0.2;                            % Max Chaos Vari Ratio
            M_t = 0.1;                                      % Chaos Slides Control Factor
            N_t = 0.3;                                      % Chaos Slides Control Factor
            offset_t = -M_t;                                   % Offset between beta and chaos.s
            A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Fassctor
            alpha_t = 1;                               % Convergence Judgement factor
            beta_select = 1;                                % QPSO Beta factor mode.
            chaos_range = [0.5 0.9];                        % Chaos range
            sig_t = 8;                                      % Convergence Judgement factor
            b1 = 1;
            b2 = 6;
            c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
        end
    end
end


%% Rastrigrin Test

for aa = 1:2
    for bb = 1:3
        for cc = 1:3
            
            clearvars -except aa bb cc
            Opt.Method = 'SHA-256';
            Opt.Format = 'HEX';
            Opt.Input = 'file';
            ExcuteTimes_t = 50;
            func_test = 'f3';
            func_name = 'Rastrigrin';
            ParticleNum_t=20*cc;             % Number of Particles
            MaxIteration_t=500+aa*500;          % Max Iterations
            dim_t=10*bb+10;                     % Dimension
            RangeL_t=-5.12;
            RangeR_t=5.12;
            XMax_t=5.12;
            MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
            MaxVariFactor = 0.2;                            % Max Chaos Vari Ratio
            M_t = 0.1;                                      % Chaos Slides Control Factor
            N_t = 0.3;                                      % Chaos Slides Control Factor
            offset_t = -M_t;                                   % Offset between beta and chaos.
            A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
            alpha_t = 1;                               % Convergence Judgement factor
            beta_select = 1;                                % QPSO Beta factor mode.
            chaos_range = [0.5 0.9];                        % Chaos range
            sig_t = 10;                                      % Convergence Judgement factor
            b1 = 2;
            b2 = 4;
            c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
        end
    end
end


%% Griewank Test

for aa = 1:2
    for bb = 1:3
        for cc = 1:3
            
            clearvars -except aa bb cc
            Opt.Method = 'SHA-256';
            Opt.Format = 'HEX';
            Opt.Input = 'file';
            ExcuteTimes_t = 50;
            func_test = 'f4';
            func_name = 'Griewank';
            ParticleNum_t=20*cc;             % Number of Particles
            MaxIteration_t=500+aa*500;          % Max Iterations
            dim_t=10*bb+10;                     % Dimension
            RangeL_t=-600;
            RangeR_t=600;
            XMax_t=600;
            MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);    % Max Chaos Iteration
            MaxVariFactor = 0.4;                            % Max Chaos Vari Ratio
            M_t = 0.1;                                      % Chaos Slides Control Factor
            N_t = 0.3;                                      % Chaos Slides Control Factor
            offset_t = -M_t;                                   % Offset between beta and chaos.
            A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
            alpha_t = 0.1;                               % Convergence Judgement factor
            beta_select = 1;                                % QPSO Beta factor mode.
            chaos_range = [0.5 0.9];                        % Chaos range
            sig_t = 10;                                    % Convergence Judgement factor
            b1 = 2;
            b2 = 4;
            c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
        end
    end
end


%% Ackley Test
% 
% for aa = 1:2
%     for bb = 1:3
%         for cc = 1:3
%             
%             clearvars -except aa bb cc
%             Opt.Method = 'SHA-256';
%             Opt.Format = 'HEX';
%             Opt.Input = 'file';
%             ExcuteTimes_t = 50;
%             func_test = 'f5';
%             func_name = 'Ackley';
%             ParticleNum_t=20*cc;             % Number of Particles
%             MaxIteration_t=500+aa*500;          % Max Iterations
%             dim_t=10*bb+10;                     % Dimension
%             RangeL_t=-32.768;
%             RangeR_t=32.768;
%             XMax_t=32.768;
%             MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
%             MaxVariFactor = 0.3;                            % Max Chaos Vari Ratio
%             M_t = 0.1;                                      % Chaos Slides Control Factor
%             N_t = 0.3;                                      % Chaos Slides Control Factor
%             offset_t = -M_t;                                   % Offset between beta and chaos.
%             A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
%             alpha_t = 0.01;                               % Convergence Judgement factor
%             beta_select = 1;                                % QPSO Beta factor mode.
%             chaos_range = [0.5 0.9];                        % Chaos range
%             sig_t = 0.2;                                    % Convergence Judgement factor
%             b1 = 2;
%             b2 = 6;
%             
%             Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
%             Ans_qpso = zeros(2,ExcuteTimes_t);
%             Ans_spso = zeros(2,ExcuteTimes_t);
%             for z = 1:ExcuteTimes_t
%                 QPSO_Dynamic_Chaos
%                 QPSO_Official
%                 spso
%                 MakeDocs
%             end
%             MakeGeneralDocs
%             
%         end
%     end
% end


%% Schaffer's F6 Test   Very Good

for aa = 1:2
    
    for cc = 1:3
        
        clearvars -except aa bb cc
        Opt.Method = 'SHA-256';
        Opt.Format = 'HEX';
        Opt.Input = 'file';
        ExcuteTimes_t = 50;
        func_test = 'f6';
        func_name = 'Schaffer"s F6';
        ParticleNum_t=20*cc;             % Number of Particles
        MaxIteration_t=500+aa*500;          % Max Iterations
        dim_t=2;                     % Dimension
        RangeL_t=-100;
        RangeR_t=100;
        XMax_t=100;
        MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
        MaxVariFactor = 0.2;                            % Max Chaos Vari Ratio
        M_t = 0.05;                                      % Chaos Slides Control Factor
        N_t = 0.1;                                      % Chaos Slides Control Factor
        offset_t = -M_t;                                   % Offset between beta and chaos.
        A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
        alpha_t = 0.1;                               % Convergence Judgement factor
        beta_select = 1;                                % QPSO Beta factor mode.
        chaos_range = [0.5 0.9];                        % Chaos range
        sig_t = 10;                                    % Convergence Judgement factor
        b1 = 1;
        b2 = 4;
        c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
    end
    
end


%% Schwefel Test   Very Good

for aa = 1:2
    for bb = 1:3
        for cc = 1:3
            
            clearvars -except aa bb cc
            Opt.Method = 'SHA-256';
            Opt.Format = 'HEX';
            Opt.Input = 'file';
            ExcuteTimes_t = 50;
            func_test = 'f7';
            func_name = 'Schwefel';
            ParticleNum_t=20*cc;             % Number of Particles
            MaxIteration_t=500+aa*500;          % Max Iterations
            dim_t=10*bb+10;                     % Dimension
            RangeL_t=-500;
            RangeR_t=500;
            XMax_t=500;
            MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
            MaxVariFactor = 0.2;                            % Max Chaos Vari Ratio
            M_t = 0.1;                                      % Chaos Slides Control Factor
            N_t = 0.3;                                      % Chaos Slides Control Factor
            offset_t = -M_t;                                   % Offset between beta and chaos.
            A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
            alpha_t = 0.1;                               % Convergence Judgement factor
            beta_select = 1;                                % QPSO Beta factor mode.
            chaos_range = [0.5 0.9];                        % Chaos range
            sig_t = 10;                                    % Convergence Judgement factor
            b1 = 2;
            b2 = 5;
            
            c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
            
        end
    end
end


%% Salomon Test
% 
% for aa = 1:2
%     for bb = 1:3
%         for cc = 1:3
%             
%             clearvars -except aa bb cc
%             Opt.Method = 'SHA-256';
%             Opt.Format = 'HEX';
%             Opt.Input = 'file';
%             ExcuteTimes_t = 50;
%             func_test = 'f8';
%             func_name = 'Salomon';
%             ParticleNum_t=20*cc;             % Number of Particles
%             MaxIteration_t=500+aa*500;          % Max Iterations
%             dim_t=10*bb+10;                     % Dimension
%             RangeL_t=-100;
%             RangeR_t=100;
%             XMax_t=100;
%             MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
%             MaxVariFactor = 0.2;                            % Max Chaos Vari Ratio
%             M_t = 0.1;                                      % Chaos Slides Control Factor
%             N_t = 0.3;                                      % Chaos Slides Control Factor
%             offset_t = -M_t;                                   % Offset between beta and chaos.
%             A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
%             alpha_t = 0.1;                               % Convergence Judgement factor
%             beta_select = 1;                                % QPSO Beta factor mode.
%             chaos_range = [0.5 0.9];                        % Chaos range
%             sig_t = 10;                                    % Convergence Judgement factor
%             b1 = 1;
%             b2 = 4;
%             
%             Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
%             Ans_qpso = zeros(2,ExcuteTimes_t);
%             Ans_spso = zeros(2,ExcuteTimes_t);
%             for z = 1:ExcuteTimes_t
%                 QPSO_Dynamic_Chaos
%                 QPSO_Official
%                 spso
%                 MakeDocs
%             end
%             MakeGeneralDocs
%             
%         end
%     end
% end

%% Shekel Test

for aa = 1:2

        for cc = 1:3
            
            clearvars -except aa bb cc
            Opt.Method = 'SHA-256';
            Opt.Format = 'HEX';
            Opt.Input = 'file';
            ExcuteTimes_t = 50;
            func_test = 'f9';
            func_name = 'Shekel';
            ParticleNum_t=20*cc;             % Number of Particles
            MaxIteration_t=500+aa*500;          % Max Iterations
            dim_t=4;                     % Dimension
            RangeL_t=0;
            RangeR_t=10;
            XMax_t=10;
            MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
            MaxVariFactor = 0.2;                            % Max Chaos Vari Ratio
            M_t = 0.1;                                      % Chaos Slides Control Factor
            N_t = 0.3;                                      % Chaos Slides Control Factor
            offset_t = -M_t;                                   % Offset between beta and chaos.
            A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
            alpha_t = 0.1;                               % Convergence Judgement factor
            beta_select = 1;                                % QPSO Beta factor mode.
            chaos_range = [0.5 0.9];                        % Chaos range
            sig_t = 10;                                      % Convergence Judgement factor
            b1 = 2;
            b2 = 5;
            c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
        end

end


%% Zakharov Test

for aa = 1:2
    for bb = 1:3
        for cc = 1:3
            
            clearvars -except aa bb cc
            Opt.Method = 'SHA-256';
            Opt.Format = 'HEX';
            Opt.Input = 'file';
            ExcuteTimes_t = 50;
            func_test = 'f10';
            func_name = 'Zakharov';
            ParticleNum_t=20*cc;             % Number of Particles
            MaxIteration_t=500+aa*500;          % Max Iterations
            dim_t=10*bb+10;                     % Dimension
            RangeL_t=-5;
            RangeR_t=10;
            XMax_t=10;
            MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
            MaxVariFactor = 0.2;                            % Max Chaos Vari Ratio
            M_t = 0.1;                                      % Chaos Slides Control Factor
            N_t = 0.3;                                      % Chaos Slides Control Factor
            offset_t = -M_t;                                   % Offset between beta and chaos.
            A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
            alpha_t = 0.1;                               % Convergence Judgement factor
            beta_select = 1;                                % QPSO Beta factor mode.
            chaos_range = [0.5 0.9];                        % Chaos range
            sig_t = 10;                                      % Convergence Judgement factor
            b1 = 1;
            b2 = 4;
            c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
        end
    end
end


%% DeJong N.4 No noise Test

for aa = 1:2

        for cc = 1:3
            
            clearvars -except aa bb cc
            Opt.Method = 'SHA-256';
            Opt.Format = 'HEX';
            Opt.Input = 'file';
            ExcuteTimes_t = 50;
            func_test = 'f11';
            func_name = 'DeJong N.4 No Noise';
            ParticleNum_t=20*cc;             % Number of Particles
            MaxIteration_t=500+aa*500;          % Max Iterations
            dim_t=4;                     % Dimension
            RangeL_t=-10;
            RangeR_t=10;
            XMax_t=10;
            MaxChaosIteration_t = 50; %Iteration_t = round(1/20 * MaxIteration_t);                        % Max Chaos Iteration = MaxChaosIteration * ParticleNum
            MaxVariFactor = 0.2;                            % Max Chaos Vari Ratio
            M_t = 0.1;                                      % Chaos Slides Control Factor
            N_t = 0.3;                                      % Chaos Slides Control Factor
            offset_t = -M_t;                                   % Offset between beta and chaos.
            A_t = round(ParticleNum_t * MaxVariFactor);     % Chaos Slides Control Factor
            alpha_t = 0.1;                               % Convergence Judgement factor
            beta_select = 1;                                % QPSO Beta factor mode.
            chaos_range = [0.5 0.9];                        % Chaos range
            sig_t = 10;                                      % Convergence Judgement factor
            b1 = 1;
            b2 = 5;
            c1 = 1;
            c2 = 1;
            Ans_dcqpso = zeros(2,ExcuteTimes_t);      % [optima;eclapsed_time];
            Ans_qpso = zeros(2,ExcuteTimes_t);
            Ans_spso = zeros(2,ExcuteTimes_t);
            Curve_All_dcqpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_qpso = zeros(ExcuteTimes_t, MaxIteration_t);
            Curve_All_spso = zeros(ExcuteTimes_t, MaxIteration_t);
            for z = 1:ExcuteTimes_t
                QPSO_Dynamic_Chaos
                QPSO_Official
                spso
                MakeDocs
            end
            MakeGeneralDocs
            
        end

end



%% Shutdown
system('shutdown ¨Cs -t 60');
