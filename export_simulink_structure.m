%% ======================================================================
% [FILE METADATA]
% - File: export_simulink_structure.m
% - Target Environment: MATLAB R2022a or newer
% - Description: 시뮬링크 모델(BuckConverter)의 구조와 각 블록의 상세 설정을 텍스트 파일로 추출/저장
% ======================================================================

function export_simulink_structure()
    model_name = 'BuckConverter';
    output_filename = 'BuckConverter_Structure.txt';
    
    fprintf('=== [시작] 시뮬링크 구조 분석 기동 ===\n');
    
    % 1. 모델 존재 여부 확인 및 로드
    if exist(model_name, 'file') ~= 4
        error('에러: 모델 파일 [%s.slx]이 현재 디렉토리에 존재하지 않습니다.', model_name);
    end
    
    fprintf('- 시뮬링크 모델 [%s.slx] 로딩 중...\n', model_name);
    load_system(model_name);
    
    % 2. 텍스트 파일 개방
    fid = fopen(output_filename, 'w', 'n', 'utf-8');
    if fid == -1
        error('에러: 출력 파일 [%s]을 생성할 수 없습니다.', output_filename);
    end
    
    % 파일 헤더 기록
    fprintf(fid, '======================================================================\n');
    fprintf(fid, ' SIMULINK MODEL STRUCTURE EXPORT: %s\n', model_name);
    fprintf(fid, ' Export Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '======================================================================\n\n');
    
    % 3. 전체 블록 검색 (Subsystems, Simscape blocks 포함)
    fprintf('- 모델 내부의 모든 블록 탐색 중...\n');
    blocks = find_system(model_name, 'LookUnderMasks', 'all', 'FollowLinks', 'on', 'Type', 'block');
    N_blocks = length(blocks);
    fprintf('- 검출된 총 블록 개수: %d개\n', N_blocks);
    
    fprintf(fid, '■ 검출된 총 블록 개수: %d개\n\n', N_blocks);
    
    % 4. 각 블록에 대한 세부 설정 및 파라미터 추출
    for i = 1:N_blocks
        block_path = blocks{i};
        block_name = get_param(block_path, 'Name');
        block_type = '';
        try
            block_type = get_param(block_path, 'BlockType');
        catch
            % 일부 Simscape 또는 특수 Subsystem은 BlockType이 없을 수 있음
            block_type = 'SubSystem / Simscape / Special Block';
        end
        
        fprintf(fid, '----------------------------------------------------------------------\n');
        fprintf(fid, '[Block %d] Path: %s\n', i, block_path);
        fprintf(fid, '  Name      : %s\n', block_name);
        fprintf(fid, '  Type      : %s\n', block_type);
        
        % 블록 타입별 중요 파라미터 추출
        fprintf(fid, '  Parameters:\n');
        
        % 모든 다이얼로그 파라미터 목록 가져오기
        try
            diag_params = get_param(block_path, 'DialogParameters');
            if ~isempty(diag_params)
                p_names = fieldnames(diag_params);
                
                % 블록 타입별로 관심 있는 주요 파라미터만 출력 (파일 크기 관리용)
                switch block_type
                    case 'Gain'
                        print_param(fid, block_path, 'Gain');
                    case 'DiscreteTransferFcn'
                        print_param(fid, block_path, 'Numerator');
                        print_param(fid, block_path, 'Denominator');
                        print_param(fid, block_path, 'SampleTime');
                    case 'DiscreteIntegrator'
                        print_param(fid, block_path, 'gainval');
                        print_param(fid, block_path, 'InitialCondition');
                    case 'Constant'
                        print_param(fid, block_path, 'Value');
                    case 'Saturate'
                        print_param(fid, block_path, 'LowerLimit');
                        print_param(fid, block_path, 'UpperLimit');
                    case 'ToWorkspace'
                        print_param(fid, block_path, 'VariableName');
                        print_param(fid, block_path, 'MaxDataPoints');
                        print_param(fid, block_path, 'SaveFormat');
                    case 'FromWorkspace'
                        print_param(fid, block_path, 'VariableName');
                    case 'ZeroOrderHold'
                        print_param(fid, block_path, 'SampleTime');
                    otherwise
                        % Simscape 블록 등은 Dialog Parameters 전체 출력 시도
                        for p_idx = 1:length(p_names)
                            p_name = p_names{p_idx};
                            % 임시로 관심 매개변수 필터링 (물리값, 게인, 저항, 인덕턴스, 커패시턴스 등)
                            if any(strcmpi(p_name, {'R_closed', 'G_open', 'Vf', 'Ron', 'Lmin', 'g', 'Cmin', 'r', 'Rmin', 'Gain', 'VariableName', 'Numerator', 'Denominator', 'SampleTime', 'LowerLimit', 'UpperLimit', 'P', 'I', 'D'}))
                                print_param(fid, block_path, p_name);
                            end
                        end
                end
            end
        catch
        end
        
        % 포트 연결 상태 기록 (PortConnectivity)
        try
            conn = get_param(block_path, 'PortConnectivity');
            fprintf(fid, '  Connectivity:\n');
            for c_idx = 1:length(conn)
                c = conn(c_idx);
                fprintf(fid, '    - Port [%s]: Type = %s', c.Type, get_port_type_str(c));
                if ~isempty(c.DstBlock)
                    fprintf(fid, ' -> Connected to Blocks: ');
                    for d_idx = 1:length(c.DstBlock)
                        try
                            dst_name = get_param(c.DstBlock(d_idx), 'Name');
                            fprintf(fid, '[%s] ', dst_name);
                        catch
                        end
                    end
                end
                if ~isempty(c.SrcBlock)
                    fprintf(fid, ' <- Sourced from Blocks: ');
                    for s_idx = 1:length(c.SrcBlock)
                        try
                            src_name = get_param(c.SrcBlock(s_idx), 'Name');
                            fprintf(fid, '[%s] ', src_name);
                        catch
                        end
                    end
                end
                fprintf(fid, '\n');
            end
        catch
        end
        
        fprintf(fid, '\n');
    end
    
    fclose(fid);
    fprintf('=== [성공] 구조 분석 완료. 출력 파일: %s ===\n\n', output_filename);
end

% 매개변수 출력 서브 함수
function print_param(fid, block_path, param_name)
    try
        val = get_param(block_path, param_name);
        % 행렬 형태는 한 줄로 깔끔하게 정리
        val_str = strtrim(regexprep(num2str(val), '\s+', ' '));
        fprintf(fid, '      %-15s : %s\n', param_name, val_str);
    catch
    end
end

% 포트 타입 문자열 변환 서브 함수
function str = get_port_type_str(c)
    if isfield(c, 'Type')
        str = c.Type;
    else
        str = 'Unknown';
    end
end
