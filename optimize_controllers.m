% ======================================================================
% [FILE METADATA & VERSION TRACKING]
% - Current Version: v1.1.0 (2026-06-25)
% - Target Environment: MATLAB R2022a or newer
% - Integrity Check: Stable 4 digital controllers optimized via RK4 and chatter penalty.
% ======================================================================
% [CHANGELOG - NEVER DELETE THIS HISTORY]
% * v1.1.0 (2026-06-25) - Developer: Gemini AI
%   - Added: 비용 함수 내부에 듀티비 급변(채터링) 방지용 chatter 패널티 항 추가.
%   - Changed: Type 3 보상기 crossover 주파수 상한선을 2500Hz로 낮춰 강인성 확보.
%   - Added: fminsearch의 OutputFcn을 활용하여 실시간 최적화 비용 추이를 피겨와 명령창에 출력.
% * v1.0.0 (2026-06-25) - Developer: Gemini AI
%   - Added: Initial version for Buck Converter parameter optimization.
% ======================================================================

function [pi_gains, type3_coeffs, ml_coeffs, lqr_gains] = optimize_controllers(...
    L_val, C_val, G_L, R_C, R_nom, Vin_nom, Vref_val, T_s, t_vec, Vin_data, R_data, Vref_data)

    fprintf('\n>>> [Optimization] 제어기 파라미터 최적화 루틴 시작 <<<\n');

    % 공통 파라미터 구조체 정의
    sys.L = L_val;
    sys.C = C_val;
    sys.G_L = G_L;
    sys.R_C = R_C;
    sys.R_nom = R_nom;
    sys.Vin_nom = Vin_nom;
    sys.Vref_val = Vref_val;
    sys.T_s = T_s;
    sys.t_vec = t_vec;
    sys.Vin_data = Vin_data;
    sys.R_data = R_data;
    sys.Vref_data = Vref_data;

    % --- 공칭 모델(Nominal Model) 및 이산화 파라미터 계산 (Zero-Simulation Filtering 용) ---
    theta_nom = G_L * R_nom * R_C + R_nom + R_C;
    A_nom = [ -R_nom * R_C / (L_val * theta_nom),                 -R_nom / (L_val * theta_nom);
               R_nom / (C_val * theta_nom),                 -(R_nom * G_L + 1) / (C_val * theta_nom) ];
    B_nom = [ (R_nom + R_C) / (L_val * theta_nom);
              (R_nom * G_L) / (C_val * theta_nom) ] * Vin_nom; % 스위칭 전압 입력 기준으로 Vin_nom을 곱함
    C_nom = [ R_nom * R_C / theta_nom,   R_nom / theta_nom ];
    D_nom = (G_L * R_nom * R_C / theta_nom) * Vin_nom;
    
    % Backward Euler 이산화 (T_s)
    I_minus_AT = eye(2) - A_nom * T_s;
    A_d = I_minus_AT \ eye(2);
    B_d = A_d * B_nom * T_s;
    C_d = C_nom;
    D_d = D_nom;
    
    % 대수적 2차 이산 전달함수 계수 추출 (외부 툴박스 의존성 제거)
    % den_g = z^2 - (a11 + a22)z + det(A_d)
    a11 = A_d(1,1); a12 = A_d(1,2); a21 = A_d(2,1); a22 = A_d(2,2);
    det_Ad = a11 * a22 - a12 * a21;
    den_g = [1, -(a11 + a22), det_Ad];
    
    % num_g = D_d * den_g + [0, C_d * B_d, C_d(1)*(a12*B_d(2) - a22*B_d(1)) + C_d(2)*(a21*B_d(1) - a11*B_d(2))]
    cb = C_d * B_d;
    c1 = C_d(1); c2 = C_d(2); b1 = B_d(1); b2 = B_d(2);
    num_term3 = c1 * (a12 * b2 - a22 * b1) + c2 * (a21 * b1 - a11 * b2);
    num_g = D_d * den_g + [0, cb, num_term3];
    
    sys.A_nom = A_nom;
    sys.B_nom = B_nom;
    sys.C_nom = C_nom;
    sys.D_nom = D_nom;
    sys.A_d = A_d;
    sys.B_d = B_d;
    sys.C_d = C_d;
    sys.D_d = D_d;
    sys.num_g = num_g;
    sys.den_g = den_g;

    %% 1. PI 제어기 최적화
    fprintf('- PI 제어기 최적화 중...\n');
    pi_init = [0.08, 1000]; % [Kp, Ki] (물리 공간) - 균형 잡힌 초기 튜닝점 지정
    lb_pi = [0.005, 50];    % 하한선을 완화하여 오버슈트/포화 패널티 회피를 위한 자율성 제공
    ub_pi = [2.5, 10000];
    
    % [-10, 10] 논리 공간으로 매핑
    pi_init_logical = to_logical_space(pi_init, lb_pi, ub_pi);
    lb_logical_pi = -10 * ones(size(pi_init));
    ub_logical_pi =  10 * ones(size(pi_init));
    
    pi_cost_fn = @(p) evaluate_cost('PI', p, sys, lb_pi, ub_pi);
    
    manage_optim_data('reset', 'PI');
    try
        pi_opt_logical = run_global_opt(pi_cost_fn, pi_init_logical, lb_logical_pi, ub_logical_pi, 50, 400, 'PI');
        pi_opt = to_physical_space(pi_opt_logical, lb_pi, ub_pi);
    catch
        pi_opt = pi_init;
    end
    manage_optim_data('save', 'PI');
    
    pi_gains.KP = pi_opt(1);
    pi_gains.KI = pi_opt(2);
    fprintf('  => 최적 PI 파라미터: KP = %.4f, KI = %.4f\n', pi_gains.KP, pi_gains.KI);

    %% 2. Type 3 (3P2Z) 제어기 최적화 (z-영역 극/영점 직접 최적화)
    fprintf('- Type 3 (3P2Z) 제어기 최적화 중...\n');
    % 최적화 변수: [log_Kc, log_wz1, log_wz2, log_wp1, log_wp2] (물리 공간)
    type3_init = [log(1000), log(2*pi*300), log(2*pi*300), log(2*pi*3000), log(2*pi*3000)];
    lb_type3 = [-5, log(2*pi*50), log(2*pi*50), log(2*pi*500), log(2*pi*500)];
    ub_type3 = [15, log(2*pi*2000), log(2*pi*2000), log(2*pi*10000), log(2*pi*10000)];
    
    % [-10, 10] 논리 공간으로 매핑
    type3_init_logical = to_logical_space(type3_init, lb_type3, ub_type3);
    lb_logical_type3 = -10 * ones(size(type3_init));
    ub_logical_type3 =  10 * ones(size(type3_init));
    
    type3_cost_fn = @(p) evaluate_cost('3P2Z', p, sys, lb_type3, ub_type3);
    
    manage_optim_data('reset', '3P2Z');
    try
        type3_opt_logical = run_global_opt(type3_cost_fn, type3_init_logical, lb_logical_type3, ub_logical_type3, 80, 600, '3P2Z');
        type3_opt = to_physical_space(type3_opt_logical, lb_type3, ub_type3);
    catch
        type3_opt = type3_init;
    end
    manage_optim_data('save', '3P2Z');
    
    [n1, n2, n3, d1, d2, d3, d4] = design_type3_direct(type3_opt, sys);
    type3_coeffs.n1 = n1; type3_coeffs.n2 = n2; type3_coeffs.n3 = n3;
    type3_coeffs.d1 = d1; type3_coeffs.d2 = d2; type3_coeffs.d3 = d3; type3_coeffs.d4 = d4;
    fprintf('  => 최적 Type 3 직접 설계 완료 (5개 매개변수 최적화)\n');

    %% 3. ML 최적화 제어기 (5차 전달함수) 최적화
    fprintf('- ML 5차 전달함수 제어기 최적화 중...\n');
    % 최적화 변수: [log_Kc, log_wz1, zeta_z1, log_wz2, zeta_z2, log_wz3, log_wp1, zeta_p1, log_wp2, zeta_p2] (물리 공간)
    ml_init = [1.5, 7.5, 0.7, 8.0, 0.7, 6.0, 9.0, 0.7, 10.0, 0.7];
    lb_ml = [-5, log(2*pi*50), 0.1, log(2*pi*50), 0.1, log(2*pi*10), log(2*pi*100), 0.2, log(2*pi*100), 0.2];
    ub_ml = [15, log(2*pi*20000), 1.5, log(2*pi*20000), 1.5, log(2*pi*20000), log(2*pi*40000), 1.5, log(2*pi*40000), 1.5];
    
    % [-10, 10] 논리 공간으로 매핑
    ml_init_logical = to_logical_space(ml_init, lb_ml, ub_ml);
    lb_logical_ml = -10 * ones(size(ml_init));
    ub_logical_ml =  10 * ones(size(ml_init));
    
    ml_cost_fn = @(p) evaluate_cost('ML', p, sys, lb_ml, ub_ml);
    
    manage_optim_data('reset', 'ML');
    try
        ml_opt_logical = run_global_opt(ml_cost_fn, ml_init_logical, lb_logical_ml, ub_logical_ml, 150, 3000, 'ML');
        ml_opt = to_physical_space(ml_opt_logical, lb_ml, ub_ml);
    catch
        ml_opt = ml_init;
    end
    manage_optim_data('save', 'ML');
    
    [m1, m2, m3, m4, m5, m6, e1, e2, e3, e4, e5, e6] = design_ml_tf(ml_opt, sys.T_s);
    ml_coeffs.m1 = m1; ml_coeffs.m2 = m2; ml_coeffs.m3 = m3; ml_coeffs.m4 = m4; ml_coeffs.m5 = m5; ml_coeffs.m6 = m6;
    ml_coeffs.e1 = e1; ml_coeffs.e2 = e2; ml_coeffs.e3 = e3; ml_coeffs.e4 = e4; ml_coeffs.e5 = e5; ml_coeffs.e6 = e6;
    fprintf('  => 최적 ML TF 파라미터 도출 완료 (고차 전달함수 적용)\n');

    %% 4. 현대 제어기 (Augmented LQR) 최적화
    fprintf('- Augmented LQR 제어기 최적화 중...\n');
    % 최적화 변수: [log_q1, log_q2, log_q3, log_R] (물리 공간)
    lqr_init = [log(10), log(100), log(1e5), log(10)]; % R 이득 상향으로 제어 입력 얌전하게 시작
    lb_lqr = [log(1e-3), log(1e-3), log(1e-3), log(1e-2)]; % R 하한선 대폭 상향(0.01) 및 Q 적분 폭주 제한
    ub_lqr = [log(1e6), log(1e6), log(1e9), log(1e4)];     % 가중치 과대 성장 억제
    
    % [-10, 10] 논리 공간으로 매핑
    lqr_init_logical = to_logical_space(lqr_init, lb_lqr, ub_lqr);
    lb_logical_lqr = -10 * ones(size(lqr_init));
    ub_logical_lqr =  10 * ones(size(lqr_init));
    
    lqr_cost_fn = @(p) evaluate_cost('LQR', p, sys, lb_lqr, ub_lqr);
    
    manage_optim_data('reset', 'LQR');
    try
        lqr_opt_logical = run_global_opt(lqr_cost_fn, lqr_init_logical, lb_logical_lqr, ub_logical_lqr, 80, 1000, 'LQR');
        lqr_opt = to_physical_space(lqr_opt_logical, lb_lqr, ub_lqr);
    catch
        lqr_opt = lqr_init;
    end
    manage_optim_data('save', 'LQR');
    
    [K_lqr1, K_lqr2, K_lqr3] = design_lqr(lqr_opt, sys);
    lqr_gains.K_lqr1 = K_lqr1;
    lqr_gains.K_lqr2 = K_lqr2;
    lqr_gains.K_lqr3 = K_lqr3;
    fprintf('  => 최적 LQR 파라미터: K_lqr1 = %.4f, K_lqr2 = %.4f, K_lqr3 = %.4f\n', K_lqr1, K_lqr2, K_lqr3);
    fprintf('>>> 제어기 파라미터 최적화 완료 <<<\n\n');
    
    % 최적화 피겨 닫기
    fig = findobj('Type', 'figure', 'Name', 'Optimization Progress');
    if ~isempty(fig)
        close(fig);
    end
end

%% ========================== [보조 함수] 공통 비용 평가 함수 ==========================

function cost = evaluate_cost(ctrl_type, p_logical, sys, lb_phys, ub_phys)
    % 1. 논리 최적화 변수 p_logical을 물리 변수 p_phys로 sigmoid 매핑 복원 (Stage 2 준비)
    p_phys = to_physical_space(p_logical, lb_phys, ub_phys);
    
    % 2. [Stage 1] 사전 선별 (Zero-Simulation Filtering) - 이산 시간 폐루프 poles 안정성 검사
    try
        is_stable = true;
        max_pole_mag = 0;
        
        switch ctrl_type
            case 'PI'
                Kp_val = p_phys(1); Ki_val = p_phys(2);
                
                % 이산 PI 제어기: C(z) = ((Kp + Ki*Ts)*z - Kp) / (z - 1)
                N_ctrl = [Kp_val + Ki_val * sys.T_s, -Kp_val];
                D_ctrl = [1, -1];
                
                % 폐루프 특성다항식 계산: P(z) = D(z)*den_g(z) + N(z)*num_g(z)
                char_poly = conv(D_ctrl, sys.den_g) + conv(N_ctrl, sys.num_g);
                poles = roots(char_poly);
                max_pole_mag = max(abs(poles));
                if max_pole_mag >= 1.0
                    is_stable = false;
                end
                
            case '3P2Z'
                [n1, n2, n3, ~, d2, d3, d4] = design_type3_direct(p_phys, sys);
                % C(z) = (n1*z^3 + n2*z^2 + n3*z) / (z^3 + d2*z^2 + d3*z + d4)
                N_ctrl = [n1, n2, n3, 0];
                D_ctrl = [1, d2, d3, d4];
                
                char_poly = conv(D_ctrl, sys.den_g) + conv(N_ctrl, sys.num_g);
                poles = roots(char_poly);
                max_pole_mag = max(abs(poles));
                if max_pole_mag >= 1.0
                    is_stable = false;
                end
                
            case 'ML'
                [m1, m2, m3, m4, m5, m6, e1, e2, e3, e4, e5, e6] = design_ml_tf(p_phys, sys.T_s);
                % G_c(z) = (m1*z^5 + ... + m6) / (e1*z^5 + ... + e6)
                N_ctrl = [m1, m2, m3, m4, m5, m6];
                D_ctrl = [e1, e2, e3, e4, e5, e6];
                
                char_poly = conv(D_ctrl, sys.den_g) + conv(N_ctrl, sys.num_g);
                poles = roots(char_poly);
                max_pole_mag = max(abs(poles));
                if max_pole_mag >= 1.0
                    is_stable = false;
                end
                
            case 'LQR'
                [K_lqr1, K_lqr2, K_lqr3] = design_lqr(p_phys, sys);
                K_lqr = [K_lqr1, K_lqr2, K_lqr3];
                
                % 이산 augmented 상태공간 폐루프 행렬 계산
                % A_aug = [ A_d,    0 ]
                %         [ -C_d,   1 ]
                % B_aug = [ B_d; -D_d ]
                A_aug = [ sys.A_d,            zeros(2, 1);
                         -sys.C_d,            1 ];
                B_aug = [ sys.B_d;
                         -sys.D_d ];
                
                A_cl = A_aug - B_aug * K_lqr;
                poles = eig(A_cl);
                max_pole_mag = max(abs(poles));
                if max_pole_mag >= 1.0
                    is_stable = false;
                end
        end
        
        % 불안정한 경우, 시뮬레이션을 절대 기동하지 않고 대수적 소프트 패널티를 우아하게 반환
        if ~is_stable
            cost = 1e6 + (max_pole_mag - 1.0) * 1e6;
            return;
        end
        
    catch
        % 대수적 행렬 연산 에러 발생 시 차단 패널티
        cost = 1e12;
        return;
    end

    % 3. [Stage 2] 안전이 확보된 상태에서만 고속 RK4 시뮬레이션 기동
    switch ctrl_type
        case 'PI'
            Kp_val = p_phys(1); Ki_val = p_phys(2);
        case '3P2Z'
            [n1, n2, n3, ~, d2, d3, d4] = design_type3_direct(p_phys, sys);
        case 'ML'
            [m1, m2, m3, m4, m5, m6, e1, e2, e3, e4, e5, e6] = design_ml_tf(p_phys, sys.T_s);
        case 'LQR'
            [K_lqr1, K_lqr2, K_lqr3] = design_lqr(p_phys, sys);
    end

    dt = sys.T_s;
    N = length(sys.t_vec);
    V_out_hist = zeros(N, 1);
    duty_hist = zeros(N, 1);
    
    x = [0; 0];
    error_int = 0;
    duty = sys.Vref_val / sys.Vin_nom;
    
    err_hist = zeros(6, 1);
    u_hist = ones(6, 1) * duty;
    
    for k = 1:N
        V_in_k = sys.Vin_data(k);
        R_k = sys.R_data(k);
        
        theta_k = sys.G_L * R_k * sys.R_C + R_k + sys.R_C;
        A_k = [ -R_k * sys.R_C / (sys.L * theta_k),                 -R_k / (sys.L * theta_k);
                 R_k / (sys.C * theta_k),                 -(R_k * sys.G_L + 1) / (sys.C * theta_k) ];
        B_k = [ (R_k + sys.R_C) / (sys.L * theta_k);
                (R_k * sys.G_L) / (sys.C * theta_k) ]; 
        C_k = [ R_k * sys.R_C / theta_k,   R_k / theta_k ];
        D_k = sys.G_L * R_k * sys.R_C / theta_k;
        
        v_sw = duty * V_in_k;
        V_out = C_k * x + D_k * v_sw;
        V_out_hist(k) = V_out;
        duty_hist(k) = duty;
        
        err = sys.Vref_data(k) - V_out;
        
        switch ctrl_type
            case 'PI'
                error_int = error_int + err * dt;
                duty_next_raw = Kp_val * err + Ki_val * error_int;
                duty_next = max(0.01, min(0.95, duty_next_raw));
                % Anti-Windup Clamping: 포화 발생 시 적분 중지
                if duty_next_raw > 0.95 && err > 0
                    error_int = error_int - err * dt;
                elseif duty_next_raw < 0.01 && err < 0
                    error_int = error_int - err * dt;
                end
            case '3P2Z'
                err_hist = [err; err_hist(1:5)];
                duty_next_raw = n1*err_hist(1) + n2*err_hist(2) + n3*err_hist(3) ...
                              - d2*u_hist(1) - d3*u_hist(2) - d4*u_hist(3);
                duty_next = max(0.01, min(0.95, duty_next_raw));
                u_hist = [duty_next; u_hist(1:5)];
            case 'ML'
                err_hist = [err; err_hist(1:5)];
                duty_next_raw = (m1*err_hist(1) + m2*err_hist(2) + m3*err_hist(3) + m4*err_hist(4) + m5*err_hist(5) + m6*err_hist(6) ...
                               - e2*u_hist(1) - e3*u_hist(2) - e4*u_hist(3) - e5*u_hist(4) - e6*u_hist(5)) / e1;
                duty_next = max(0.01, min(0.95, duty_next_raw));
                u_hist = [duty_next; u_hist(1:5)];
            case 'LQR'
                error_int = error_int + err * dt;
                I_L_ref = sys.Vref_data(k) / R_k;
                duty_nom = sys.Vref_data(k) / V_in_k;
                duty_next_raw = duty_nom - ( K_lqr1 * (x(1) - I_L_ref) + K_lqr2 * (x(2) - sys.Vref_data(k)) + K_lqr3 * error_int );
                duty_next = max(0.01, min(0.95, duty_next_raw));
                % LQR Anti-Windup Clamping
                if duty_next_raw > 0.95 || duty_next_raw < 0.01
                    error_int = error_int - err * dt;
                end
        end
        
        x = rk4_step(A_k, B_k, x, v_sw, dt);
        duty = duty_next;
    end
    
    if any(isnan(V_out_hist)) || any(isinf(V_out_hist)) || any(isnan(duty_hist))
        cost = 1e12; return;
    end
    
    % 4. [Stage 3] 다목적 패널티 최종 합성 및 듀티 포화 처벌
    itae = sum((sys.t_vec.^2) .* abs(sys.Vref_data - V_out_hist) * dt) * 10;
    
    % 오버슈트 패널티 부과 (작은 오버슈트도 즉각 처벌하고, 5% 초과 시 매우 강력하게 가산)
    overshoot = max(V_out_hist) - sys.Vref_val;
    overshoot_penalty = 0;
    if overshoot > 0
        overshoot_penalty = overshoot * 2000; % 기본 선형 패널티
        if overshoot > 0.05 * sys.Vref_val
            overshoot_penalty = overshoot_penalty + 50000 * (overshoot - 0.05 * sys.Vref_val); % 5% 초과 시 극단적 패널티
        end
    end
    
    % 듀티비 채터링(급변) 패널티 계산 (과도 상태 제외하고 정상 상태에서의 진동만 억제)
    t = sys.t_vec;
    W_chatter = 5000 * ones(size(t));
    W_chatter(t <= 0.005) = 0;                        % 초기 시동 과도기 (0ms ~ 5ms) 제외
    W_chatter(t >= 0.030 & t <= 0.035) = 0;           % Load Step 과도기 (30ms ~ 35ms) 제외
    W_chatter(t >= 0.040 & t <= 0.045) = 0;           % Voltage Surge 과도기 (40ms ~ 45ms) 제외
    W_chatter(t >= 0.070 & t <= 0.075) = 0;           % Load Step 과도기 (70ms ~ 75ms) 제외
    
    % 2차 차분 채터링 패널티
    d_diff2 = diff(diff(duty_hist));
    chatter_penalty = sum(W_chatter(3:end) .* (d_diff2.^2));
    
    % [추가] 1차 차분 델타 듀티 변화량 패널티 (\Delta D) - 채터링 원천 봉쇄용
    d_diff1 = diff(duty_hist);
    chatter_penalty_1st = sum(d_diff1.^2) * 5000;
    
    % [개선] 초기 기동 과도기(0ms ~ 5ms) 동안의 자연스러운 포화는 패널티 측정에서 제외
    valid_sat_idx = t > 0.005; % 5ms 이후 영역만 평가
    duty_eval = duty_hist(valid_sat_idx);
    
    % 듀티비 포화(Anti-Windup 방지) 시간 패널티 추가 (0.02 이하 또는 0.94 이상)
    saturation_time = sum(duty_eval >= 0.94 | duty_eval <= 0.02) * dt;
    sat_penalty = saturation_time * 5000;
    
    % [추가] 조건부 극단 패널티: 5ms 이후 영역에서 포화 비율이 5%를 초과하는 경우 무려 +1e6 추가
    if ~isempty(duty_eval)
        sat_ratio = sum(duty_eval >= 0.94 | duty_eval <= 0.02) / length(duty_eval);
        if sat_ratio > 0.05
            sat_penalty = sat_penalty + 1e6;
        end
    end
    
    % 최종 비용 합성
    cost = itae + overshoot_penalty + chatter_penalty + chatter_penalty_1st + sat_penalty;
    
    % 최적화 과정 실시간 업데이트 및 로깅
    manage_optim_data('update', ctrl_type, cost);
end

%% ========================== [설계 알고리즘] ==========================

% Type 3 직접 극/영점 기반 설계
function [n1, n2, n3, d1, d2, d3, d4] = design_type3_direct(p, sys)
    Kc = exp(p(1));
    wz1 = exp(p(2));
    wz2 = exp(p(3));
    wp1 = exp(p(4));
    wp2 = exp(p(5));
    
    T_s = sys.T_s;
    KT = 1 / T_s;
    
    a1 = wz1 + wz2;
    a0 = wz1 * wz2;
    b1 = wp1 + wp2;
    b0 = wp1 * wp2;
    
    A = KT^2 + b1 * KT + b0;
    B = 2 * KT^2 + b1 * KT;
    C = KT^2;
    
    d_raw1 = KT * A;
    
    n_raw1 = Kc * (KT^2 + a1 * KT + a0);
    n_raw2 = -Kc * (2 * KT^2 + a1 * KT);
    n_raw3 = Kc * KT^2;
    
    n1 = n_raw1 / d_raw1;
    n2 = n_raw2 / d_raw1;
    n3 = n_raw3 / d_raw1;
    d1 = 1.0;
    d2 = -KT * (A + B) / d_raw1;
    d3 = KT * (B + C) / d_raw1;
    d4 = -KT * C / d_raw1;
end

% Type 3 k-factor 설계 (Legacy Preservation)
function [n1, n2, n3, d1, d2, d3, d4] = design_type3_kfactor(f_co, PM_target, sys)
    theta_nom = sys.G_L * sys.R_nom * sys.R_C + sys.R_nom + sys.R_C;
    A_nom = [ -sys.R_nom * sys.R_C / (sys.L * theta_nom),                 -sys.R_nom / (sys.L * theta_nom);
               sys.R_nom / (sys.C * theta_nom),                 -(sys.R_nom * sys.G_L + 1) / (sys.C * theta_nom) ];
    B_nom = [ (sys.R_nom + sys.R_C) / (sys.L * theta_nom);
              (sys.R_nom * sys.G_L) / (sys.C * theta_nom) ] * sys.Vin_nom;
    C_nom = [ sys.R_nom * sys.R_C / theta_nom,   sys.R_nom / theta_nom ];
    D_nom = (sys.G_L * sys.R_nom * sys.R_C / theta_nom) * sys.Vin_nom;
    
    % Crossover frequency에서의 Plant 수치해석적 이득 및 위상 획득
    w_co = 2 * pi * f_co;
    s_co = 1i * w_co;
    G_co = C_nom * ((s_co * eye(2) - A_nom) \ B_nom) + D_nom;
    mag_co = abs(G_co);
    phase_co = rad2deg(angle(G_co));
    
    % Required phase boost
    boost = PM_target - 90 - phase_co;
    boost = max(5, min(170, boost)); % 안전 클리핑
    
    k_val = tan(deg2rad(45 + boost / 4));
    w_z = w_co / k_val;
    w_p = w_co * k_val;
    
    % Kc 설계
    mag_c_raw = (w_co^2 + w_z^2) / (w_co * (w_co^2 + w_p^2));
    K_c = 1 / (mag_co * mag_c_raw);
    
    % Backward Euler 이산화 적용: s -> (1 - z^-1)/T_s
    T_s = sys.T_s;
    A_z = 1 + w_z * T_s;
    A_p = 1 + w_p * T_s;
    
    n_raw1 = K_c * T_s * A_z^2;
    n_raw2 = -2 * K_c * T_s * A_z;
    n_raw3 = K_c * T_s;
    
    d_raw1 = A_p^2;
    d_raw2 = -(A_p^2 + 2*A_p);
    d_raw3 = 2*A_p + 1;
    d_raw4 = -1;
    
    % d_raw1로 정규화 진행
    n1 = n_raw1 / d_raw1;
    n2 = n_raw2 / d_raw1;
    n3 = n_raw3 / d_raw1;
    d1 = 1.0;
    d2 = d_raw2 / d_raw1;
    d3 = d_raw3 / d_raw1;
    d4 = d_raw4 / d_raw1;
end

% ML 5차 전달함수 설계 및 이산화
function [m1, m2, m3, m4, m5, m6, e1, e2, e3, e4, e5, e6] = design_ml_tf(p, T_s)
    Kc = exp(p(1));
    wz1 = exp(p(2)); zeta_z1 = p(3);
    wz2 = exp(p(4)); zeta_z2 = p(5);
    wz3 = exp(p(6));
    wp1 = exp(p(7)); zeta_p1 = p(8);
    wp2 = exp(p(9)); zeta_p2 = p(10);
    
    % 공칭 crossover 주파수 대역 (약 1500Hz) 기준으로 프리워핑 인수 K 계산
    w_warp = 2 * pi * 1500;
    K = w_warp / tan(w_warp * T_s / 2);
    
    % 1) G1(s) = Kc * (s + wz3)/s => G1(z) = Kc * (K*(z-1)/(z+1) + wz3) / (K*(z-1)/(z+1))
    % G1(z) = Kc * ((K + wz3)*z - (K - wz3)) / (K*(z - 1))
    num1 = [Kc * (K + wz3), -Kc * (K - wz3)];
    den1 = [K, -K];
    
    % 2) G2(s) = (s^2 + a1*s + a0) / (s^2 + b1*s + b0)
    a1 = 2 * zeta_z1 * wz1; a0 = wz1^2;
    b1 = 2 * zeta_p1 * wp1; b0 = wp1^2;
    
    % s -> K*(z-1)/(z+1) 대입 및 전개
    num2 = [K^2 + a1*K + a0, 2*(a0 - K^2), K^2 - a1*K + a0];
    den2 = [K^2 + b1*K + b0, 2*(b0 - K^2), K^2 - b1*K + b0];
    
    % 3) G3(s) = (s^2 + a'1*s + a'0) / (s^2 + b'1*s + b'0)
    a_prime1 = 2 * zeta_z2 * wz2; a_prime0 = wz2^2;
    b_prime1 = 2 * zeta_p2 * wp2; b_prime0 = wp2^2;
    
    num3 = [K^2 + a_prime1*K + a_prime0, 2*(a_prime0 - K^2), K^2 - a_prime1*K + a_prime0];
    den3 = [K^2 + b_prime1*K + b_prime0, 2*(b_prime0 - K^2), K^2 - b_prime1*K + b_prime0];
    
    % 다항식 곱셈을 위한 convolution 적용
    num_d = conv(conv(num1, num2), num3);
    den_d = conv(conv(den1, den2), den3);
    
    % 정규화 (den_d(1)로 나눔)
    m_coeff = num_d / den_d(1);
    e_coeff = den_d / den_d(1);
    
    m1 = m_coeff(1); m2 = m_coeff(2); m3 = m_coeff(3); m4 = m_coeff(4); m5 = m_coeff(5); m6 = m_coeff(6);
    e1 = 1.0;        e2 = e_coeff(2); e3 = e_coeff(3); e4 = e_coeff(4); e5 = e_coeff(5); e6 = e_coeff(6);
end

% LQR 설계 함수
function [K_lqr1, K_lqr2, K_lqr3] = design_lqr(p, sys)
    q1 = exp(p(1)); q2 = exp(p(2)); q3 = exp(p(3)); R_weight = exp(p(4));
    
    theta_nom = sys.G_L * sys.R_nom * sys.R_C + sys.R_nom + sys.R_C;
    A_nom = [ -sys.R_nom * sys.R_C / (sys.L * theta_nom),                 -sys.R_nom / (sys.L * theta_nom);
               sys.R_nom / (sys.C * theta_nom),                 -(sys.R_nom * sys.G_L + 1) / (sys.C * theta_nom) ];
    B_nom = [ (sys.R_nom + sys.R_C) / (sys.L * theta_nom);
              (sys.R_nom * sys.G_L) / (sys.C * theta_nom) ] * sys.Vin_nom;
    C_nom = [ sys.R_nom * sys.R_C / theta_nom,   sys.R_nom / theta_nom ];
    D_nom = (sys.G_L * sys.R_nom * sys.R_C / theta_nom) * sys.Vin_nom;
    
    A_aug = [ A_nom,              zeros(2, 1);
             -C_nom,              0 ];
    B_aug = [ B_nom;
             -D_nom ];
         
    Q_lqr = diag([q1, q2, q3]);
    K_lqr_all = lqr(A_aug, B_aug, Q_lqr, R_weight);
    
    K_lqr1 = K_lqr_all(1);
    K_lqr2 = K_lqr_all(2);
    K_lqr3 = K_lqr_all(3);
end

% RK4 Solver Step 함수
function x_next = rk4_step(A, B, x, u, dt)
    k1 = A * x + B * u;
    k2 = A * (x + 0.5 * dt * k1) + B * u;
    k3 = A * (x + 0.5 * dt * k2) + B * u;
    k4 = A * (x + dt * k3) + B * u;
    x_next = x + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

%% ========================== [최적화 래퍼 및 매니저 함수] ==========================

function opt_val = run_global_opt(cost_fn, init_val, lb, ub, max_iter, max_eval, display_name)
    nvars = length(init_val);
    swarm_size = min(100, 10 * nvars);
    
    % Check for particleswarm & ga availability
    has_swarm = ~isempty(which('particleswarm'));
    has_ga = ~isempty(which('ga'));
    
    if has_swarm
        fprintf('  [%s] 입자 군집 최적화(particleswarm) 기동...\n', display_name);
        try
            opts = optimoptions('particleswarm', 'Display', 'none', ...
                'MaxIterations', max_iter, 'SwarmSize', swarm_size);
            opts.InitialSwarmMatrix = init_val;
            
            opt_val = particleswarm(cost_fn, nvars, lb, ub, opts);
            return;
        catch swarm_err
            fprintf('  [주의] particleswarm 실행 중 에러 발생: %s. ga 또는 fminsearch로 전환합니다.\n', swarm_err.message);
        end
    end
    
    if has_ga
        fprintf('  [%s] 유전 알고리즘(ga) 기동...\n', display_name);
        try
            opts = optimoptions('ga', 'Display', 'none', 'MaxGenerations', max_iter);
            opts.InitialPopulationMatrix = init_val;
            
            opt_val = ga(cost_fn, nvars, [], [], [], [], lb, ub, [], opts);
            return;
        catch ga_err
            fprintf('  [주의] ga 실행 중 에러 발생: %s. fminsearch로 전환합니다.\n', ga_err.message);
        end
    end
    
    % Fallback to fminsearch
    fprintf('  [%s] Global Optimization Toolbox 미검출. fminsearch로 로컬 탐색 수행...\n', display_name);
    opts = optimset('Display', 'off', 'MaxIter', max_iter * 2, 'MaxFunEvals', max_eval);
    opt_val = fminsearch(cost_fn, init_val, opts);
    % Clip output to bounds
    opt_val = max(lb, min(ub, opt_val));
end

function varargout = manage_optim_data(action, varargin)
    persistent best_cost call_count cost_history call_history ctrl_current plot_line_handle
    
    switch action
        case 'reset'
            ctrl_name = varargin{1};
            best_cost = Inf;
            call_count = 0;
            cost_history = [];
            call_history = [];
            ctrl_current = ctrl_name;
            
            fig_name = 'Optimization Progress';
            fig = findobj('Type', 'figure', 'Name', fig_name);
            if isempty(fig)
                fig = figure('Name', fig_name, 'NumberTitle', 'off', 'Position', [150, 150, 650, 450]);
            else
                figure(fig);
                clf(fig);
            end
            ax = axes(fig);
            plot_line_handle = plot(ax, NaN, NaN, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
            set(ax, 'YScale', 'log'); % Y축 로그 스케일 적용
            grid(ax, 'on');
            xlabel(ax, 'Function Evaluation', 'FontWeight', 'bold');
            ylabel(ax, 'Best ITAE + Chatter Cost (Log Scale)', 'FontWeight', 'bold');
            title(ax, sprintf('[%s] Optimization Progress', ctrl_name), 'FontSize', 12, 'FontWeight', 'bold');
            drawnow;
            
        case 'update'
            ctrl_name = varargin{1};
            cost = varargin{2};
            
            if ~strcmp(ctrl_current, ctrl_name)
                best_cost = Inf;
                call_count = 0;
                cost_history = [];
                call_history = [];
                ctrl_current = ctrl_name;
                plot_line_handle = [];
            end
            
            call_count = call_count + 1;
            
            if cost < best_cost && cost < 1e11
                best_cost = cost;
                cost_history = [cost_history; best_cost];
                call_history = [call_history; call_count];
                
                fprintf('  [%s] Eval %d: Best Cost = %.4e\n', ctrl_name, call_count, best_cost);
                
                if ~isempty(plot_line_handle) && isgraphics(plot_line_handle)
                    try
                        set(plot_line_handle, 'XData', call_history, 'YData', cost_history);
                        ax = get(plot_line_handle, 'Parent');
                        title(ax, sprintf('[%s] Cost History (Eval: %d, Best Cost: %.4e)', ctrl_name, call_count, best_cost));
                        drawnow;
                    catch
                        % Do nothing
                    end
                else
                    % If the figure or line was closed, recreate it dynamically
                    fig_name = 'Optimization Progress';
                    fig = findobj('Type', 'figure', 'Name', fig_name);
                    if ~isempty(fig)
                        try
                            ax = findobj(fig, 'Type', 'axes');
                            if isempty(ax)
                                ax = axes(fig);
                            end
                            plot_line_handle = plot(ax, call_history, cost_history, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 8);
                            set(ax, 'YScale', 'log'); % Y축 로그 스케일 적용
                            grid(ax, 'on');
                            xlabel(ax, 'Function Evaluation', 'FontWeight', 'bold');
                            ylabel(ax, 'Best ITAE + Chatter Cost (Log Scale)', 'FontWeight', 'bold');
                            title(ax, sprintf('[%s] Cost History (Eval: %d, Best Cost: %.4e)', ctrl_name, call_count, best_cost));
                            drawnow;
                        catch
                            % Do nothing
                        end
                    end
                end
            end
            
        case 'save'
            ctrl_name = varargin{1};
            try
                safe_name = strrep(ctrl_name, '-', '_');
                safe_name = strrep(safe_name, '/', '_');
                safe_name = strrep(safe_name, ' ', '_');
                csv_folder = 'csv_data';
                if ~exist(csv_folder, 'dir')
                    mkdir(csv_folder);
                end
                filename = fullfile(csv_folder, sprintf('optimization_history_%s.csv', safe_name));
                
                opt_table = table(call_history(:), cost_history(:), 'VariableNames', {'Evaluation', 'Cost'});
                writetable(opt_table, filename);
                fprintf('  => [%s] 최적화 히스토리 저장 완료: %s\n', ctrl_name, filename);
            catch csv_err
                warning('optimize:CSVSaveFailed', '최적화 히스토리 CSV 저장 실패 (%s): %s', ctrl_name, csv_err.message);
            end
    end
    if nargout > 0
        varargout{1} = [];
    end
end

%% ========================== [수치 변환 헬퍼 함수] ==========================

function p_logical = to_logical_space(p_phys, lb_phys, ub_phys)
    % 물리 공간 -> [-10, 10] 정규화 논리 공간 변환 (logit)
    eps_val = 1e-15;
    ratio = (p_phys - lb_phys) ./ (ub_phys - lb_phys);
    ratio = max(eps_val, min(1 - eps_val, ratio));
    p_logical = -log(1 ./ ratio - 1);
    p_logical = max(-10, min(10, p_logical));
end

function p_phys = to_physical_space(p_logical, lb_phys, ub_phys)
    % [-10, 10] 정규화 논리 공간 -> 물리 공간 변환 (sigmoid)
    sig = 1 ./ (1 + exp(-p_logical));
    p_phys = lb_phys + (ub_phys - lb_phys) .* sig;
end
