% ======================================================================
% [FILE METADATA & VERSION TRACKING]
% - Current Version: v1.0.0 (2026-06-25)
% - Target Environment: MATLAB R2022a or newer
% - Integrity Check: All 4 controllers are digital and optimized via RK4 simulation.
% ======================================================================
% [CHANGELOG - NEVER DELETE THIS HISTORY]
% * v1.0.0 (2026-06-25) - Developer: Gemini AI
%   - Added: Initial version for Buck Converter parameter optimization.
% ======================================================================

function [pi_gains, type3_coeffs, ml_coeffs, lqr_gains] = optimize_controllers(...
    L_val, C_val, G_L, R_C, R_nom, Vin_nom, Vref_val, T_s, t_vec, Vin_data, R_data, Vref_data)

    fprintf('\n>>> [Optimization] 제어기 파라미터 최적화 루틴 시작 <<<\n');
    opt_options = optimset('Display', 'off', 'MaxIter', 150, 'MaxFunEvals', 250);

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

    %% 1. PI 제어기 최적화
    fprintf('- PI 제어기 최적화 중...\n');
    pi_init = [0.08, 500]; % [Kp, Ki]
    pi_cost_fn = @(p) evaluate_pi_cost(p, sys);
    try
        pi_opt = fminsearch(pi_cost_fn, pi_init, opt_options);
        pi_opt(1) = max(0.001, min(2.0, pi_opt(1)));
        pi_opt(2) = max(0.1, min(10000, pi_opt(2)));
    catch
        pi_opt = pi_init;
    end
    pi_gains.KP = pi_opt(1);
    pi_gains.KI = pi_opt(2);
    fprintf('  => 최적 PI 파라미터: KP = %.4f, KI = %.4f\n', pi_gains.KP, pi_gains.KI);

    %% 2. Type 3 (3P2Z) 제어기 최적화 (k-factor 접근법)
    fprintf('- Type 3 (3P2Z) 제어기 최적화 중...\n');
    % 최적화 변수: [f_co, PM_target]
    type3_init = [5000, 60]; 
    type3_cost_fn = @(p) evaluate_type3_cost(p, sys);
    try
        type3_opt = fminsearch(type3_cost_fn, type3_init, opt_options);
        type3_opt(1) = max(500, min(25000, type3_opt(1)));
        type3_opt(2) = max(30, min(85, type3_opt(2)));
    catch
        type3_opt = type3_init;
    end
    % 최적 설계값 적용하여 계수 산출
    [n1, n2, n3, d1, d2, d3, d4] = design_type3_kfactor(type3_opt(1), type3_opt(2), sys);
    type3_coeffs.n1 = n1; type3_coeffs.n2 = n2; type3_coeffs.n3 = n3;
    type3_coeffs.d1 = d1; type3_coeffs.d2 = d2; type3_coeffs.d3 = d3; type3_coeffs.d4 = d4;
    fprintf('  => 최적 Type 3 설계: f_co = %.1f Hz, PM_target = %.1f deg\n', type3_opt(1), type3_opt(2));

    %% 3. ML 최적화 제어기 (5차 전달함수) 최적화
    fprintf('- ML 5차 전달함수 제어기 최적화 중...\n');
    % 최적화 변수: [log_Kc, log_wz1, zeta_z1, log_wz2, zeta_z2, log_wz3, log_wp1, zeta_p1, log_wp2, zeta_p2]
    % 초기값 세팅 (Crossover 주파수 및 댐핑 계수 기반 설계점)
    ml_init = [11.5, 9.4, 0.7, 9.8, 0.7, 8.0, 11.4, 0.7, 12.0, 0.7];
    ml_cost_fn = @(p) evaluate_ml_cost(p, sys);
    try
        ml_opt = fminsearch(ml_cost_fn, ml_init, opt_options);
    catch
        ml_opt = ml_init;
    end
    % 최적 설계값 적용하여 계수 산출
    [m1, m2, m3, m4, m5, m6, e1, e2, e3, e4, e5, e6] = design_ml_tf(ml_opt, sys.T_s);
    ml_coeffs.m1 = m1; ml_coeffs.m2 = m2; ml_coeffs.m3 = m3; ml_coeffs.m4 = m4; ml_coeffs.m5 = m5; ml_coeffs.m6 = m6;
    ml_coeffs.e1 = e1; ml_coeffs.e2 = e2; ml_coeffs.e3 = e3; ml_coeffs.e4 = e4; ml_coeffs.e5 = e5; ml_coeffs.e6 = e6;
    fprintf('  => 최적 ML TF 파라미터 도출 완료 (고차 전달함수 적용)\n');

    %% 4. 현대 제어기 (Augmented LQR) 최적화
    fprintf('- Augmented LQR 제어기 최적화 중...\n');
    % 최적화 변수: [log_q1, log_q2, log_q3, log_R]
    lqr_init = [log(10), log(100), log(1e7), log(1)];
    lqr_cost_fn = @(p) evaluate_lqr_cost(p, sys);
    try
        lqr_opt = fminsearch(lqr_cost_fn, lqr_init, opt_options);
    catch
        lqr_opt = lqr_init;
    end
    % 최적 설계값 적용하여 피드백 게인 산출
    [K_lqr1, K_lqr2, K_lqr3] = design_lqr(lqr_opt, sys);
    lqr_gains.K_lqr1 = K_lqr1;
    lqr_gains.K_lqr2 = K_lqr2;
    lqr_gains.K_lqr3 = K_lqr3;
    fprintf('  => 최적 LQR 파라미터: K_lqr1 = %.4f, K_lqr2 = %.4f, K_lqr3 = %.4f\n', K_lqr1, K_lqr2, K_lqr3);
    fprintf('>>> 제어기 파라미터 최적화 완료 <<<\n\n');
end

%% ========================== [보조 함수] 최적화 비용 함수 및 설계 루틴 ==========================

% 1. PI 비용 함수
function cost = evaluate_pi_cost(p, sys)
    Kp = p(1); Ki = p(2);
    if Kp < 0 || Ki < 0
        cost = 1e12;
        return;
    end
    
    dt = sys.T_s;
    N = length(sys.t_vec);
    V_out_hist = zeros(N, 1);
    
    x = [0; 0];
    error_int = 0;
    duty = sys.Vref_val / sys.Vin_nom;
    
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
        
        err = sys.Vref_data(k) - V_out;
        error_int = error_int + err * dt;
        duty_next = Kp * err + Ki * error_int;
        duty_next = max(0.01, min(0.95, duty_next));
        
        % RK4 Step
        x = rk4_step(A_k, B_k, x, v_sw, dt);
        duty = duty_next;
    end
    
    if any(isnan(V_out_hist)) || any(isinf(V_out_hist))
        cost = 1e12;
        return;
    end
    
    % ITAE 비용 계산
    cost = sum(sys.t_vec .* abs(sys.Vref_data - V_out_hist) * dt);
    
    % 과도 상태 이후의 overshoot에 대한 과도 패널티 추가
    overshoot = max(V_out_hist) - sys.Vref_val;
    if overshoot > 0.05 * sys.Vref_val
        cost = cost + 100 * (overshoot - 0.05 * sys.Vref_val);
    end
end

% 2. Type 3 비용 함수
function cost = evaluate_type3_cost(p, sys)
    f_co = p(1); PM_target = p(2);
    if f_co < 100 || f_co > 40000 || PM_target < 10 || PM_target > 89
        cost = 1e12;
        return;
    end
    
    [n1, n2, n3, d1, d2, d3, d4] = design_type3_kfactor(f_co, PM_target, sys);
    if any(isnan([n1 n2 n3 d1 d2 d3 d4])) || any(isinf([n1 n2 n3 d1 d2 d3 d4]))
        cost = 1e12;
        return;
    end
    
    dt = sys.T_s;
    N = length(sys.t_vec);
    V_out_hist = zeros(N, 1);
    
    x = [0; 0];
    duty = sys.Vref_val / sys.Vin_nom;
    u_hist = [duty, duty, duty];
    err_hist = [0, 0, 0];
    
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
        
        err = sys.Vref_data(k) - V_out;
        err_hist = [err, err_hist(1:2)];
        
        duty_next = n1*err_hist(1) + n2*err_hist(2) + n3*err_hist(3) ...
                  - d2*u_hist(1) - d3*u_hist(2) - d4*u_hist(3);
        duty_next = max(0.01, min(0.95, duty_next));
        u_hist = [duty_next, u_hist(1:2)];
        
        x = rk4_step(A_k, B_k, x, v_sw, dt);
        duty = duty_next;
    end
    
    if any(isnan(V_out_hist)) || any(isinf(V_out_hist))
        cost = 1e12;
        return;
    end
    
    cost = sum(sys.t_vec .* abs(sys.Vref_data - V_out_hist) * dt);
    overshoot = max(V_out_hist) - sys.Vref_val;
    if overshoot > 0.05 * sys.Vref_val
        cost = cost + 100 * (overshoot - 0.05 * sys.Vref_val);
    end
end

% 3. ML 비용 함수
function cost = evaluate_ml_cost(p, sys)
    % 제약 조건 검사: 댐핑 계수 음수 방지 및 물리적 범위 제한
    zeta_z1 = p(3); zeta_z2 = p(5); zeta_p1 = p(8); zeta_p2 = p(10);
    if zeta_z1 < 0.05 || zeta_z1 > 2.0 || zeta_z2 < 0.05 || zeta_z2 > 2.0 || ...
       zeta_p1 < 0.05 || zeta_p1 > 2.0 || zeta_p2 < 0.05 || zeta_p2 > 2.0
        cost = 1e12;
        return;
    end
    
    [m1, m2, m3, m4, m5, m6, e1, e2, e3, e4, e5, e6] = design_ml_tf(p, sys.T_s);
    if any(isnan([m1 m2 m3 m4 m5 m6 e1 e2 e3 e4 e5 e6])) || any(isinf([m1 m2 m3 m4 m5 m6 e1 e2 e3 e4 e5 e6]))
        cost = 1e12;
        return;
    end
    
    dt = sys.T_s;
    N = length(sys.t_vec);
    V_out_hist = zeros(N, 1);
    
    x = [0; 0];
    duty = sys.Vref_val / sys.Vin_nom;
    u_hist = ones(5, 1) * duty;
    err_hist = zeros(6, 1);
    
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
        
        err = sys.Vref_data(k) - V_out;
        err_hist = [err; err_hist(1:5)];
        
        duty_next = (m1*err_hist(1) + m2*err_hist(2) + m3*err_hist(3) + m4*err_hist(4) + m5*err_hist(5) + m6*err_hist(6) ...
                   - e2*u_hist(1) - e3*u_hist(2) - e4*u_hist(3) - e5*u_hist(4) - e6*u_hist(5)) / e1;
        duty_next = max(0.01, min(0.95, duty_next));
        u_hist = [duty_next; u_hist(1:4)];
        
        x = rk4_step(A_k, B_k, x, v_sw, dt);
        duty = duty_next;
    end
    
    if any(isnan(V_out_hist)) || any(isinf(V_out_hist))
        cost = 1e12;
        return;
    end
    
    cost = sum(sys.t_vec .* abs(sys.Vref_data - V_out_hist) * dt);
    overshoot = max(V_out_hist) - sys.Vref_val;
    if overshoot > 0.05 * sys.Vref_val
        cost = cost + 100 * (overshoot - 0.05 * sys.Vref_val);
    end
end

% 4. LQR 비용 함수
function cost = evaluate_lqr_cost(p, sys)
    q1 = exp(p(1)); q2 = exp(p(2)); q3 = exp(p(3)); R_weight = exp(p(4));
    
    % LQR 제어 법칙 설계
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
    
    try
        K_lqr_all = lqr(A_aug, B_aug, Q_lqr, R_weight);
        K_lqr1 = K_lqr_all(1);
        K_lqr2 = K_lqr_all(2);
        K_lqr3 = K_lqr_all(3);
    catch
        cost = 1e12;
        return;
    end
    
    dt = sys.T_s;
    N = length(sys.t_vec);
    V_out_hist = zeros(N, 1);
    
    x = [0; 0];
    error_int = 0;
    duty = sys.Vref_val / sys.Vin_nom;
    
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
        
        err = sys.Vref_data(k) - V_out;
        error_int = error_int + err * dt;
        
        I_L_ref = sys.Vref_data(k) / R_k;
        duty_nom = sys.Vref_data(k) / V_in_k;
        
        duty_next = duty_nom - ( K_lqr1 * (x(1) - I_L_ref) + K_lqr2 * (x(2) - sys.Vref_data(k)) + K_lqr3 * error_int );
        duty_next = max(0.01, min(0.95, duty_next));
        
        x = rk4_step(A_k, B_k, x, v_sw, dt);
        duty = duty_next;
    end
    
    if any(isnan(V_out_hist)) || any(isinf(V_out_hist))
        cost = 1e12;
        return;
    end
    
    cost = sum(sys.t_vec .* abs(sys.Vref_data - V_out_hist) * dt);
    overshoot = max(V_out_hist) - sys.Vref_val;
    if overshoot > 0.05 * sys.Vref_val
        cost = cost + 100 * (overshoot - 0.05 * sys.Vref_val);
    end
end

%% ========================== [설계 알고리즘] ==========================

% Type 3 k-factor 설계
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
    
    % 분모/분자 각 컴포넌트별 Backward Euler 이산화 적용
    % 1) G1(s) = Kc * (s + wz3)/s => G1(z) = Kc * ((1 + wz3*Ts) - z^-1) / (1 - z^-1)
    num1 = [Kc * (1 + wz3 * T_s), -Kc];
    den1 = [1, -1];
    
    % 2) G2(s) = (s^2 + a1*s + a0) / (s^2 + b1*s + b0)
    a1 = 2 * zeta_z1 * wz1; a0 = wz1^2;
    b1 = 2 * zeta_p1 * wp1; b0 = wp1^2;
    num2 = [1 + a1*T_s + a0*T_s^2, -(2 + a1*T_s), 1];
    den2 = [1 + b1*T_s + b0*T_s^2, -(2 + b1*T_s), 1];
    
    % 3) G3(s) = (s^2 + a'1*s + a'0) / (s^2 + b'1*s + b'0)
    a_prime1 = 2 * zeta_z2 * wz2; a_prime0 = wz2^2;
    b_prime1 = 2 * zeta_p2 * wp2; b_prime0 = wp2^2;
    num3 = [1 + a_prime1*T_s + a_prime0*T_s^2, -(2 + a_prime1*T_s), 1];
    den3 = [1 + b_prime1*T_s + b_prime0*T_s^2, -(2 + b_prime1*T_s), 1];
    
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
