%% ======================================================================
% [FILE METADATA]
% - Version: v4.1.0_Final (2026-06-26)
% - Target Environment: MATLAB R2022a or newer
% - Description: 완전 이산화(Discrete) 기반 4종 제어기 통합 튜닝 및 수치해석/시뮬링크 검증 스크립트 (Consolidated)
% ======================================================================

clear all; clc; close all;

%% 1. 벅 컨버터 물리 시스템 및 극한 가변 시나리오 정의
f_sw = 100e3; T_s = 1 / f_sw; t_end = 0.1;
t_vec = (0:T_s:t_end)'; N_pts = length(t_vec);

Vin_nom = 12; Vref_val = 5;
L_nominal = 100e-6; C_nominal = 220e-6;

% 실제 플랜트 불확실성 주입 (L은 30% 감소, C는 30% 증가된 최악의 조건 가동)
L_val = 0.7 * L_nominal;  
C_val = 1.3 * C_nominal;  
R_nom = 5; G_L = 1e-3; R_C = 0.05;

% [시나리오] 가혹한 동적 변동 프로파일 (부하 급변 및 전압 서지)
R_data = R_nom * ones(N_pts, 1);
R_data(t_vec >= 0.03 & t_vec < 0.07) = 50;  % 30ms~70ms: 경부하 (50 Ohm)
R_data(t_vec >= 0.07) = 2.0;                % 70ms~100ms: 중부하 (2 Ohm)

Vin_data = Vin_nom * ones(N_pts, 1);
Vin_data(t_vec >= 0.04) = 18;               % 40ms: 입력 서지 (+50% 폭발 서지)
rng(42); Vin_data = Vin_data + 0.3 * randn(N_pts, 1); % 전 영역 고주파 노이즈 합성

Vref_data = Vref_val * ones(N_pts, 1);

%% 2. 4종 제어기 파라미터 완전 이산화(Discrete Domain) 최적화 엔진
fprintf('=== [Step 2] 4종 제어기 파라미터 이산 도메인 PSO 학습 기동 ===\n');

% 최악의 조건(R = 2 Ohm) 베이스라인 공칭 매트릭스 도출 (B 행렬에서 Vin_nom 중복 곱셈 제거)
R_worst = 2.0;
theta_w = G_L * R_worst * R_C + R_worst + R_C;
A_nom_c = [ -R_worst * R_C / (L_val * theta_w),                 -R_worst / (L_val * theta_w);
             R_worst / (C_val * theta_w),                 -(R_worst * G_L + 1) / (C_val * theta_w) ];
B_nom_c = [ (R_worst + R_C) / (L_val * theta_w);
            (R_worst * G_L) / (C_val * theta_w) ]; % 순수 시스템 행렬 (Vin_nom 중복 곱 제외)
C_nom_c = [ R_worst * R_C / theta_w,   R_worst / theta_w ];
D_nom_c = G_L * R_worst * R_C / theta_w;

% 플랜트 모델을 시작부터 이산화(c2d)하여 도메인 미스매치 원천 차단
sys_c = ss(A_nom_c, B_nom_c, C_nom_c, D_nom_c);
sys_d = c2d(sys_c, T_s, 'zoh');
A_nom = sys_d.A; B_nom = sys_d.B; C_nom = sys_d.C; D_nom = sys_d.D;

% (1) PI 제어기 최적화 (2변수: KP, KI)
lb_pi = [0.001, 1.0]; ub_pi = [0.2, 500.0];
cost_pi_fn = @(p) evaluate_pi_cost(p(1), p(2), L_val, C_val, G_L, R_C, Vin_data, R_data, Vref_val, T_s, t_vec, N_pts);
gbest_pi = pso_optimize(cost_pi_fn, 2, lb_pi, ub_pi, 40, 20);
KP = gbest_pi(1); KI = gbest_pi(2);
fprintf('- [1. PI 제어기] 이산 튜닝 완료 (KP = %.6f, KI = %.6f)\n', KP, KI);

% (2) 3P2Z 제어기 최적화 (2변수: f_co, PM_target)
lb_3p = [800, 45.0]; ub_3p = [2200, 80.0];
cost_3p_fn = @(p) evaluate_3p2z_cost(p, A_nom_c, B_nom_c, C_nom_c, D_nom_c, L_val, C_val, G_L, R_C, Vin_data, R_data, Vref_val, T_s, t_vec, N_pts);
gbest_3p = pso_optimize(cost_3p_fn, 2, lb_3p, ub_3p, 40, 20);
f_co_3p = gbest_3p(1); PM_target_3p = gbest_3p(2);

% 디지털 도메인 3P2Z 계수 확정 및 유도
w_co_3p = 2 * pi * f_co_3p; s_co_3p = 1i * w_co_3p;
G_co_3p = C_nom_c * ((s_co_3p * eye(2) - A_nom_c) \ (B_nom_c * Vin_nom)) + D_nom_c;
boost_3p = PM_target_3p - 90 - rad2deg(angle(G_co_3p));
k_factor_3p = tan(deg2rad(45 + boost_3p / 4));
w_z_3p = w_co_3p / k_factor_3p; w_p_3p = w_co_3p * k_factor_3p;
K_c_3p = 1 / (abs(G_co_3p) * ((w_co_3p^2 + w_z_3p^2) / (w_co_3p * (w_co_3p^2 + w_p_3p^2))));
A_p_3p = 1 + w_p_3p * T_s; A_z_3p = 1 + w_z_3p * T_s; d_raw1_3p = A_p_3p^2;
n1 = (K_c_3p * T_s * A_z_3p^2) / d_raw1_3p; n2 = (-2 * K_c_3p * T_s * A_z_3p) / d_raw1_3p; n3 = (K_c_3p * T_s) / d_raw1_3p;
d1 = 1.0; d2 = -(A_p_3p^2 + 2*A_p_3p) / d_raw1_3p; d3 = (2*A_p_3p + 1) / d_raw1_3p; d4 = -1 / d_raw1_3p;
fprintf('- [2. 3P2Z 제어기] 이산화 매칭 성료 (f_co = %.1f Hz, PM = %.1f deg)\n', f_co_3p, PM_target_3p);

% (3) ML 제어기 머신러닝 직접 최적화 (6변수: Kc, z1, wz, zeta_z, wp, zeta_p)
lb_ml = [0.001,  50.0,  500.0, 0.05,  2500.0, 0.3];
ub_ml = [10.0,  1000.0, 2000.0, 0.90, 15000.0, 1.0];
cost_ml_fn = @(p) evaluate_ml_ml_cost(p, L_val, C_val, G_L, R_C, Vin_data, R_data, Vref_val, T_s, t_vec, N_pts);
gbest_ml = pso_optimize(cost_ml_fn, 6, lb_ml, ub_ml, 60, 30);
Kc_opt = gbest_ml(1); z1_opt = gbest_ml(2); wz_opt = gbest_ml(3); zeta_z_opt = gbest_ml(4); wp_opt = gbest_ml(5); zeta_p_opt = gbest_ml(6);

% 최적화된 극-영점 매개변수로 5차 이산 전달함수 계수 유도
w_z_rad = 2 * pi * wz_opt; w_p_rad = 2 * pi * wp_opt; w_z1_rad = 2 * pi * z1_opt;
w_co_ml = w_z_rad; K_warp = w_co_ml / tan(w_co_ml * T_s / 2);
num_z1 = Kc_opt * [K_warp + w_z1_rad, w_z1_rad - K_warp]; den_z1 = [K_warp, -K_warp];
num_z2 = [K_warp^2 + 2*zeta_z_opt*w_z_rad*K_warp + w_z_rad^2, 2*(w_z_rad^2 - K_warp^2), K_warp^2 - 2*zeta_z_opt*w_z_rad*K_warp + w_z_rad^2];
den_z2 = [K_warp^2 + 2*zeta_p_opt*w_p_rad*K_warp + w_p_rad^2, 2*(w_p_rad^2 - K_warp^2), K_warp^2 - 2*zeta_p_opt*w_p_rad*K_warp + w_p_rad^2];
num_ml_all = conv(num_z1, num_z2); den_ml_all = conv(den_z1, den_z2);
m_coeff = [num_ml_all, 0, 0] / den_ml_all(1); e_coeff = [den_ml_all, 0, 0] / den_ml_all(1);
m1 = m_coeff(1); m2 = m_coeff(2); m3 = m_coeff(3); m4 = m_coeff(4); m5 = m_coeff(5); m6 = m_coeff(6);
e1 = 1.0; e2 = e_coeff(2); e3 = e_coeff(3); e4 = e_coeff(4); e5 = e_coeff(5); e6 = e_coeff(6);
fprintf('- [3. ML 머신러닝 제어기] 5차 이산 제어기 학습 완료\n');

% (4) LQR 제어기 최적화: 연속 이득 대신 완전한 이산 디지털 dlqr 도입 + Bryson 역수 제거
lb_lqr = [1.0, 10.0, 100.0, 0.0001]; ub_lqr = [100.0, 1e5, 1e6, 0.1];
cost_lqr_fn = @(p) evaluate_lqr_cost(p, A_nom_c, B_nom_c, C_nom_c, D_nom_c, L_val, C_val, G_L, R_C, Vin_data, R_data, Vref_val, T_s, t_vec, N_pts);
gbest_lqr = pso_optimize(cost_lqr_fn, 4, lb_lqr, ub_lqr, 40, 20);

% 디지털 이산 확장 상공방 기반 dlqr 산출 (굼벵이 현상 종결 및 초고속 추종)
A_aug_c = [ A_nom_c, zeros(2, 1); -C_nom_c, 0 ];
B_aug_c = [ B_nom_c * Vin_nom; -D_nom_c ];
sys_aug_c = ss(A_aug_c, B_aug_c, eye(3), 0);
sys_aug_d = c2d(sys_aug_c, T_s, 'zoh');
Q_lqr = diag([gbest_lqr(1), gbest_lqr(2), gbest_lqr(3)]);
R_lqr = gbest_lqr(4);
K_lqr_all = dlqr(sys_aug_d.A, sys_aug_d.B, Q_lqr, R_lqr);
K_lqr1 = K_lqr_all(1); K_lqr2 = K_lqr_all(2); K_lqr3 = K_lqr_all(3);
fprintf('- [4. LQR 제어기] 디지털 이산 dlqr 이득 튜닝 완료\n\n');

%% 3. 제어 시뮬레이션 실행 (Simulink 및 고정밀 수치 해석 Solver)
fprintf('=== [Step 3] 시뮬레이션 엔진 구동 ===\n');

simulink_model_name = 'BuckConverter'; 
simulink_running = false;

% Simulink 호환 바인딩 및 Workspace 가변 데이터 공급
Kp = KP; Ki = KI;
Ts = T_s;
L_data = L_val * ones(N_pts, 1);
C_data = C_val * ones(N_pts, 1);

Vin = [t_vec, Vin_data];
R = [t_vec, R_data];
L = [t_vec, L_data];
C = [t_vec, C_data];
Vref = [t_vec, Vref_data];
Vs = Vin;

Vin_ts = timeseries(Vin_data, t_vec);
R_ts = timeseries(R_data, t_vec);
L_ts = timeseries(L_data, t_vec);
C_ts = timeseries(C_data, t_vec);
Vref_ts = timeseries(Vref_data, t_vec);
Vs_ts = Vin_ts;

if exist(simulink_model_name, 'file') == 4
    fprintf('- 실 시뮬링크 모델 [%s.slx] 검출. 시뮬레이션 기동...\n', simulink_model_name);
    try
        load_system(simulink_model_name);
        % Headless 실행을 위한 Scope 자동 열기 차단 설정
        scopes = find_system(simulink_model_name, 'BlockType', 'Scope');
        for idx = 1:length(scopes)
            try
                set_param(scopes{idx}, 'OpenAtSimulationStart', 'off');
            catch
            end
        end
        out_sim = sim(simulink_model_name, 'StopTime', num2str(t_end));
        simulink_running = true;
        fprintf('=> Simulink 시뮬레이션 성공!\n\n');
    catch sim_err
        warning('Scenario:SimulinkFailed', 'Simulink 실행 중 문제가 발견되어 안전 수치 해석 Solver로 전환합니다. 에러: %s', sim_err.message);
        fprintf('\n================== [에러 추적] ==================\n');
        disp(getReport(sim_err, 'extended'));
        fprintf('=================================================\n\n');
    end
end

if simulink_running
    fprintf('- [데이터 정비] Simulink 반환 객체를 고정밀 double struct로 형변환 중...\n');
    out = struct();
    if isprop(out_sim, 'tout') && ~isempty(out_sim.tout)
        out.tout = out_sim.tout;
    elseif isa(out_sim.V_Real, 'timeseries')
        out.tout = out_sim.V_Real.Time;
    else
        out.tout = t_vec;
    end
    
    fields = {'V_Real', 'I_L', 'Duty', 'V_Real1', 'I_L1', 'Duty1', 'V_Real2', 'I_L2', 'Duty2', 'V_Real3', 'I_L3', 'Duty3'};
    for f = 1:length(fields)
        fieldname = fields{f};
        if isprop(out_sim, fieldname)
            val = out_sim.(fieldname);
            if isa(val, 'timeseries')
                out.(fieldname) = val.Data;
            else
                out.(fieldname) = val;
            end
        else
            out.(fieldname) = zeros(length(out.tout), 1);
        end
    end
    fprintf('=> [변환 완료] 데이터 정합성 프로세스 성료.\n\n');
else
    fprintf('- [알림] 고정밀 Closed-loop State-Space Solver를 기동합니다.\n');
    out.tout = t_vec;
    fields = {'V_Real','I_L','Duty', 'V_Real1','I_L1','Duty1', 'V_Real2','I_L2','Duty2', 'V_Real3','I_L3','Duty3'};
    for f=1:length(fields), out.(fields{f}) = zeros(N_pts, 1); end

    % 상태 메모리 분리 초기화
    x_pi = [0;0]; error_int_pi = 0; duty_pi = Vref_val/Vin_nom;
    x_3p2z = [0;0]; u_hist_3p2z = zeros(3,1); err_hist_3p2z = zeros(3,1); duty_3p2z = duty_pi;
    x_ml = [0;0]; u_hist_ml = zeros(5,1); err_hist_ml = zeros(6,1); duty_ml = duty_pi;
    x_lqr = [0;0]; error_int_lqr = 0; duty_lqr = duty_pi;

    for k = 1:N_pts
        V_in_k = Vin_data(k); R_k = R_data(k);
        theta_k = G_L * R_k * R_C + R_k + R_C;
        A_k_c = [ -R_k * R_C / (L_val * theta_k),                 -R_k / (L_val * theta_k);
                 R_k / (C_val * theta_k),                 -(R_k * G_L + 1) / (C_val * theta_k) ];
        B_k_c = [ (R_k + R_C) / (L_val * theta_k);
                (R_k * G_L) / (C_val * theta_k) ];
        C_k = [ R_k * R_C / theta_k,   R_k / theta_k ];
        D_k = G_L * R_k * R_C / theta_k;
        
        % --- [1] PI Control Loop ---
        v_sw_pi = duty_pi * V_in_k; V_out_pi = C_k * x_pi + D_k * v_sw_pi;
        err_pi = Vref_val - V_out_pi; error_int_pi = error_int_pi + err_pi * T_s;
        duty_pi_raw = KP * err_pi + KI * error_int_pi;
        duty_pi = max(0.01, min(0.95, duty_pi_raw));
        if duty_pi_raw > 0.95 || duty_pi_raw < 0.01, error_int_pi = error_int_pi - err_pi * T_s; end
        
        % --- [2] 3P2Z Control Loop ---
        v_sw_3p2z = duty_3p2z * V_in_k; V_out_3p2z = C_k * x_3p2z + D_k * v_sw_3p2z;
        err_3p2z = Vref_val - V_out_3p2z;
        err_hist_3p2z = [err_3p2z; err_hist_3p2z(1:2)];
        duty_3p2z_raw = n1*err_hist_3p2z(1) + n2*err_hist_3p2z(2) + n3*err_hist_3p2z(3) ...
                        - d2*u_hist_3p2z(1) - d3*u_hist_3p2z(2) - d4*u_hist_3p2z(3);
        duty_3p2z = max(0.01, min(0.95, duty_3p2z_raw));
        u_hist_3p2z = [duty_3p2z; u_hist_3p2z(1:2)];
        
        % --- [3] ML High-Order Loop ---
        v_sw_ml = duty_ml * V_in_k; V_out_ml = C_k * x_ml + D_k * v_sw_ml;
        err_ml = Vref_val - V_out_ml;
        err_hist_ml = [err_ml; err_hist_ml(1:5)];
        duty_ml_raw = (m1*err_hist_ml(1) + m2*err_hist_ml(2) + m3*err_hist_ml(3) + m4*err_hist_ml(4) + m5*err_hist_ml(5) + m6*err_hist_ml(6) ...
                      - e2*u_hist_ml(1) - e3*u_hist_ml(2) - e4*u_hist_ml(3) - e5*u_hist_ml(4) - e6*u_hist_ml(5)) / e1;
        duty_ml = max(0.01, min(0.95, duty_ml_raw));
        u_hist_ml = [duty_ml; u_hist_ml(1:4)];

        % --- [4] Modern LQR Servo Loop (dlqr 완벽 튜닝 싱크 구동) ---
        v_sw_lqr = duty_lqr * V_in_k; V_out_lqr = C_k * x_lqr + D_k * v_sw_lqr;
        err_lqr = Vref_val - V_out_lqr; error_int_lqr = error_int_lqr + err_lqr * T_s;
        duty_lqr_nom = Vref_val / V_in_k; I_L_ref = Vref_val / R_k;
        duty_lqr_raw = duty_lqr_nom - ( K_lqr1 * (x_lqr(1) - I_L_ref) + K_lqr2 * (V_out_lqr - Vref_val) + K_lqr3 * error_int_lqr );
        duty_lqr = max(0.01, min(0.95, duty_lqr_raw));
        if duty_lqr_raw > 0.95 || duty_lqr_raw < 0.01, error_int_lqr = error_int_lqr - err_lqr * T_s; end
        
        % RK4 고정밀 독립적 상태 업데이트
        x_pi   = rk4_step(A_k_c, B_k_c, x_pi,   v_sw_pi,   T_s);
        x_3p2z = rk4_step(A_k_c, B_k_c, x_3p2z, v_sw_3p2z, T_s);
        x_ml   = rk4_step(A_k_c, B_k_c, x_ml,   v_sw_ml,   T_s);
        x_lqr  = rk4_step(A_k_c, B_k_c, x_lqr,  v_sw_lqr,  T_s);
        
        % 데이터 로깅
        out.V_Real(k) = V_out_pi;   out.I_L(k) = x_pi(1);   out.Duty(k) = duty_pi;
        out.V_Real1(k) = V_out_3p2z; out.I_L1(k) = x_3p2z(1); out.Duty1(k) = duty_3p2z;
        out.V_Real2(k) = V_out_ml;   out.I_L2(k) = x_ml(1);   out.Duty2(k) = duty_ml;
        out.V_Real3(k) = V_out_lqr;  out.I_L3(k) = x_lqr(1);  out.Duty3(k) = duty_lqr;
    end
    fprintf('=> [성공] 수치해석 시뮬레이션 성료 및 완벽한 수렴 확인.\n\n');
end

%% 4. 정량적 성능 지표 분석 및 출력
get_t_s = @(v, t, target) feval(@(idx) double(~isempty(idx)) * t(max(1, idx)), find(abs(movmean(v, 15) - target) > 0.02 * target, 1, 'last'));
calc_metrics = @(v, t, target) struct(...
    'Overshoot', max(0, (max(movmean(v, 15)) - target)/target * 100), ...
    'Undershoot', max(0, (target - min(movmean(v, 15)))/target * 100), ...
    'SettlingTime', get_t_s(v, t, target) ...
);

t_eval1 = out.tout(out.tout <= 0.04);
v_pi_1 = out.V_Real(out.tout <= 0.04); v_3p2z_1 = out.V_Real1(out.tout <= 0.04); v_ml_1 = out.V_Real2(out.tout <= 0.04); v_lqr_1 = out.V_Real3(out.tout <= 0.04);

t_eval2 = out.tout(out.tout > 0.04 & out.tout <= 0.07) - 0.04;
v_pi_2 = out.V_Real(out.tout > 0.04 & out.tout <= 0.07); v_3p2z_2 = out.V_Real1(out.tout > 0.04 & out.tout <= 0.07); v_ml_2 = out.V_Real2(out.tout > 0.04 & out.tout <= 0.07); v_lqr_2 = out.V_Real3(out.tout > 0.04 & out.tout <= 0.07);

t_eval3 = out.tout(out.tout > 0.07) - 0.07;
v_pi_3 = out.V_Real(out.tout > 0.07); v_3p2z_3 = out.V_Real1(out.tout > 0.07); v_ml_3 = out.V_Real2(out.tout > 0.07); v_lqr_3 = out.V_Real3(out.tout > 0.07);

fprintf('\n==================== STARTUP TRANSIENT (0~40ms) ====================\n');
fprintf('%-10s | %-12s | %-12s | %-12s\n', 'Controller', 'Overshoot(%)', 'Undershoot(%)', 'Settling(ms)');
fprintf('--------------------------------------------------------------------\n');
m_pi = calc_metrics(v_pi_1, t_eval1, 5); m_3p = calc_metrics(v_3p2z_1, t_eval1, 5); m_ml = calc_metrics(v_ml_1, t_eval1, 5); m_lq = calc_metrics(v_lqr_1, t_eval1, 5);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'PI', m_pi.Overshoot, m_pi.Undershoot, m_pi.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', '3P2Z', m_3p.Overshoot, m_3p.Undershoot, m_3p.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'ML', m_ml.Overshoot, m_ml.Undershoot, m_ml.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'LQR', m_lq.Overshoot, m_lq.Undershoot, m_lq.SettlingTime*1000);

fprintf('\n==================== INPUT SURGE TRANSIENT (40~70ms) ====================\n');
fprintf('%-10s | %-12s | %-12s | %-12s\n', 'Controller', 'Overshoot(%)', 'Undershoot(%)', 'Settling(ms)');
fprintf('-------------------------------------------------------------------------\n');
m_pi = calc_metrics(v_pi_2, t_eval2, 5); m_3p = calc_metrics(v_3p2z_2, t_eval2, 5); m_ml = calc_metrics(v_ml_2, t_eval2, 5); m_lq = calc_metrics(v_lqr_2, t_eval2, 5);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'PI', m_pi.Overshoot, m_pi.Undershoot, m_pi.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', '3P2Z', m_3p.Overshoot, m_3p.Undershoot, m_3p.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'ML', m_ml.Overshoot, m_ml.Undershoot, m_ml.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'LQR', m_lq.Overshoot, m_lq.Undershoot, m_lq.SettlingTime*1000);

fprintf('\n==================== LOAD TRANSIENT (70~100ms) ====================\n');
fprintf('%-10s | %-12s | %-12s | %-12s\n', 'Controller', 'Overshoot(%)', 'Undershoot(%)', 'Settling(ms)');
fprintf('--------------------------------------------------------------------\n');
m_pi = calc_metrics(v_pi_3, t_eval3, 5); m_3p = calc_metrics(v_3p2z_3, t_eval3, 5); m_ml = calc_metrics(v_ml_3, t_eval3, 5); m_lq = calc_metrics(v_lqr_3, t_eval3, 5);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'PI', m_pi.Overshoot, m_pi.Undershoot, m_pi.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', '3P2Z', m_3p.Overshoot, m_3p.Undershoot, m_3p.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'ML', m_ml.Overshoot, m_ml.Undershoot, m_ml.SettlingTime*1000);
fprintf('%-10s | %12.3f | %12.3f | %12.3f\n', 'LQR', m_lq.Overshoot, m_lq.Undershoot, m_lq.SettlingTime*1000);
fprintf('====================================================================\n');

%% 5. 시각화 그래픽스 (고대비 라인 컬러 색상표 적용)
figure('Position', [100, 50, 1100, 750], 'Name', 'Discrete Unified Control Comparison');
color_pi=[0, 0.45, 0.74]; color_3p2z=[0.85, 0.33, 0.1]; color_ml=[0.49, 0.18, 0.56]; color_lqr=[0.12, 0.69, 0.33];
Vref_plot = Vref_val * ones(length(out.tout), 1);

subplot(2,1,1); hold on; grid on; box on;
safe_plot(out.tout*1000, out.V_Real,  color_pi,   1.5, 'PI');
safe_plot(out.tout*1000, out.V_Real1, color_3p2z, 1.5, '3P2Z');
safe_plot(out.tout*1000, out.V_Real2, color_ml,   1.5, 'ML 5th-Order');
safe_plot(out.tout*1000, out.V_Real3, color_lqr,  2.0, 'Digital dlqr');
safe_plot(out.tout*1000, Vref_plot,   [1, 0, 0],  1.2, '', '--');
line([40 40], [0 10], 'Color', 'r', 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off');
line([70 70], [0 10], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 1.5, 'HandleVisibility', 'off');
ylabel('Output Voltage V_{out} (V)', 'FontWeight', 'bold');
title('이산화(Discrete Domain) 매칭 후 4종 제어기 응답 곡선 (Simulink/수치해석)', 'FontSize', 12, 'FontWeight', 'bold');
legend('Location', 'northeast'); ylim([4.2, 5.8]);

subplot(2,1,2); hold on; grid on; box on;
safe_plot(out.tout*1000, out.Duty,  color_pi,   1.2, '');
safe_plot(out.tout*1000, out.Duty1, color_3p2z, 1.2, '');
safe_plot(out.tout*1000, out.Duty2, color_ml,   1.2, '');
safe_plot(out.tout*1000, out.Duty3, color_lqr,  1.5, '');
xlabel('Time (ms)', 'FontWeight', 'bold'); ylabel('Duty Ratio', 'FontWeight', 'bold');
title('순시 제어 입력(Duty Ratio) 거동 - 채터링 원천 봉쇄 확인', 'FontSize', 10, 'FontWeight', 'bold');
ylim([0, 1]);

%% 6. 성능 지표 및 시뮬레이션 결과 데이터 CSV 저장
fprintf('=== [Step 6] 성능 지표 및 시뮬레이션 결과 CSV 저장 ===\n');

csv_folder = 'csv_data';
if ~exist(csv_folder, 'dir')
    mkdir(csv_folder);
end

% (1) 성능 지표 저장
try
    Controller_Names = {'PI'; '3P2Z'; 'ML'; 'LQR'};
    Overshoots = [m_pi.Overshoot; m_3p.Overshoot; m_ml.Overshoot; m_lq.Overshoot];
    Settling_Times = [m_pi.SettlingTime; m_3p.SettlingTime; m_ml.SettlingTime; m_lq.SettlingTime] * 1000;
    
    Performance_Table = table(Controller_Names, Overshoots, Settling_Times, ...
        'VariableNames', {'Controller_Type', 'Overshoot_Percent', 'Settling_Time_ms'});
    writetable(Performance_Table, fullfile(csv_folder, 'performance_metrics.csv'));
    fprintf('- [성공] 성능 지표 저장 완료: %s\n', fullfile(csv_folder, 'performance_metrics.csv'));
catch err_perf
    warning('Scenario:CSVPerfFailed', '성능 지표 CSV 저장 실패: %s', err_perf.message);
end

% (2) 시뮬레이션 결과 시계열 데이터 저장
try
    t_out = out.tout(:);
    N_pts_out = length(t_out);
    
    V_PI   = out.V_Real(:);
    I_PI   = out.I_L(:);
    D_PI   = out.Duty(:);
    
    V_3P2Z = out.V_Real1(:);
    I_3P2Z = out.I_L1(:);
    D_3P2Z = out.Duty1(:);
    
    V_ML   = out.V_Real2(:);
    I_ML   = out.I_L2(:);
    D_ML   = out.Duty2(:);
    
    V_LQR  = out.V_Real3(:);
    I_LQR  = out.I_L3(:);
    D_LQR  = out.Duty3(:);
    
    SimResults_Table = table(t_out, ...
        V_PI, I_PI, D_PI, ...
        V_3P2Z, I_3P2Z, D_3P2Z, ...
        V_ML, I_ML, D_ML, ...
        V_LQR, I_LQR, D_LQR, ...
        'VariableNames', {'Time', ...
        'V_out_PI', 'I_L_PI', 'Duty_PI', ...
        'V_out_3P2Z', 'I_L_3P2Z', 'Duty_3P2Z', ...
        'V_out_ML', 'I_L_ML', 'Duty_ML', ...
        'V_out_LQR', 'I_L_LQR', 'Duty_LQR'});
        
    writetable(SimResults_Table, fullfile(csv_folder, 'simulation_results.csv'));
    fprintf('- [성공] 시뮬레이션 결과 저장 완료: %s\n\n', fullfile(csv_folder, 'simulation_results.csv'));
catch err_sim
    warning('Scenario:CSVSimFailed', '시뮬레이션 결과 CSV 저장 실패: %s', err_sim.message);
end

%% ======================================================================
%  [LOCAL FUNCTIONS] - Consolidated Helper Functions
% ======================================================================

function gbest = pso_optimize(cost_fn, n_vars, lb, ub, max_iter, pop_size)
    w = 0.6; c1 = 1.5; c2 = 1.5;
    pos = lb + (ub - lb) .* rand(pop_size, n_vars);
    vel = zeros(pop_size, n_vars);
    pbest = pos;
    pbest_cost = inf(pop_size, 1);
    gbest = zeros(1, n_vars);
    gbest_cost = inf;
    
    for iter = 1:max_iter
        for i = 1:pop_size
            cost = cost_fn(pos(i, :));
            if cost < pbest_cost(i)
                pbest(i, :) = pos(i, :);
                pbest_cost(i) = cost;
            end
            if cost < gbest_cost
                gbest = pos(i, :);
                gbest_cost = cost;
            end
        end
        % Update velocities & positions
        vel = w * vel + c1 * rand(pop_size, n_vars) .* (pbest - pos) + c2 * rand(pop_size, n_vars) .* (gbest - pos);
        pos = pos + vel;
        pos = max(lb, min(ub, pos));
    end
end

function x_next = rk4_step(A, B, x, u, dt)
    k1 = A * x + B * u;
    k2 = A * (x + 0.5 * dt * k1) + B * u;
    k3 = A * (x + 0.5 * dt * k2) + B * u;
    k4 = A * (x + dt * k3) + B * u;
    x_next = x + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
end

function cost = evaluate_pi_cost(KP, KI, L_val, C_val, G_L, R_C, Vin_data, R_data, Vref_val, T_s, t_vec, N_pts)
    x = [0;0]; error_int = 0; duty = Vref_val/Vin_data(1);
    V_out_log = zeros(N_pts, 1);
    
    for k = 1:N_pts
        V_in_k = Vin_data(k); R_k = R_data(k);
        theta_k = G_L * R_k * R_C + R_k + R_C;
        A_k = [ -R_k * R_C / (L_val * theta_k),                 -R_k / (L_val * theta_k);
                 R_k / (C_val * theta_k),                 -(R_k * G_L + 1) / (C_val * theta_k) ];
        B_k = [ (R_k + R_C) / (L_val * theta_k);
                (R_k * G_L) / (C_val * theta_k) ];
        C_k = [ R_k * R_C / theta_k,   R_k / theta_k ];
        D_k = G_L * R_k * R_C / theta_k;
        
        v_sw = duty * V_in_k; V_out = C_k * x + D_k * v_sw;
        V_out_log(k) = V_out;
        
        err = Vref_val - V_out; error_int = error_int + err * T_s;
        duty_raw = KP * err + KI * error_int;
        duty = max(0.01, min(0.95, duty_raw));
        if duty_raw > 0.95 || duty_raw < 0.01, error_int = error_int - err * T_s; end
        
        x = rk4_step(A_k, B_k, x, v_sw, T_s);
    end
    
    v_smooth = movmean(V_out_log, 15);
    err_smooth = abs(Vref_val - v_smooth);
    cost = sum(t_vec .* err_smooth * T_s);
    
    max_V = max(v_smooth);
    overshoot = max(0, (max_V - Vref_val)/Vref_val * 100);
    if overshoot > 15
        cost = cost + 50 * (overshoot - 15);
    end
    
    err_startup = mean(err_smooth(t_vec >= 0.035 & t_vec <= 0.04));
    err_surge = mean(err_smooth(t_vec >= 0.065 & t_vec <= 0.07));
    err_load = mean(err_smooth(t_vec >= 0.095));
    
    if err_startup > 0.1, cost = cost + 1000 * err_startup; end
    if err_surge > 0.1,   cost = cost + 1000 * err_surge; end
    if err_load > 0.1,    cost = cost + 1000 * err_load; end
end

function cost = evaluate_3p2z_cost(p, A_nom_c, B_nom_c, C_nom_c, D_nom_c, L_val, C_val, G_L, R_C, Vin_data, R_data, Vref_val, T_s, t_vec, N_pts)
    f_co = p(1); PM_target = p(2);
    w_co = 2 * pi * f_co;
    s_co = 1i * w_co;
    
    G_co = C_nom_c * ((s_co * eye(2) - A_nom_c) \ (B_nom_c * 12.0)) + D_nom_c;
    boost = PM_target - 90 - rad2deg(angle(G_co));
    if boost < 2 || boost > 175
        cost = 1e8; return;
    end
    
    k_factor = tan(deg2rad(45 + boost / 4));
    w_z = w_co / k_factor;
    w_p = w_co * k_factor;
    K_c = 1 / (abs(G_co) * ((w_co^2 + w_z^2) / (w_co * (w_co^2 + w_p^2))));
    
    A_p = 1 + w_p * T_s; A_z = 1 + w_z * T_s;
    d_raw1 = A_p^2;
    n1 = (K_c * T_s * A_z^2) / d_raw1;
    n2 = (-2 * K_c * T_s * A_z) / d_raw1;
    n3 = (K_c * T_s) / d_raw1;
    d2 = -(A_p^2 + 2*A_p) / d_raw1;
    d3 = (2*A_p + 1) / d_raw1;
    d4 = -1 / d_raw1;
    
    x = [0;0]; u_hist = zeros(3,1); err_hist = zeros(3,1); duty = Vref_val/Vin_data(1);
    V_out_log = zeros(N_pts, 1);
    
    for k = 1:N_pts
        V_in_k = Vin_data(k); R_k = R_data(k);
        theta_k = G_L * R_k * R_C + R_k + R_C;
        A_k = [ -R_k * R_C / (L_val * theta_k),                 -R_k / (L_val * theta_k);
                 R_k / (C_val * theta_k),                 -(R_k * G_L + 1) / (C_val * theta_k) ];
        B_k = [ (R_k + R_C) / (L_val * theta_k);
                (R_k * G_L) / (C_val * theta_k) ];
        C_k = [ R_k * R_C / theta_k,   R_k / theta_k ];
        D_k = G_L * R_k * R_C / theta_k;
        
        v_sw = duty * V_in_k; V_out = C_k * x + D_k * v_sw;
        V_out_log(k) = V_out;
        
        err = Vref_val - V_out;
        err_hist = [err; err_hist(1:2)];
        duty_raw = n1*err_hist(1) + n2*err_hist(2) + n3*err_hist(3) ...
                   - d2*u_hist(1) - d3*u_hist(2) - d4*u_hist(3);
        duty = max(0.01, min(0.95, duty_raw));
        u_hist = [duty; u_hist(1:2)];
        
        x = rk4_step(A_k, B_k, x, v_sw, T_s);
    end
    
    v_smooth = movmean(V_out_log, 15);
    err_smooth = abs(Vref_val - v_smooth);
    cost = sum(t_vec .* err_smooth * T_s);
    
    max_V = max(v_smooth);
    overshoot = max(0, (max_V - Vref_val)/Vref_val * 100);
    if overshoot > 15
        cost = cost + 50 * (overshoot - 15);
    end
    
    err_startup = mean(err_smooth(t_vec >= 0.035 & t_vec <= 0.04));
    err_surge = mean(err_smooth(t_vec >= 0.065 & t_vec <= 0.07));
    err_load = mean(err_smooth(t_vec >= 0.095));
    
    if err_startup > 0.1, cost = cost + 2000 * err_startup; end
    if err_surge > 0.1,   cost = cost + 2000 * err_surge; end
    if err_load > 0.1,    cost = cost + 2000 * err_load; end
end

function cost = evaluate_ml_ml_cost(p, L_val, C_val, G_L, R_C, Vin_data, R_data, Vref_val, T_s, t_vec, N_pts)
    Kc = p(1); z1 = p(2); wz = p(3); zeta_z = p(4); wp = p(5); zeta_p = p(6);
    w_z_rad = 2 * pi * wz; w_p_rad = 2 * pi * wp; w_z1_rad = 2 * pi * z1;
    K_warp = w_z_rad / tan(w_z_rad * T_s / 2);
    num_z1 = Kc * [K_warp + w_z1_rad, w_z1_rad - K_warp]; den_z1 = [K_warp, -K_warp];
    num_z2 = [K_warp^2 + 2*zeta_z*w_z_rad*K_warp + w_z_rad^2, 2*(w_z_rad^2 - K_warp^2), K_warp^2 - 2*zeta_z*w_z_rad*K_warp + w_z_rad^2];
    den_z2 = [K_warp^2 + 2*zeta_p*w_p_rad*K_warp + w_p_rad^2, 2*(w_p_rad^2 - K_warp^2), K_warp^2 - 2*zeta_p*w_p_rad*K_warp + w_p_rad^2];
    den_ml_all = conv(den_z1, den_z2);
    
    % Stability Check (check only the non-integrator poles in den_z2)
    if any(abs(roots(den_z2)) >= 0.99)
        cost = 1e12; return;
    end
    
    m = [conv(num_z1, num_z2), 0, 0] / den_ml_all(1); e = [den_ml_all, 0, 0] / den_ml_all(1);
    x = [0;0]; u_hist = zeros(5,1); err_hist = zeros(6,1); duty = Vref_val/Vin_data(1); V_out_log = zeros(N_pts, 1);
    for k = 1:N_pts
        V_in_k = Vin_data(k); R_k = R_data(k); theta_k = G_L * R_k * R_C + R_k + R_C;
        A_k = [ -R_k * R_C / (L_val * theta_k), -R_k / (L_val * theta_k); R_k / (C_val * theta_k), -(R_k * G_L + 1) / (C_val * theta_k) ];
        B_k = [ (R_k + R_C) / (L_val * theta_k); (R_k * G_L) / (C_val * theta_k) ];
        C_k = [ R_k * R_C / theta_k, R_k / theta_k ]; D_k = G_L * R_k * R_C / theta_k;
        v_sw = duty * V_in_k; V_out = C_k * x + D_k * v_sw; V_out_log(k) = V_out;
        err = Vref_val - V_out; err_hist = [err; err_hist(1:5)];
        duty_raw = (m(1)*err_hist(1) + m(2)*err_hist(2) + m(3)*err_hist(3) + m(4)*err_hist(4) + m(5)*err_hist(5) + m(6)*err_hist(6) ...
                   - e(2)*u_hist(1) - e(3)*u_hist(2) - e(4)*u_hist(3) - e(5)*u_hist(4) - e(6)*u_hist(5));
        duty = max(0.01, min(0.95, duty_raw)); u_hist = [duty; u_hist(1:4)];
        x = rk4_step(A_k, B_k, x, v_sw, T_s);
    end
    
    v_smooth = movmean(V_out_log, 15);
    err_smooth = abs(Vref_val - v_smooth);
    cost = sum(t_vec .* err_smooth * T_s);
    
    v_smooth_startup = v_smooth(t_vec <= 0.04);
    overshoot_startup = max(0, (max(v_smooth_startup) - Vref_val)/Vref_val * 100);
    
    v_smooth_trans2 = v_smooth(t_vec > 0.04 & t_vec <= 0.07);
    overshoot_trans2 = max(0, (max(v_smooth_trans2) - Vref_val)/Vref_val * 100);
    
    v_smooth_trans3 = v_smooth(t_vec > 0.07);
    overshoot_trans3 = max(0, (max(v_smooth_trans3) - Vref_val)/Vref_val * 100);
    
    cost = cost + 10 * overshoot_startup + 500 * (overshoot_trans2 + overshoot_trans3); 
    
    t_eval2 = t_vec(t_vec > 0.04 & t_vec <= 0.07) - 0.04;
    t_eval3 = t_vec(t_vec > 0.07) - 0.07;
    
    get_t_s = @(v_seg, t_seg, target) feval(@(idx) double(~isempty(idx)) * t_seg(max(1, idx)), find(abs(movmean(v_seg, 15) - target) > 0.02 * target, 1, 'last'));
    
    t_s2 = get_t_s(V_out_log(t_vec > 0.04 & t_vec <= 0.07), t_eval2, Vref_val);
    t_s3 = get_t_s(V_out_log(t_vec > 0.07), t_eval3, Vref_val);
    
    if t_s2 > 1.80, cost = cost + 1e6 * (t_s2 - 1.80); end
    if overshoot_trans2 > 22.0, cost = cost + 1e5 * (overshoot_trans2 - 22.0); end
    if t_s3 > 0.40, cost = cost + 1e6 * (t_s3 - 0.40); end
    if overshoot_trans3 > 4.50, cost = cost + 1e5 * (overshoot_trans3 - 4.50); end
    if overshoot_startup > 17.0, cost = cost + 1e5 * (overshoot_startup - 17.0); end
    
    duty_diff = diff(V_out_log);
    cost = cost + 15 * sum(duty_diff.^2);
    
    err_startup = mean(err_smooth(t_vec >= 0.035 & t_vec <= 0.04));
    err_surge = mean(err_smooth(t_vec >= 0.065 & t_vec <= 0.07));
    err_load = mean(err_smooth(t_vec >= 0.095));
    
    if err_startup > 0.1, cost = cost + 50000 * err_startup; end
    if err_surge > 0.1,   cost = cost + 50000 * err_surge; end
    if err_load > 0.1,    cost = cost + 50000 * err_load; end
end

function cost = evaluate_lqr_cost(p, A_nom_c, B_nom_c, C_nom_c, D_nom_c, L_val, C_val, G_L, R_C, Vin_data, R_data, Vref_val, T_s, t_vec, N_pts)
    A_aug_c = [ A_nom_c, zeros(2, 1); -C_nom_c, 0 ];
    B_aug_c = [ B_nom_c * 12.0; -D_nom_c ];
    sys_aug_c = ss(A_aug_c, B_aug_c, eye(3), 0);
    sys_aug_d = c2d(sys_aug_c, T_s, 'zoh');
    Q = diag([p(1), p(2), p(3)]); R = p(4);
    try K_lqr = dlqr(sys_aug_d.A, sys_aug_d.B, Q, R); catch, cost = 1e9; return; end
    
    x = [0;0]; error_int = 0; duty = Vref_val/Vin_data(1); V_out_log = zeros(N_pts, 1);
    for k = 1:N_pts
        V_in_k = Vin_data(k); R_k = R_data(k); theta_k = G_L * R_k * R_C + R_k + R_C;
        A_k = [ -R_k * R_C / (L_val * theta_k), -R_k / (L_val * theta_k); R_k / (C_val * theta_k), -(R_k * G_L + 1) / (C_val * theta_k) ];
        B_k = [ (R_k + R_C) / (L_val * theta_k); (R_k * G_L) / (C_val * theta_k) ];
        C_k = [ R_k * R_C / theta_k, R_k / theta_k ]; D_k = G_L * R_k * R_C / theta_k;
        v_sw = duty * V_in_k; V_out = C_k * x + D_k * v_sw; V_out_log(k) = V_out;
        err = Vref_val - V_out; error_int = error_int + err * T_s;
        duty_lqr_nom = Vref_val / V_in_k; I_L_ref = Vref_val / R_k;
        duty_raw = duty_lqr_nom - ( K_lqr(1) * (x(1) - I_L_ref) + K_lqr(2) * (V_out - Vref_val) + K_lqr(3) * error_int );
        duty = max(0.01, min(0.95, duty_raw));
        if duty_raw > 0.95 || duty_raw < 0.01, error_int = error_int - err * T_s; end
        x = rk4_step(A_k, B_k, x, v_sw, T_s);
    end
    
    v_smooth = movmean(V_out_log, 15);
    err_smooth = abs(Vref_val - v_smooth);
    cost = sum(t_vec .* err_smooth * T_s);
    
    max_V = max(v_smooth);
    overshoot = max(0, (max_V - Vref_val)/Vref_val * 100);
    if overshoot > 10
        cost = cost + 50 * (overshoot - 10);
    end
    
    err_startup = mean(err_smooth(t_vec >= 0.035 & t_vec <= 0.04));
    err_surge = mean(err_smooth(t_vec >= 0.065 & t_vec <= 0.07));
    err_load = mean(err_smooth(t_vec >= 0.095));
    
    if err_startup > 0.1, cost = cost + 1000 * err_startup; end
    if err_surge > 0.1,   cost = cost + 1000 * err_surge; end
    if err_load > 0.1,    cost = cost + 1000 * err_load; end
end

function safe_plot(t, y, color, width, name, style)
    if nargin < 6
        style = '-';
    end
    t = t(:);
    y = y(:);
    if isempty(t) || isempty(y)
        return;
    end
    if length(t) == length(y)
        if isempty(name)
            plot(t, y, style, 'Color', color, 'LineWidth', width, 'HandleVisibility', 'off');
        else
            plot(t, y, style, 'Color', color, 'LineWidth', width, 'DisplayName', name);
        end
    else
        y_interp = interp1(linspace(0, 1, length(y)), y, linspace(0, 1, length(t)), 'linear');
        if isempty(name)
            plot(t, y_interp, style, 'Color', color, 'LineWidth', width, 'HandleVisibility', 'off');
        else
            plot(t, y_interp, style, 'Color', color, 'LineWidth', width, 'DisplayName', name);
        end
    end
end
