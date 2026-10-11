"""MATLAB/Simulink export adapter for conductor state-space models."""
from pathlib import Path
import numpy as np

def save_mat(out_path, *, A, B, C, D, V, S, sigma, mu, K_SIBC, c_1,
              tau_foster, g_foster, n_foster, n_warburg, n_integrator,
              description="Conductor Y(s) Foster pole-residue LTI state-space."):
    """Save state-space + metadata to MATLAB-loadable .mat.

    Also writes a sibling .m helper script for Simulink loading.
    """
    from scipy.io import savemat
    out_path = Path(out_path)
    savemat(str(out_path), {
        "A": A, "B": B, "C": C, "D": D,
        "V_conductor": V, "S_conductor": S,
        "sigma": sigma, "mu": mu,
        "K_SIBC": K_SIBC, "c_1": c_1,
        "tau_foster": np.asarray(tau_foster),
        "g_foster": np.asarray(g_foster),
        "n_foster": n_foster, "n_warburg": n_warburg,
        "n_integrator": n_integrator,
        "description": description,
    })

    m_script = out_path.with_suffix(".m")
    with open(m_script, "w", encoding="utf-8") as fp:
        fp.write(f"""% {out_path.name} loader for Simulink
load('{out_path.name}');

% State-space LTI
sys = ss(A, B, C, D);
disp(sys);
fprintf('  states: %d (%d Foster + %d Warburg + %d integrator)\\n', ...
        size(A, 1), n_foster, n_warburg, n_integrator);
fprintf('  conductor: V = %.2f mm^3, S = %.2f mm^2, sigma = %.2e\\n', ...
        V_conductor*1e9, S_conductor*1e6, sigma);

figure;
bodeplot(sys);
title('Conductor Y(s) frequency response');

% Time-domain step response example
T = 0:1e-6:5e-3;
u_step = 0.1 * (T > 1e-3);
du_step = gradient(u_step, T(2)-T(1));
[y_out, t_out] = lsim(sys, du_step, T);
figure; plot(t_out, y_out);
xlabel('time (s)'); ylabel('Y(t) admittance');
title('Step response');
""")
    return str(out_path), str(m_script)

