"""
Extract panel parameter definitions from calc_*.py argparse and generate panel_registry.json.

Usage:
    python sync_registry.py              # Generate panel_registry.json
    python sync_registry.py --check      # Verify registry matches code

This is the "reverse" direction: code -> registry.
The "forward" direction (registry -> code) is handled by the LLM via MCP tools.
"""

import argparse
import ast
import json
import os
import sys


_this_dir = os.path.dirname(os.path.abspath(__file__))

# Japanese labels and physics descriptions for each parameter.
# This is the domain knowledge that enables natural language <-> CLI mapping.
PARAM_JA = {
    # --- Common ---
    "vol": {"ja": ".vol メッシュファイル", "physics": "Netgen Vol 形式。材料・境界ラベル・高次要素を含む"},
    "vol_file": {"ja": ".vol メッシュファイル", "physics": "Netgen Vol 形式"},
    "output": {"ja": "出力 JSON ファイル", "physics": ""},
    "msh_output": {"ja": "GMSH 出力パス", "physics": "場の可視化用 .msh v4.1"},
    "msh-output": {"ja": "GMSH 出力パス", "physics": "場の可視化用 .msh v4.1"},

    # --- Frequency / Material ---
    "frequency": {"ja": "周波数", "physics": "渦電流の浸透深さ delta = sqrt(2*rho/(omega*mu)) を決める"},
    "sigma": {"ja": "導電率", "physics": "材料の電気伝導率。銅: 5.8e7, 鉄: 2e6 [S/m]"},
    "mu_r": {"ja": "比透磁率", "physics": "mu = mu_0 * mu_r。非磁性体: 1, 鉄: 100-5000"},
    "mu-r": {"ja": "比透磁率", "physics": "mu = mu_0 * mu_r"},
    "H_amplitude": {"ja": "固定子等価磁界振幅", "physics": "回転子に印加する一様背景磁界 H [A/m]"},
    "field_angle_deg": {"ja": "背景磁界角", "physics": "固定した全体座標系での背景磁界方向 [deg]"},
    "rotor_angle_start_deg": {"ja": "回転子開始角", "physics": "機械角スイープの開始角 [deg]"},
    "rotor_angle_stop_deg": {"ja": "回転子終了角", "physics": "機械角スイープの終了角 [deg]"},
    "rotor_angle_steps": {"ja": "回転角点数", "physics": "開始角と終了角を含む機械角評価点数"},
    "maxwell_radius": {"ja": "Maxwell応力評価半径", "physics": "回転子を囲み固定子を含まない空隙円の半径 [m]"},
    "stack_length": {"ja": "積厚", "physics": "2D単位長さトルクを実機トルクへ変換する軸方向長さ [m]"},
    "energy_delta_deg": {"ja": "共エネルギー差分角", "physics": "固定電流共エネルギーの中心差分半幅 [deg]"},
    "circle_points": {"ja": "空隙円積分点数", "physics": "Maxwell応力円周積分の等間隔点数"},
    "center_x": {"ja": "回転中心X", "physics": "回転子と空隙円の共通回転中心X [m]"},
    "center_y": {"ja": "回転中心Y", "physics": "回転子と空隙円の共通回転中心Y [m]"},
    "eta": {"ja": "H行列許容係数", "physics": "Charge Gram H行列の幾何学的admissibility係数"},
    "material": {"ja": "材料プリセット", "physics": "steel/copper/aluminum の物性値を自動設定"},
    "bh_file": {"ja": "BH カーブファイル", "physics": "2列: H[A/m] B[T]。非線形透磁率"},
    "bh-file": {"ja": "BH カーブファイル", "physics": "2列: H[A/m] B[T]"},

    # --- BEM (calc_inductance) ---
    "source": {"ja": "電流注入面", "physics": "BEM で単位電流を注入する境界面"},
    "sink": {"ja": "電流引出面", "physics": "BEM で電流を引き出す境界面"},
    "fes_order": {"ja": "基底関数次数", "physics": "0=RWG (BEM), 1-3=高次 (FEM HCurl)"},
    "fes-order": {"ja": "基底関数次数", "physics": "0=RWG, 1+=高次"},
    "coil": {"ja": "コイル材料名", "physics": "BEM をコイル表面のみに限定するフィルタ"},
    "coil_sigma": {"ja": "コイル導電率", "physics": "コイル抵抗 R = L/(sigma*A) の計算用"},
    "coil-sigma": {"ja": "コイル導電率", "physics": "コイル抵抗計算用 [S/m]"},
    "workpiece": {"ja": "ワークピース境界", "physics": "加熱対象の表面ラベル (prefix match)"},
    "wp_loop_dof": {"ja": "ループ電流自由度", "physics": "種数1ワークピースの短絡周回電流とLenz遮蔽を解く"},
    "wp-loop-dof": {"ja": "ループ電流自由度", "physics": "種数1ワークピースの短絡周回電流とLenz遮蔽を解く"},

    # --- Impedance model ---
    "impedance_model": {"ja": "インピーダンスモデル", "physics": "esim=非線形セル問題, dowell=解析式, sibc=線形"},
    "impedance-model": {"ja": "インピーダンスモデル", "physics": "esim/dowell/sibc"},
    "impedance": {"ja": "インピーダンスモデル", "physics": "sibc=線形SIBC, esim=非線形ESIM"},
    "half_thickness": {"ja": "半厚さ", "physics": "ESIM セル問題の厚さ方向サイズ [m]"},
    "half-thickness": {"ja": "半厚さ", "physics": "ESIM 厚さ [m]"},
    "esim_geometry": {"ja": "ESIM 形状", "physics": "cylinder=ベッセル関数, slab=双曲線関数"},
    "esim-geometry": {"ja": "ESIM 形状", "physics": "cylinder/slab"},

    # --- FEM (calc_fem_kelvin) ---
    "current": {"ja": "コイル電流", "physics": "ソース電流 [A]。J = I / A_cross で電流密度に変換"},
    "a_coil": {"ja": "コイル線半径", "physics": "J_theta フォールバック用 [m]"},
    "a-coil": {"ja": "コイル線半径", "physics": "[m]"},
    "r_wp": {"ja": "ワークピース半径", "physics": "ESIM 半厚さの推定に使用 [m]"},
    "r-wp": {"ja": "ワークピース半径", "physics": "[m]"},
    "max_iter": {"ja": "最大反復回数", "physics": "Karl 反復 (ESIM 非線形) の上限"},
    "max-iter": {"ja": "最大反復回数", "physics": ""},
    "solver": {"ja": "線形ソルバー", "physics": "pardiso=直接法, bddc=反復法, ams=Compact AMS+COCR"},
    "reg": {"ja": "ゲージ正則化", "physics": "curl-curl 系のカーネル除去。reg*nu0*u*v*dx を追加"},
    "shift_eps": {"ja": "シフト前処理", "physics": "AMS/ICCG 用。前処理のみに eps*nu*u*v*dx を追加"},
    "shift-eps": {"ja": "シフト前処理", "physics": "AMS/ICCG 用"},
    "nthreads": {"ja": "スレッド数", "physics": "TaskManager 並列数。0=自動, 1=シングル"},

    # --- BEM-SIBC (calc_heating_bem) ---
    "coil_radius": {"ja": "コイル半径", "physics": "Biot-Savart フィラメントコイル [m]"},
    "coil-radius": {"ja": "コイル半径", "physics": "[m]"},
    "coil_current": {"ja": "コイル電流", "physics": "[A]"},
    "coil-current": {"ja": "コイル電流", "physics": "[A]"},
    "gap_deg": {"ja": "ギャップ角度", "physics": "コイルの開口角 [deg]"},
    "gap-deg": {"ja": "ギャップ角度", "physics": "[deg]"},
    "h1_order": {"ja": "H1 次数", "physics": "スカラー BIE の基底関数次数"},
    "h1-order": {"ja": "H1 次数", "physics": ""},
    "wp_label": {"ja": "WP 境界ラベル", "physics": "ワークピース表面のラベル名"},
    "wp-label": {"ja": "WP 境界ラベル", "physics": ""},
    "tol": {"ja": "収束判定", "physics": "Karl 反復の Z_s 変化率閾値"},

    # --- Accelerator (calc_accel_magnet) ---
    "coil_script": {"ja": "コイル定義スクリプト", "physics": "build_coil() -> CoilBuilder を含む .py"},
    "coil-script": {"ja": "コイル定義スクリプト", "physics": ""},
    "formulation": {"ja": "定式化", "physics": "omega=スカラーH1, a=ベクトルHCurl"},
    "hys_file": {"ja": "ヒステリシスファイル", "physics": "Play モデル .hys"},
    "hys-file": {"ja": "ヒステリシスファイル", "physics": ""},
    "n_steps": {"ja": "準静的ステップ数", "physics": "ヒステリシス: 0->I->0 の分割数"},
    "n-steps": {"ja": "準静的ステップ数", "physics": ""},
    "newton": {"ja": "Newton 法", "physics": "True=Newton-Raphson, False=Picard 反復"},
    "relax": {"ja": "緩和係数", "physics": "0=フルステップ, 0.3=30%減衰"},

    # --- HDiv-VIM (calc_accel_hdiv) ---
    "ima": {"ja": "IMA 対称性", "physics": "'+x-z'=1/4モデル。符号: 平行=+, 垂直=-"},
    "hdiv_order": {"ja": "HDiv 次数", "physics": "1=BDM1, 2=BDM2"},
    "unit_scale": {"ja": "座標スケール", "physics": "Cubit座標→メートル。0.001=mm→m"},
    "unit-scale": {"ja": "座標スケール", "physics": ""},

    # --- PEEC (calc_pcb_peec) ---
    "inp": {"ja": "FastHenry 入力", "physics": ".inp ファイル (ノード・セグメント定義)"},
    "freq_min": {"ja": "最小周波数", "physics": "周波数スイープ下限 [Hz]"},
    "freq-min": {"ja": "最小周波数", "physics": "[Hz]"},
    "freq_max": {"ja": "最大周波数", "physics": "周波数スイープ上限 [Hz]"},
    "freq-max": {"ja": "最大周波数", "physics": "[Hz]"},
    "n_freq": {"ja": "周波数点数", "physics": "対数等間隔"},
    "n-freq": {"ja": "周波数点数", "physics": ""},
    "solver_method": {"ja": "ソルバー", "physics": "0=LU, 1=BiCGSTAB, 2=HACApK"},
    "solver-method": {"ja": "ソルバー", "physics": ""},
    "spice_output": {"ja": "SPICE 出力", "physics": "ネットリスト .sp ファイル"},
    "spice-output": {"ja": "SPICE 出力", "physics": ""},

    # --- Volume / Mesh ---
    "order": {"ja": "曲線次数", "physics": "メッシュ高次要素次数 (1-5)"},
    "maxh": {"ja": "最大メッシュサイズ", "physics": "[m]"},
    "max_order": {"ja": "最大次数", "physics": "p-convergence の上限"},
    "max-order": {"ja": "最大次数", "physics": ""},
    "cub5": {"ja": "Cubit モデル", "physics": ".cub5 ファイル (メッシュ生成に Cubit 必要)"},
    "step": {"ja": "STEP 形状", "physics": "OCC 経由で Netgen メッシュ生成"},

    # --- Stream-function volume (calc_streamfunction_volume) ---
    "target_bz": {"ja": "目標軸方向磁場 Bz", "physics": "ボア中心の一様 Bz [T]"},
    "target-bz": {"ja": "目標軸方向磁場 Bz", "physics": "[T]"},
    "foliation": {"ja": "葉層化 mu", "physics": "radial = 入れ子円筒 mu=r"},
    "n_leaves": {"ja": "mu 葉層数", "physics": "径方向シェル数 (円筒葉)"},
    "n-leaves": {"ja": "mu 葉層数", "physics": ""},
    "fes_order": {"ja": "lambda H1 次数", "physics": "体積ストリーム関数 lambda の H1 多項式次数"},
    "fes-order": {"ja": "lambda H1 次数", "physics": ""},
    "aca_eps": {"ja": "ACA+ 許容誤差", "physics": "応答行列 ACA+ 分解の停止許容誤差"},
    "aca-eps": {"ja": "ACA+ 許容誤差", "physics": ""},
    "rings_per_span": {"ja": "リング数/lambda 範囲", "physics": "等電流ワイヤ密度 (大きいほど細かい)"},
    "rings-per-span": {"ja": "リング数/lambda 範囲", "physics": ""},
    "n_targets": {"ja": "軸方向ターゲット点数", "physics": "ボア軸上の磁場制約点数"},
    "n-targets": {"ja": "軸方向ターゲット点数", "physics": ""},
    "nphi": {"ja": "等高線 方位サンプル数", "physics": "ワイヤ抽出の方位方向サンプル数"},
    "nz": {"ja": "等高線 軸方向サンプル数", "physics": "ワイヤ抽出の軸方向サンプル数"},
}

# Parameters such as ``--order`` have panel-specific meanings.  Apply these
# overrides after the generic CLI-name metadata above.
PANEL_PARAM_JA = {
    "motor_hdiv_reduced": {
        "order": {"ja": "HDiv 次数", "physics": "1=BDM1, 2=BDM2"},
    },
}


# Panel definitions: script -> panel metadata
PANELS = {
    "inductance": {
        "script": "calc_inductance.py",
        "function": "run_inductance",
        "ja_name": "Inductance (vacuum + weak-coupled BEM-SIBC)",
        "ja_description": "コイル STEP -> L_coil (vacuum) または L_coil + ΔL (weak-coupled).  --coil-solver peec | bem-a で切替。--vol 指定で workpiece 弱結合 (Telegen φ・B back-reaction)。",
        "method": "PEEC perimeter filaments OR BEM-A Weggler EFIE saddle (RWG); workpiece scalar BIE + SIBC Robin",
        "command_builder": "radia.ih_design:IHDesignSpec.build_command",
    },
    "fem_coilmesh": {
        "script": "calc_fem_coilmesh.py",
        "function": "solve_fem_coilmesh",
        "ja_name": "FEM A-V (IH、コイル体積メッシュ)",
        "ja_description": "A-V 定式化 (HCurl + H1 phi on coil) + source/sink Dirichlet lift + 体積積分電流抽出 + wp SIBC Robin + Kelvin。gapped torus 専用。L, P_wp, P_coil を出力。",
        "method": "FEM A-V (HCurl A + H1 phi on coil, phi Dirichlet source/sink)",
        "command_builder": "radia.ih_design:IHDesignSpec.build_command",
    },
    "volume": {
        "script": "calc_volume.py",
        "function": "calculate_volume_vol",
        "ja_name": "体積・面積積分",
        "ja_description": "材料別体積、境界別面積を NGSolve Integrate で計算",
        "method": "NGSolve Integrate(CF(1), mesh)",
        "command_builder": None,
    },
    "accel_magnet": {
        "script": "calc_accel_magnet.py",
        "function": "solve_accel",
        "ja_name": "加速器電磁石 FEM",
        "ja_description": "Omega-reduced/A 定式化 + Kelvin + BH 非線形 + ヒステリシス",
        "method": "FEM (H1 Omega / HCurl A-field)",
        "command_builder": "radia.em_design:EMDesignSpec.build_command",
    },
    "accel_hdiv": {
        "script": "calc_accel_hdiv.py",
        "function": "solve_hdiv",
        "ja_name": "加速器電磁石 HDiv-VIM",
        "ja_description": "HDiv-VIM + IMA ラベル読取 + HACApK charge Gram",
        "method": "FEEC HDiv-VIM (NGSolve-aligned mesh-backed soft iron)",
        "command_builder": "radia.em_design:EMDesignSpec.build_command",
    },
    "pcb_peec": {
        "script": "calc_pcb_peec.py",
        "function": "solve_peec",
        "ja_name": "PCB PEEC インピーダンス",
        "ja_description": "FastHenry 入力 → PEEC 行列 → 周波数スイープ → SPICE 出力",
        "method": "PEEC (Loop-Star + MNA)",
        "command_builder": "radia.pcb_design:PCBDesignSpec.build_command",
    },
    "streamfunction": {
        "script": "calc_streamfunction.py",
        "function": "run",
        "ja_name": "ストリーム関数コイル設計",
        "ja_description": "コイル面 .vol + 評価領域 .vol + 目標磁場 (CoefficientFunction 式; スカラー→Bz / 3成分→ベクトルB) から表面電流ストリーム関数 ψ を設計。design (均一度・最大電流密度) / pareto (フロント α掃引) / manufacture (一筆書き+板金+STEP+PEEC) の3モード。",
        "method": "FE-direct H1 psi on a surface .vol + surface Biot-Savart + folded-Tikhonov (ACA+TSVD RegularizedTSVD)",
        "command_builder": "radia.streamfunction_design:StreamFunctionDesignSpec.build_command",
    },
    "streamfunction_volume": {
        "script": "calc_streamfunction_volume.py",
        "function": "run",
        "ja_name": "体積ストリーム関数コイル設計 (3D)",
        "ja_description": "中空円筒導体 .vol + 目標軸方向 Bz から、葉層化カレント・クレブシュ体積ストリーム関数 J=grad(lambda)xgrad(mu) (mu=r 固定) を線形最小ノルム逆問題 (ACA+TSVD) で解き、等電流の巻線ワイヤに等高線抽出。ワイヤの Biot-Savart を Radia と二重コード照合。radia_streamfunction パネルの 'Volume 3D' モード。",
        "method": "Foliated current-Clebsch volume stream function (H1 lambda on a conductor volume) + ACA+TSVD least-norm + equal-current contour extraction + two-codebase Biot-Savart check",
        "command_builder": "radia.streamfunction_design:StreamFunctionDesignSpec.build_command",
    },
    "motor_transient": {
        "script": "calc_motor_transient.py",
        "function": "run",
        "ja_name": "モータ過渡解析 (Lange-Henrotte-Hameyer)",
        "ja_description": "モータ .vol から、非線形 FE + 回路 ODE 連成で過渡応答を解く。PM 回転と Arkkio トルクを含む。radia_motor パネルの 'Transient' モード。",
        "method": "Nonlinear FE (A-formulation) + circuit ODE co-simulation, PM rotation, Arkkio torque (Lange-Henrotte-Hameyer 2009)",
        "command_builder": "radia.motor_design:MotorDesignSpec.build_command",
    },
    "motor_lamination": {
        "script": "calc_motor_lamination.py",
        "function": "run",
        "ja_name": "積層鉄心 等価材料 (Hollaus)",
        "ja_description": "積層鉄心の 1D セル問題で等価材料テーブルを構築し (cell)、それを .vol 上の global FE に適用 (global/full)。radia_motor パネルの 'Lamination' モード。",
        "method": "Hollaus effective-material homogenization (1D cell problem) + global HCurl FE with the effective-material table",
        "command_builder": "radia.motor_design:MotorDesignSpec.build_command",
    },
    "motor_hdiv_reduced": {
        "script": "calc_motor_hdiv_reduced.py",
        "function": "run",
        "ja_name": "HDiv-VIM 縮約リラクタンスモータ",
        "ja_description": "回転子のみの2D .volで対称HDiv-VIM Gramを一度だけ構築し、固定子等価MMFに対する回転角スイープ、Maxwell応力、磁化体積トルク、共エネルギー差分を照合する。radia_motorパネルの'HDiv Reduced'モード。",
        "method": "Planar BDM1 HDiv-VIM reduced reluctance motor with cached symmetric charge Gram and three-way torque check",
        "command_builder": "radia.motor_design:MotorDesignSpec.build_command",
    },
}


def _extract_argparse_from_file(filepath):
    """Parse a calc_*.py file and extract argparse add_argument calls.

    Uses AST parsing to find parser.add_argument(...) calls and extract
    the argument name, type, default, help, choices, and required flag.
    """
    with open(filepath, "r", encoding="utf-8") as f:
        source = f.read()

    tree = ast.parse(source)
    args = []

    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        # Match parser.add_argument(...) or p.add_argument(...)
        func = node.func
        if not (isinstance(func, ast.Attribute) and func.attr == "add_argument"):
            continue

        # First positional arg is the flag name
        if not node.args:
            continue
        flag = ast.literal_eval(node.args[0])
        if not flag.startswith("--"):
            continue

        param_name = flag.lstrip("-").replace("-", "_")
        info = {"cli": flag, "name": param_name}

        for kw in node.keywords:
            key = kw.arg
            try:
                val = ast.literal_eval(kw.value)
            except (ValueError, TypeError):
                # Non-literal (e.g., function reference for type=)
                if key == "type":
                    if isinstance(kw.value, ast.Name):
                        val = kw.value.id  # "int", "float", "str"
                    else:
                        val = "str"
                else:
                    continue

            if key in ("default", "help", "choices", "required", "type"):
                info[key] = val

        # Add Japanese label and physics if available
        for lookup_key in [flag.lstrip("-"), param_name]:
            if lookup_key in PARAM_JA:
                info["ja"] = PARAM_JA[lookup_key]["ja"]
                info["physics"] = PARAM_JA[lookup_key]["physics"]
                break

        args.append(info)

    return args


def build_registry():
    """Build the complete panel registry from code."""
    registry = {"_meta": {
        "generated_by": "sync_registry.py",
        "description": "Panel parameter registry for Radia-NGSolve. "
                       "Source of truth for LLM <-> panel mapping.",
        "ja_description": "Radia-NGSolve パネル定義。"
                          "自然言語 <-> CLI パラメータの双方向変換に使用。",
    }, "panels": {}}

    for panel_id, meta in PANELS.items():
        script_path = os.path.join(_this_dir, meta["script"])
        if not os.path.exists(script_path):
            print(f"SKIP: {meta['script']} not found", file=sys.stderr)
            continue

        params = _extract_argparse_from_file(script_path)
        panel_param_ja = PANEL_PARAM_JA.get(panel_id, {})
        for param in params:
            override = panel_param_ja.get(param["name"])
            if override is not None:
                param.update(override)

        registry["panels"][panel_id] = {
            "script": meta["script"],
            "function": meta["function"],
            "ja_name": meta["ja_name"],
            "ja_description": meta["ja_description"],
            "method": meta["method"],
            "command_builder": meta["command_builder"],
            "params": params,
        }

    return registry


def main():
    parser = argparse.ArgumentParser(
        description="Generate panel_registry.json from calc_*.py argparse")
    parser.add_argument("--check", action="store_true",
                        help="Check existing registry matches code")
    parser.add_argument("--output", default="",
                        help="Output path (default: panel_registry.json)")
    args = parser.parse_args()

    registry = build_registry()

    output_path = args.output or os.path.join(_this_dir, "panel_registry.json")

    if args.check:
        if not os.path.exists(output_path):
            print(f"FAIL: {output_path} does not exist")
            sys.exit(1)
        with open(output_path, "r", encoding="utf-8") as f:
            existing = json.load(f)
        # Compare panels (ignore _meta)
        for pid in registry["panels"]:
            if pid not in existing.get("panels", {}):
                print(f"FAIL: panel '{pid}' not in registry")
                sys.exit(1)
            code_params = {p["cli"] for p in registry["panels"][pid]["params"]}
            reg_params = {p["cli"] for p in existing["panels"][pid]["params"]}
            if code_params != reg_params:
                added = code_params - reg_params
                removed = reg_params - code_params
                if added:
                    print(f"FAIL: {pid}: params in code but not registry: {added}")
                if removed:
                    print(f"FAIL: {pid}: params in registry but not code: {removed}")
                sys.exit(1)
        print("OK: registry matches code")
        sys.exit(0)

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(registry, f, indent=2, ensure_ascii=False)
    print(f"Generated: {output_path}")

    # Summary
    for pid, panel in registry["panels"].items():
        n = len(panel["params"])
        print(f"  {pid}: {panel['ja_name']} ({n} params)")


if __name__ == "__main__":
    main()
