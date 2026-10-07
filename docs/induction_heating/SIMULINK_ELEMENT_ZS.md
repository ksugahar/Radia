# Simulinkで要素別の表面インピーダンスを指定する

P1・弱連成の閉じたワーク表面に、三角形ごとの固定した複素Zs
（単位Ω）を指定できます。係数は面の重み付き剛性と局所発熱の両方に
使います。平均Zsで電磁場を解く方法ではありません。
genus-0に加え、対応する一つのz軸貫通穴では循環電流を含めて解きます。
穴を持つ場合は平面三角形の密行列経路を使い、最大7000表面頂点です。
P0 行列だけでも約14000面で1.46 GiBを要し、組立時の一時コピーには
その数倍のメモリを見込んでください。8192面の試験では行列512 MiBに対し、
P0 組立プロセスのピークは4.15 GiBでした（通常 BIE のメモリは別途必要）。

## 学生向けの手順

1. `sibc`境界を持つワークのNetgen `.vol`を用意します。周波数を決め、
   次のコマンドで、その表面に対応したJSONのひな形を作ります。

   ```powershell
   python -m radia.surface_impedance C:\temp\ih-case\workpiece.vol --label sibc --frequency 1000 --real-ohm 0.001 --imag-ohm 0.001 --output C:\temp\ih-case\zs-template.json
   ```

2. MATLABで表面三角形の重心から領域を選び、Zsを設定します。例の値は
   入力方法の説明用です。実材料のZsは別途決定・検証してください。

   ```matlab
   table = jsondecode(fileread('C:\temp\ih-case\zs-template.json'));
   Zs = 0.001*(1+1i)*ones(size(table.centroids_m,1),1);
   Zs(table.centroids_m(:,3)>0) = 0.002*(1+1i);
   radia.simulink.writeIHPanelImpedance( ...
       'C:\temp\ih-case\zs-template.json', 'C:\temp\ih-case\zs.json', Zs);
   ```

3. `radia.simulink.openIH()`でIHモデルを開き、**Geometry Update**を開きます。
   ワークとコイル、周波数、材料・熱物性を設定し、**Element Zs JSON**に
   `C:\temp\ih-case\zs.json`を指定します。**Rebuild now**で電磁解と熱源を
   再構築します。コイルSTEPはPEEC、コイル`.vol`はBEM-Aを使います。
   ひな形とモデルの周波数は一致させてください。

4. 出力の設定JSONを保存します。単体IHモデルではそのまま実行できます。
   共通ライブラリの**Induction Heating**ブロックでは、その設定JSONを
   `IHDesignSpec configuration`から読み込み、電流・角度・周囲温度を接続します。
   q_surf、体積熱源の積分電力、温度分布を確認してください。

5. `radia-ih-fem-crosscheck`スキルで同条件の独立FEM解と比較します。
   電力収支だけが合う結果を、物理的に正しいという判断には使いません。

JSONは`mesh.Elements(BND)`順です。`surface_sha256`、単位、周波数は変更せず、
`real_ohm`と`imag_ohm`だけを編集してください。実部は非負、全値は有限に
します。メッシュを変更・再分割したらひな形から作り直します。異なるメッシュ、
配列長、周波数はエラーです。

## 適用範囲と検証

Zsはシミュレーション中に固定します。電流変更に対して磁場は比例、発熱は
電流の二乗に比例する線形応答です。温度や磁場に応じた非線形Zsの逐次更新は
このルートでは行いません。Zsファイルの内容が変われば次の更新で再構築し、
ジョブ内のコピーとSHA-256を記録します。強連成、P2電磁場、複数の穴や
未対応のハンドル形状は拒否します。一つの対応する穴には、平面三角形と
密行列BEMを使う循環電流の拡張が必要です。Htと発熱はその循環成分を含み、
Zsの不連続面でも面ごとの重みを使います。非線形ESIMのgenus-1対応は
この固定値指定とは別です。熱FEMの次数とは別の条件です。

設定JSONの`independent_fem_validation.status=not-performed`は独立FEM比較が
自動実行されていないことを表します。スキルによるケースごとの比較報告を
併せて保存してください。同じZsを使うFEM-SIBC比較は離散化の検証であり、
任意に指定したZsの材料モデルとしての妥当性を保証しません。

指定した要素別Zsが電磁境界条件になります。sigma / mu の入力はこのモードでは
参照情報であり、指定Zsから物理的な表皮深さを逆算しません。出力の
`skin_depth_wp_mm` は null、`impedance_model` は `specified-panel` です。
`writeIHPanelImpedance` は全要素の順序付き複素ベクトルを受け取ります。
座標・材料領域からの要素選択は呼び出し側で `centroids_m` を使って行います。

## 既存経路を通す実行例

`validation_test/induction_heating/run_panel_zs_student_workflow.py` は、既存の
形状・電磁アセンブラを使って小さなワークとコイルを生成し、ひな形 CLI →
`writeIHPanelImpedance` → Geometry Update → native Eddy/Thermal を通します。
公式 MATLAB MCP と Simulink Agentic Toolkit の `model_edit` でマスクを設定し、
モデルはジョブ用ディレクトリに生成します。新しいソルバーや追跡対象の SLX を
追加する例ではありません。

固定した Radia wheel/native runtime、互換 MEX、`cubit_mesh_export`、MATLAB、
公式 MCP/Toolkit、Python MCP SDK と pytest がある LAB で実行します。
既存 MATLAB がある場合は追加起動せず、明示したセッションを再利用する必要があります。
この実行例は MATLAB が存在しない場合だけ、公式 MCP が所有する1セッションを起動します。
各オプションには、そのホストの実在する絶対パスを指定します。
ソースツリーをジョブ内に転送して実行してください。固定 wheel の native 部分を
使い、転送した Python アダプタをジョブ内の bootstrap で読み込みます。
editable install は不要です。MEX は検証対象の native ソースから作成し、
`radia_mex.mexw64.build.json` を同じフォルダに置いてください。

```powershell
$server = Join-Path $env:USERPROFILE '.matlab\agentic-toolkits\bin\matlab-mcp-server.exe'
$tools = Join-Path $env:USERPROFILE '.matlab\agentic-toolkits\simulink\tools\tools.json'
$mex = Read-Host 'Compatible native MATLAB directory (absolute path)'
python validation_test/induction_heating/run_panel_zs_student_workflow.py --job C:\temp\panel-zs-student --matlab-mcp $server --toolkit-tools $tools --mex-dir $mex
```

パス表記の違い、入力不変時の再構築省略と workspace 再読込、保存・再オープン、
Zs 変更時の再構築と snapshot hash、Zs 指定を消した場合の uniform モードを確認します。
solve が記録した要素順の複素 Zs 配列を指定値と照合し、非一様指定時と
uniform 時の電力に差が出ることも確認します。差の符号は判定条件にしません。
成功時の `result.json` と失敗時の `failure.json`、公式ツール呼出しの
`tool-events.json` はジョブディレクトリに保存します。電磁・熱の積分電力保存と
native 出力の有限性は接続の検証です。独立 FEM 比較は別途必要で、設定の
`independent_fem_validation.status` は `not-performed` のままです。
この自前例題の簡潔な検証結果は
`validation_test/induction_heating/results/panel_zs_student_workflow.json` に保存しています。


## Python/CLI の固定面別 Zs・強結合

Python/CLI の `bem-a` と `peec` 強結合は、genus-0/1、平面 P1、
`intree-dense` の workpiece に対して固定 `--panel-zs-file` を使用できる。
各面の値は局所発熱、完全反作用、genus-1 のループ電流に直接用いられる。
面積平均 Zs は報告用であり、計算には代用しない。値が全て完全一致する
配列はスカラー演算に一致し、微小でも値が異なる配列は重み付き経路を使う。
値・周波数の変更時は分解を更新し、同じ形状に対する準備済み幾何演算を再利用する。

PEEC 強結合では `--no-peec-proximity` が必要。非線形 strong ESIM、P2、
面別 Zs の HACApK、複数ハンドルは未対応。この追加は Python/CLI の電磁経路であり、
本書の Simulink 熱連成手順は引き続き弱結合用。strong の native 熱連成認証を意味しない。
