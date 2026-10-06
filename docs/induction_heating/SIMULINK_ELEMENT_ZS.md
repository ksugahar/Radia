# Simulinkで要素別の表面インピーダンスを指定する

P1・弱連成・閉じたgenus-0ワーク表面に、三角形ごとの固定した複素Zs
（単位Ω）を指定できます。係数は面の重み付き剛性と局所発熱の両方に
使います。平均Zsで電磁場を解く方法ではありません。

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
ジョブ内のコピーとSHA-256を記録します。強連成、P2電磁場、ハンドル電流を
要するgenus-1以上の表面では、この指定を拒否します。熱FEMの次数とは別です。

設定JSONの`independent_fem_validation.status=not-performed`は独立FEM比較が
自動実行されていないことを表します。スキルによるケースごとの比較報告を
併せて保存してください。同じZsを使うFEM-SIBC比較は離散化の検証であり、
任意に指定したZsの材料モデルとしての妥当性を保証しません。

指定した要素別Zsが電磁境界条件になります。sigma / mu の入力はこのモードでは
参照情報であり、指定Zsから物理的な表皮深さを逆算しません。出力の
`skin_depth_wp_mm` は null、`impedance_model` は `specified-panel` です。
`writeIHPanelImpedance` は全要素の順序付き複素ベクトルを受け取ります。
座標・材料領域からの要素選択は呼び出し側で `centroids_m` を使って行います。
