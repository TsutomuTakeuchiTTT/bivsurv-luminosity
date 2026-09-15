# 第22段階: 連続親分布からの初期統合pilot

120標本の生成・点推定・親信頼集合の真値所属と,
代表3標本の多限界生存射影を保存しています.
第21段階のLaTeX記録も同梱します.

## 再実行

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit_stage22.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python project_stage22.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python refine_projections_stage22.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_stage22.py
```

Windowsでは同じ環境変数を利用するシェルの方法で設定してください.
Python, NumPy, SciPyを使います. ネット接続や外部データは不要です.
ソースのversionはoutputs/stage22_sampling_results.jsonに記録しています.
数値ソルバーの提案は環境によって変わり得ます.

- audit_stage22.py: 個体生成, 公開記録, 学習予測, 点推定, 評価.
- project_stage22.py: 無改変の第20段階エンジンによる初期証拠.
- refine_projections_stage22.py: 無改変の第21段階エンジンによる再検証.
- verify_stage22.py: 保存配列, 積分, 予測と証拠を読み取り専用で検査.
- outputs/observed_inputs: 推定側への公開入力.
- outputs/latent_evaluation_only: Dを含む評価専用の潜在標本.
- outputs/fits: 生のEM履歴, 予測, 真値所属の保証区間.
- outputs/refined_projection_certificates: 代表3標本の最終証拠.

120標本の親集合所属と,3標本だけの鋭い射影を区別してください.
親の不可視質量は元のEMでは初期値を保ちます.
数値停止を大域MLEの保証へ読み替えていません.
