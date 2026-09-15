# 第25段階: 大域尤度ギャップの証拠

`THEORY_ja.md` が導出と適用条件, `SUMMARY_ja.md` が結果.
`stage24_record.tex` は前段階の追加用LaTeX記録.

## 再実行

Python 3とNumPy/SciPyを使用する. 入力は全て同梱済み.
次をこのディレクトリで実行する.

```bash
OPENBLAS_NUM_THREADS=1 python global_bound.py --tol 1/1000000 --max-nodes 501 --outdir outputs/main_certificates
OPENBLAS_NUM_THREADS=1 python global_bound.py --tol 1/1000000 --max-nodes 501 --outdir outputs/rerun_certificates
OPENBLAS_NUM_THREADS=1 python global_bound.py --case coarse_d3__lambda0_N400_r004 --tol 1/1000000 --max-nodes 1 --outdir outputs/budget_control
OPENBLAS_NUM_THREADS=1 python global_bound.py --case coarse_d3__lambda0_N400_r004 --tol 1/100000000 --max-nodes 501 --outdir outputs/strict_control
OPENBLAS_NUM_THREADS=1 python verify_stage25.py
OPENBLAS_NUM_THREADS=1 python audit_stage25.py
```

主ケースは平均対数尤度ギャップ1e-6を全13件で達成.
予算対照と高精度対照は未達成が正しい保存結果である.
公開入力・元候補のスナップショットは `inputs/`.
祖先コード・元入力・元証拠は `reference_stage24/` に保存した.
新しい探索で元候補質量を選び直すことはない.

## 主な出力

- `outputs/main_certificates/`: 全13件の大域上界と被覆木.
- `outputs/global_gap_comparison.csv`: 一階残差, 大域ギャップ, 80GEM後のギャップ.
- `outputs/independent_verification.json`: 最適化を呼ばない別表示の検証.
- `outputs/audit_results.json`: 全件再実行の一致と算術対照.
- `outputs/same_likelihood_family.json`: 同じ尤度を持つ三つの親分布.
- `outputs/stage24_recheck.json`: 前段階の1040更新の再検査結果.

有理数対数を扱う `reference_stage24/vendor/exact.py` と観測幾何コードは既存のまま使用.
独立検証では新探索の `tangent_coeff`, `polytope`, LPやSLSQPを呼ばない.
ただし共有する有理数対数実装と観測幾何の定義まで独立に再発明したわけではない.
入力や証拠を変更すると検証は失敗する. 厳密なMLE認定や統計的な1e-6誤差ではない.
