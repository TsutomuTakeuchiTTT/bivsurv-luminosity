# Stage24再現パッケージ

Stage23で明示尾部クラスを外れた13設計・標本ケースを扱います.
これらは8種類の既存母標本に由来します. 新しい乱数生成はありません.

## 実行

Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0の環境で検証しました.
標準ライブラリのfractions.Fractionとdecimalも使います.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python audit_stage24.py --outdir new_outputs
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_stage24.py --outdir new_outputs
```

既定は各ケース80回の有理数検証付きGEMです.
原結果は`outputs/`, 再実行の一致検査は`outputs/reproducibility.json`です.
`outputs/certificates/`には各ケースの80更新と, 別経路の境界候補を保存します.

`verify_stage24.py`は最適化を呼びません.
記録セルの構築と学習2EMの定義は再計算します.
M-stepの提案関数は呼ばず, KKT乗数から解を再構成して検査します.

## 出力の区別

- `GEM_steps`: 可行な初期測度からの80回の増加証明. 固定80回での収束は主張しません.
- `independent_candidate`: 同じ制約付き尤度を別に最適化した候補.
  可行性・80更新後からの尤度増加・一次方向残差上界を確認します.
  この候補をGEMの一更新と呼んだり, 大域MLEと認定したりはしません.
- `comparison.csv`: 旧クラス外候補, 80GEM候補, 別経路の制約付き候補の診断.

観測入力と旧適合は`inputs/`にそのまま保存し,
源ファイルのハッシュは`PROVENANCE.json`に記録しています.
`vendor/`の三ソースは以前の検証パッケージから無改変で再利用しています.
全段階を呼ぶ依存関係は不要で, このフォルダだけで上の再実行ができます.

`stage23_record.tex`は今回の追記用LaTeX本文です.
原稿全体やBibTeXは変更していません.
