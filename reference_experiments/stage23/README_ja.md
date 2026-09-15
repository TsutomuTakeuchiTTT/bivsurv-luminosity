# 第23段階: 保存母標本の観測深度・記録解像度の対比較

Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0, SymPyを用いた保存実行です.
入力は第22段階の120母標本をそのままコピーしています. 新しい乱数は生成しません.
全データの再観測結果は720条件ですが, 独立標本数は120です.

```
python depth_study.py --outdir outputs
python population_fibers.py --outdir outputs
python audit_controls.py --outdir outputs
python verify_stage23.py --outdir outputs
```

`depth_study.py` は学習予測・点推定・真値所属を計算します.
`population_fibers.py` は別の母集団識別問題を解き, 有理数の主双対証拠を保存します.
`verify_stage23.py` はLP・非線型最適化・乱数生成・数値EM再適合を行わずに検査します.
学習2EMは定義から再計算します.

- `outputs/summary.csv`: 36設定の各20反復の集約値.
- `outputs/trials.csv`: 720観測設計評価と2160参考生存値.
- `outputs/population_ranges.csv`: 54個の母集団識別範囲.
- `outputs/population_certificates/`: 108端点親分布と全列双対証拠.
- `outputs/independent_verification.json`: 独立検証結果.
- `outputs/structural_controls.json`: 記録情報の保持/喪失と尾部制約違反.
- `outputs/reproducibility.json`: 新しい二度の実行の比較.
- `THEORY_ja.md`: 観測設計, 初期測度, 主結果の証明, 適用範囲.
- `stage22_record.tex`: 前段階の追記用LaTeX. 環境・参照を静的検査; コンパイル未実施.

真の支持・密度・D個数・親総数は評価に使いますが, 条件付き尤度の適合には使いません.
原EMは尾部制約を強制する更新ではありません. クラス外の13候補を修復せずに記録しています.
719回の真値採択は有限反復の実現結果であり, 被覆確率の厳密値ではありません.
母集団識別範囲は信頼区間ではありません. 720件の鋭い標本信頼射影は今回追加していません.
