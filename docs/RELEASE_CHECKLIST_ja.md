# bivsurv-luminosity: 初回公開の手順

## 確定した事項

新規リポジトリ名は`bivsurv-luminosity`, ソフトウェアライセンスは`BSD-3-Clause`, 著作権者表示は`Tsutomu T. Takeuchi`です. 最初の公開タグは`v0.1.0rc1`とし, GitHubではPrereleaseとして扱います. 旧Privateリポジトリの履歴・可視性・ファイルは変更しません.

この準備版の作成時には, リモート作成, push, GitHub Actionsの実行, GitHub Releaseの発行, Zenodo登録のいずれも行っていません. リポジトリの所有者とURLは, 実際のGitHub接続または作成画面で確認します.

## 1. 配布物を確認する

ZIPを新しいフォルダへ展開し, `LICENSE`と`pyproject.toml`が直接見える`bivsurv-luminosity`フォルダへ移動します. ZIPファイルを一個だけリポジトリへ置くのではなく, 展開したソースを公開します.

```bash
python -m pip install -e ".[test]"
python scripts/check_release.py
python -m pytest -q
```

保存済みの数値結果を最初から再実行する必要はありません. 全実験を再生成する場合はREADMEに記載した別出力先を用いてください. 主ライブラリはPython 3.10以上を宣言し, 本準備版はPython 3.13.5で検査します. `requirements-pinned.txt`は既存実験の保存環境であり, 他のPython環境へそのまま適合するとは限りません.

## 2. 新しいPublicリポジトリを作る

GitHubの新規作成画面で正しい所有者を確認し, 名前を`bivsurv-luminosity`, 可視性をPublicにします. 既に同名が存在する場合は上書きせず, 内容を確認して作業を止めます. Descriptionには次を使えます.

```text
D-completion bivariate survival estimation for two-band luminosity catalogues, with reproducible benchmarks and validation code.
```

README, .gitignore, LICENSEはこの一式に含まれるため, 作成画面では追加しません. GitHub DesktopやGitを使ってフォルダ全体をコミット・pushします. 数千ファイルがあるので, ブラウザで一度に全てドラッグする方法には依存しないでください.

GitHub CLIを使う場合の例です. 新しく展開したフォルダでのみ実行し, 正しいアカウントでログインしていることを先に確認してください. `gh repo create`が既存リポジトリなどの理由で失敗した場合は, 自動的に別のremoteへpushしないでください.

```bash
gh auth status
git init -b main
git add .
git commit -m "Prepare licensed reproducibility release 0.1.0rc1"
gh repo create bivsurv-luminosity --public --source=. --remote=origin --push --description "D-completion bivariate survival estimation for two-band luminosity catalogues, with reproducible benchmarks and validation code."
```

Gitの著者設定が未登録の場合は, 本人のGitHub設定に対応する氏名とメールをローカルに設定してからコミットします. 認証トークンや秘密鍵をファイルに貼り付けてコミットしないでください. この手順は新規リポジトリだけが対象で, 旧リポジトリのremoteは再利用しません.

## 3. 実際のURLと公開情報を揃える

リポジトリ作成後, 実在するURLを`CITATION.cff`の`repository-code`と必要に応じたREADMEリンクへ追加します. 公開前の説明も, 実際の公開状態に合わせて更新します. `date-released`は実際のリリース日が決まったときだけ追加し, DOIは発行されるまで記入しません. 引用情報へ架空の識別子を入れないでください.

文書を変更したら, 意図した差分を確認してチェックサムを再生成し, メタデータ更新をコミット・pushします.

```bash
python scripts/check_release.py --write
python scripts/check_release.py
python -m pytest -q
git add CITATION.cff README.md README_ja.md SHA256SUMS
git commit -m "Record verified publication metadata"
git push
```

## 4. CIを確認してから初回プレリリースを発行する

GitHubのActionsで`Tests and smoke benchmark`が成功したことを確認します. CIはPython 3.13上の単体テストと小標本の実験であり, 全480適合の再実行や新しい大域証明ではありません. 失敗した場合はリリース前に原因を確認します.

GitHubのReleasesで新タグ`v0.1.0rc1`, Target=`main`, Prereleaseを指定し, `RELEASE_NOTES.md`の内容を説明に使います. 公開済みの同名タグがあれば動かしたり削除したりせず, 次の版番号を使ってください. リリース作成前に, 必要なメタデータ更新が含まれるコミットを選びます.

CLIを用いる場合は, 検査済みコミットでローカルタグを作ってpushした後に発行します.

```bash
git tag -a v0.1.0rc1 -m "First public release candidate"
git push origin v0.1.0rc1
gh release create v0.1.0rc1 --verify-tag --prerelease --title "bivsurv-luminosity 0.1.0rc1" --notes-file RELEASE_NOTES.md
```

## 5. 論文へのコード公開情報

論文のCode availabilityは, リポジトリとリリースが実在することを確認した後に追加します. 永続アーカイブを利用する場合は, 実際に発行された版固有DOIを後から記載します. 現時点では, ソフトウェアが公開済みまたはDOI取得済みであるという文章を論文へ挿入しません.

## 参照実装についての留意点

`reference_experiments/stage22/resumption_checks/compare_runs.py`は旧作業領域を前提とする保存用の比較補助です. 新しい再現比較には`scripts/verify_visual.py`を使います. 高速APIは整合した有限記録境界に特化し, 任意の非整合限界や測定誤差を扱う汎用APIへ変更されたわけではありません. 新実験の数値停止診断と, 既存Stage 25の別13候補に対する大域証拠を区別して公開します.

## 公式手順

- リポジトリ作成: https://docs.github.com/en/repositories/creating-and-managing-repositories/creating-a-new-repository
- GitHub CLI作成: https://cli.github.com/manual/gh_repo_create
- 引用ファイル: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
- リリース作成: https://cli.github.com/manual/gh_release_create
- ライセンス本文: https://spdx.org/licenses/BSD-3-Clause.html
