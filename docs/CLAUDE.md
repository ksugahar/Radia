## コーディング規約
- インデント: タブ文字
- スクリプト形式: PowerShell (.ps1) を使用（.cmd は古い形式）

## PyPI 配布に必要な重要ファイル - 削除禁止

以下のファイルは PyPI パッケージの配布に必要なため、削除してはいけません：

1. **Build_PyPI_Package.ps1** - PyPI パッケージビルドスクリプト
   - 旧: build_pypi_package.cmd（廃止）
   - source distribution (sdist) と wheel distribution (bdist_wheel) をビルド
   - 使用方法: `.\Build_PyPI_Package.ps1`

2. **setup.py** - パッケージセットアップ設定
   - パッケージメタデータと依存関係を定義
   - バイナリ拡張モジュール (.pyd) のパッケージング

3. **pyproject.toml** - モダンな Python プロジェクト設定
   - PEP 518/621 準拠
   - ビルドシステムとプロジェクトメタデータ

4. **MANIFEST.in** - 配布ファイルの包含ルール
   - ソース配布に含めるファイルを指定

5. **LICENSE** - ライセンステキスト
   - LGPL-2.1 + 元の RADIA BSD-style ライセンス

6. **COPYRIGHT.txt** - 元の Radia 著作権表示
   - ESRF (1997-2018) の著作権を維持
   - 絶対に削除しないこと

## ローカル開発ファイル（.gitignore に含む）

以下のファイルはローカル環境のみで使用し、リポジトリには含めません：

- **Publish_to_PyPI.ps1** - PyPI アップロードスクリプト（認証情報を含む可能性）
- **publish_to_pypi.cmd** - 旧形式（廃止予定）
- **CLAUDE.md** - プロジェクト固有の開発メモ

