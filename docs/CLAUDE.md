## コーディング規約
- インデント: タブ文字
- スクリプト形式: PowerShell (.ps1) を使用（.cmd は古い形式）

## PyPI 配布に必要な重要ファイル - 削除禁止

以下のファイルは PyPI パッケージの配布に必要なため、削除してはいけません：

1. **setup.py** - パッケージセットアップ設定
   - パッケージメタデータと依存関係を定義
   - バイナリ拡張モジュール (.pyd) のパッケージング

2. **pyproject.toml** - モダンな Python プロジェクト設定
   - PEP 518/621 準拠
   - ビルドシステムとプロジェクトメタデータ

3. **MANIFEST.in** - 配布ファイルの包含ルール
   - ソース配布に含めるファイルを指定

4. **LICENSE** - ライセンステキスト
   - LGPL-2.1 + 元の RADIA BSD-style ライセンス

5. **COPYRIGHT.txt** - 元の Radia 著作権表示
   - ESRF (1997-2018) の著作権を維持
   - 絶対に削除しないこと

## ローカル開発ファイル（.gitignore に含む）

以下のファイルはローカル環境のみで使用し、リポジトリには含めません：

- **Publish_to_PyPI.ps1** - PyPI ビルド・アップロード統合スクリプト（ローカルのみ）
  - 旧: build_pypi_package.cmd, publish_to_pypi.cmd（廃止）
  - ビルドとアップロードを統合した単一スクリプト
  - .gitignore に含まれているため Git にコミットされない
  - 使用方法:
    ```powershell
    $env:PYPI_TOKEN = "your-token-here"
    .\Publish_to_PyPI.ps1
    ```
  - **セキュリティ**: API トークンは環境変数として設定（スクリプトにハードコードしない）
- **Build_PyPI_Package.ps1** - 廃止（Publish_to_PyPI.ps1 に統合）
- **CLAUDE.md** - プロジェクト固有の開発メモ

## PyPI パッケージ公開ワークフロー

1. **ビルド**: `.\Build.ps1` でコアモジュールをビルド（必須）
2. **公開**: API トークンを設定して `.\Publish_to_PyPI.ps1` を実行（ローカルのみ）
   ```powershell
   $env:PYPI_TOKEN = "your-token-here"  # PyPI API トークンを設定
   .\Publish_to_PyPI.ps1                 # ビルドとアップロードを統合して実行
   ```

