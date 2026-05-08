# KOspec Quicklook Pipeline

低分散分光器用のquicklookです。Aポジ-Bポジの差分を取り、2D spectrum画像と簡易1D spectrumを自動生成します。

## 機能

### 現在のバージョン (v1.0.0)

- **A-B差分**: Aポジション、Bポジションの2枚の画像の差分でdark/sky成分を相殺
- **2D spectrum可視化**: PNG出力で2D spectrum画像を確認
- **Trace推定**: 空間方向profileから自動的にtrace位置を推定
- **Aperture extraction**: 固定幅apertureを使用した1Dスペクトル抽出
- **Sky subtraction**: Aperture周辺領域からのsky推定と差分
- **1D spectrum出力**: テキスト形式とPNG両方に保存
- **波長較正**: Angstrom単位に統一した指定5次多項式を使用 (f(x) = 10314.3294 - 22.5327275*x + 4.81183127E-2*x^2 - 6.39126818E-5*x^3 + 4.59967054E-8*x^4 - 1.35717307E-11*x^5)
- **赤方偏移対応**: `--z` オプションで主要な輝線を赤方偏移させて表示
- **OBJECTキーワードベースのペアリング**: ファイル名ではなくFITSヘッダーのOBJECTキーワードから天体名を取得してA-Bペアを形成
- **ロバストなエラーハンドリング**: 壊れたFITSやエラーがあってもパイプライン全体は止まらない

## 使用方法

### 基本的な使い方

```bash
python main.py
```

デフォルトは `./spectra/` から入力を読み込み、`./quicklook/` に出力します。

### オプション

```bash
python main.py --help
```

主要なオプション：
- `--spectra-dir`: 入力FITS directory (デフォルト: ./spectra)
- `--output-dir`: 出力directory (デフォルト: ./quicklook)
- `--pattern-a`: Aフレームパターン (デフォルト: _A.fits)
- `--pattern-b`: Bフレームパターン (デフォルト: _B.fits)
- `--aperture`: Aperture幅（ピクセル、デフォルト: 10）
- `--z`: 赤方偏移（輝線表示用、デフォルト: 0）
- `-v, --verbose`: 詳細なログ出力

### 例

```bash
# 赤方偏移z=0.05で処理
python main.py --z 0.05

# Aperture幅を15ピクセルに変更
python main.py --aperture 15 -v
```

### Quicklook + 実験的強度校正

`main_all.py` は、安定版 quicklook の実行後に `experimental_flux_calibration/` の強度校正を続けて実行する wrapper です。

```bash
python main_all.py \
  --z 0.05 \
  --standard-object HR4468 \
  --flux-targets SN2026acd \
  --save-template \
  --smooth-window 75
```

デフォルトでは、標準星は `HR4468`、強度校正出力は `./flux_calibrated/` です。
露出時間は FITS header の `EXPTIME`、差分 airmass 補正は `AIRMASS` から読みます。
target には quicklook summary と同じ体裁で、強度校正済みスペクトルまで含めた `{object}_summary_all.png` も生成します。
airmass 補正を切る場合は `--no-airmass-correction`、観測所固有の大気減光曲線を使う場合は `--extinction-curve` を指定します。

## 入力ファイル形式

- **入力ディレクトリ**: `./spectra/`
- **ファイル名パターン**: `{object_name}_A.fits` と `{object_name}_B.fits`
- **形式**: 標準的なFITS 2D画像ファイル

例：
```
spectra/
  ├── star001_A.fits
  ├── star001_B.fits
  ├── star002_A.fits
  ├── star002_B.fits
  └── ...
```

## 出力ファイル

`./quicklook/` ディレクトリに以下が生成されます：

- `{object}_2d.png`: 2D spectrum画像（traceとaperture表示）
- `{object}_1d.txt`: 1D spectrum (テキスト形式, wavelength と flux)
- `{object}_1d.png`: 1D spectrum プロット（赤方偏移時は輝線表示）

## エラーハンドリング

パイプラインは以下のエラーに対してロバストです：

- 壊れたFITS file
- Truncated file
- Header error
- Shape mismatch (A-Bのサイズが異なる)
- Trace推定失敗
- 1D extraction失敗

これらのエラーが発生した場合、該当のオブジェクトをスキップしますが、パイプライン全体は継続します。
1D extraction に失敗しても 2D PNG は保存されます。

## 将来的な拡張

以下の機能に対応する構造になっています：

- [ ] `config.yaml`による設定ファイル管理
- [ ] Arc line fitting による正確な波長較正
- [ ] `summary.html` による処理結果の一括表示
- [ ] Standard star flux calibration
- [ ] Bad pixel mask対応
- [ ] Slit function fitting
- [ ] Optimal extraction

## トラブルシューティング

### 処理が途中で止まる場合
- 詳細なログを見る: `-v` オプションを使用
- ファイル形式を確認: FITSファイルが正しく読めるか確認

### Trace推定が失敗する場合
- S/N比が低い可能性
- Aperture幅を調整してみる

### 輝線表示がされない場合
- `--z` オプションでredshiftを正しく指定しているか確認
- 波長範囲が輝線の波長を含んでいるか確認

## インストール

```bash
# 必要なパッケージをインストール
pip install -r requirements.txt
```

## 著者

KOspec pipeline development team

## ライセンス

MIT License
