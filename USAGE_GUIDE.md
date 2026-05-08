# KOspec Quicklook Pipeline - 詳細ガイド

## 概要

このパイプラインは、分光観測直後にraw FITSファイルから2D spectrum画像と1D spectrumを自動生成し、観測が成功したかを迅速に確認するためのツールです。

## システム要件

- Python 3.7 以上
- 依存パッケージ: numpy, scipy, astropy, matplotlib

## インストール

```bash
# プロジェクトディレクトリに移動
cd KOspec_pipeline

# 依存パッケージをインストール
pip install -r requirements.txt
```

## 基本的な使用方法

### 最小構成での実行

```bash
python3 main.py
```

このコマンドは以下を実行します：
1. `./spectra/` ディレクトリから FITS ペアを検索
2. A-B差分を計算
3. 2D spectrum と 1D spectrum を生成
4. 結果を `./quicklook/` に保存

### よく使うオプション

#### 赤方偏移を指定して実行

```bash
python3 main.py --z 0.05
```

指定した赤方偏移に対応した主要輝線の位置を1Dスペクトルプロットに縦線で表示します。

#### Aperture幅を変更

```bash
python3 main.py --aperture 15
```

デフォルトは10ピクセルです。S/N比が低い場合は幅を広げることで改善される可能性があります。

#### 詳細なログを表示

```bash
python3 main.py -v
```

処理の詳細を確認できます。トラブルシューティング時に便利です。

#### ディレクトリを指定

```bash
python3 main.py --spectra-dir /data/my_observation --output-dir /data/quicklook_results
```

### すべてのオプション

```bash
python3 main.py --help
```

## ファイル構造

### 入力ファイル形式

入力FITSファイルは以下の命名規則に従う必要があります：

```
./spectra/
  ├── {object_name}_A.fits    # A-position frame
  ├── {object_name}_B.fits    # B-position frame
  └── ...
```

A と B のペアが自動的にマッチングされます。

### 出力ファイル

処理結果は `./quicklook/` に以下の形式で保存されます：

```
./quicklook/
  ├── {object}_2d.png         # 2D spectrum画像 (trace とaperture表示)
  ├── {object}_1d.txt         # 1D spectrum (wavelength, flux, sky)
  ├── {object}_1d.png         # 1D spectrum プロット
  └── ...
```

#### ファイルの説明

- **`{object}_2d.png`**: 2D spectrumの画像表示
  - 背景：スペクトル強度（RdBu_r カラーマップ、対称スケール）
  - 黄色の破線：Trace位置
  - 緑色の矩形：抽出用Aperture

- **`{object}_1d.txt`**: テキスト形式の1D spectrum
  - 3列：Wavelength(Å), Flux(ADU), Sky(ADU)
  - ヘッダーに処理パラメータを記載

- **`{object}_1d.png`**: 1D spectrum のプロット
  - 青線：抽出スペクトル
  - 赤破線：Background成分
  - 縦線：輝線位置（赤方偏移指定時）

## 処理の詳細

### A-B差分

各オブジェクトフレームについて：

```
Processed_Image = A_Position - B_Position
```

これにより、dark current と sky background をキャンセルします。

### Trace推定

1. Positive部分のみを取得（A-Bから出現するネガティブスペクトルを除外）
2. 空間方向に統合して1Dプロファイルを作成
3. Gaussian fitting で trace 中心位置を推定

手動で trace 位置を修正することはできないため、信号が弱い場合は実行失敗する可能性があります。

### Aperture Extraction

Trace位置から固定幅の aperture で スペクトルを抽出します：

```
Spectrum_1D = Sum(Image_2D[trace-half_width : trace+half_width, :])
```

デフォルト aperture 幅：10 ピクセル

### Sky Subtraction

Aperture の上下（sky_offset ピクセル離れた場所）から sky level を推定し差分します。

```
Final_Spectrum = Spectrum_1D - Sky_Level
```

## エラーハンドリング

パイプラインは以下のエラーに対してロバストです：

### 処理を続行するエラー

- **壊れたFITSファイル**：該当オブジェクトをスキップ
- **Truncated file**：該当オブジェクトをスキップ
- **Header error**：該当オブジェクトをスキップ
- **Shape mismatch (A-B)** ：該当オブジェクトをスキップ

### 段階的な失敗処理

- **A-B差分失敗** → オブジェクト全体をスキップ
- **2D plot失敗** → 1D処理に進む（2D PNG は生成されない）
- **Trace推定失敗** → 1D extraction をスキップ、1D PNG は生成されない
- **1D extraction失敗** → 1D PNG は生成されない、2D PNG のみ保存

## 波長較正

現在のバージョンでは、Angstrom単位に統一した6次多項式の波長解法を使用しています：

```
Wavelength = c0 + c1*x + c2*x^2 + c3*x^3 + c4*x^4 + c5*x^5 + c6*x^6
```

デフォルト：
- c0 = 10439.47468362323
- c1 = -24.91445451495574
- c2 = 6.468775026217139E-2
- c3 = -1.183580177752003E-4
- c4 = 1.365974323268898E-7
- c5 = -8.73649357565434E-11
- c6 = 2.335710049781674E-14

### 将来の拡張

- Arc line fitting による精密波長較正
- 高次多項式解法
- 設定ファイルによる波長解法の指定

## 赤方偏移付き出力

`--z` オプションで赤方偏移を指定すると、観測波長系での主要輝線位置が1Dプロットに表示されます。

### 表示される輝線

デフォルト：
- H-alpha (6563 Å)
- H-beta (4861 Å)
- H Balmer series (3835, 3889, 3970, 4102, 4340, 4861, 6563 Å)
- He I optical lines (3889, 4026, 4471, 4713, 4922, 5016, 5876, 6678, 7065 Å)
- He II optical lines (4200, 4542, 4686, 5411, 6560 Å)
- Na I D (5890 Å)
- [O I] 6300 (6300 Å)
- Ca II NIR triplet (8498, 8542, 8662 Å)

### インタラクティブな輝線選択

quicklook が出力した 1D テキストスペクトルに対して、別スクリプトで表示する輝線を切り替えられます。

```bash
python interactive_lines.py quicklook/SN2026acd_1d.txt --z 0.01
python interactive_lines.py quicklook/SN2026acd_1d.txt --lines H_Balmer He_I He_II
```

通常の quicklook / `main_all.py` でも表示する line ID を指定できます。

```bash
python main_all.py --z 0.01 --line-list H_Balmer He_I He_II C_III C_IV N_III N_IV N_V
```

`main.py` / `main_all.py` はデフォルトで `spectra/` 内の全 A-B ペアを処理します。特定天体だけ処理する場合は `--objects` を指定します。

```bash
python main.py --objects SN2026acd
python main_all.py --objects HR4468 SN2026acd
```

### カスタム輝線

カスタム輝線リストを作成する場合は、`calibration.py` の `EMISSION_LINES` 辞書を編集してください。

## トラブルシューティング

### 問題：処理が途中で止まる

**原因の確認**
```bash
python3 main.py -v 2>&1 | tail -50
```

詳細なエラーメッセージを確認します。

### 問題：Trace推定が失敗する

可能な原因：
1. 信号が弱い（S/N比が低い）
2. Slit misalignment
3. オブザーバー位置が大きくずれている

対策：
- Aperture 幅を広げてS/Nを改善：`--aperture 15`
- 観測パラメータを確認

### 問題：1D spectrumに期待した輝線が見えない

原因：
1. 赤方偏移を誤って指定している
2. 波長範囲外
3. 信号が弱い

確認方法：
```bash
# 赤方偏移を確認
python3 main.py --z <correct_z_value>
```

### 問題：Sky subtraction後にノーマ値になっている

原因：
1. Aperture width が小さすぎる
2. Sky領域の見積もりが不正確

対策：
- Aperture 幅を広げる
- Sky_offset パラメータを調整（コード内）

## 出力例の解釈

### 2D Spectrum PNG

- 色が濃い部分：スペクトル強度が高い
- 赤い部分：ネガティブ値（A-Bの差分）
- 青い部分：ポジティブ値

### 1D Spectrum PNG

- 急激な上昇：輝線の検出
- 滑らかな曲線：連続光の検出
- ノーマルな平坦部分：ノイズまたはsky

## パイプライン拡張への準備

現在のコード構造は以下の拡張に対応しています：

1. **Arc line fitting** (`calibration.py` 拡張)
   - 参照スペクトル追加
   - 2次元波長マッピング

2. **Standard star flux calibration**
   - Sensitivity curve 導入
   - Magnitude calibration

3. **Summary HTML** (`visualization.py` 拡張)
   - 一括処理結果表示
   - 統計情報

4. **Config YAML** (config.yaml)
   - パラメータの外部管理
   - 観測サイト別設定

## パフォーマンス

- 処理時間：オブジェクトあたり ～1秒
- メモリ使用量：1枚のFITSあたり ～10 MB
- 大量処理：100枚程度なら ～2-3分

## ライセンス

MIT License

## 更新履歴

### v0.1.0 (Initial Release)
- A-B差分処理
- 2D spectrum 出力
- Trace推定と Aperture extraction
- 1D spectrum テキスト/PNG出力
- 赤方偏移対応
- ロバストなエラーハンドリング
