# KOspec Quicklook Pipeline - 実装仕様書

## プロジェクト概要

Kagoshima Optical spectroglaph (KOspec) 用のquicklook pipelineシステムです。観測直後にraw FITSから2D spectrum画像と簡易1D spectrumを自動生成し、スペクトルが取れているか、輝線が見えているかを迅速に確認するツールです。

## 実装内容

### v0.1.0 完了機能

#### 1. 前処理モジュール (`preprocessing.py`)
- ✓ A-Bポジション画像の自動ペアリング
- ✓ A-B差分によるdark/sky相殺
- ✓ ファイルマッチングの自動化
- ✓ エラーハンドリング機構

#### 2. FITS読み込みモジュール (`loader.py`)
- ✓ 標準FITSファイルの読み込み
- ✓ 複数HDUへの対応
- ✓ データ型チェック
- ✓ 壊れたFITS/truncated fileの検出
- ✓ Header error の検出

#### 3. スペクトル抽出モジュール (`extraction.py`)
- ✓ 空間方向profileからのtrace推定
- ✓ 複数trace推定手法（peak, center_of_mass, gaussian_fit）
- ✓ 固定幅apertureでの1D抽出
- ✓ Sky subtraction (aperture周辺の値から推定)
- ✓ 処理情報の自動記録

#### 4. 波長較正モジュール (`calibration.py`)
- ✓ 線形波長解法
- ✓ 多項式波長解法の枠組み
- ✓ 主要輝線データベース（17種類以上の輝線）
- ✓ 赤方偏移対応の輝線変換
- ✓ デフォルト波長解の自動生成

#### 5. 可視化モジュール (`visualization.py`)
- ✓ 2D spectrum PNG出力（colormap, auto-scale, trace/aperture overlay）
- ✓ 1D spectrum PNG出力（格子, legend）
- ✓ 1D spectrum テキスト出力（wavelength/flux/sky）
- ✓ 赤方偏移での輝線マーキング
- ✓ 自動ファイル管理

#### 6. ユーティリティモジュール (`utils.py`)
- ✓ Trace推定アルゴリズム
- ✓ Aperture extraction 関数
- ✓ Sky level推定関数
- ✓ Sky subtraction 関数
- ✓ ロギング設定

#### 7. メインスクリプト (`main.py`)
- ✓ パイプラインオーケストレーション
- ✓ コマンドラインインターフェース
- ✓ エラーハンドリング（失敗オブジェクトのスキップ）
- ✓ 処理結果の詳細サマリー出力
- ✓ 段階的な処理フロー

### エラーハンドリング

#### 処理を続行するエラー

| エラー | 処理 | 結果 |
|-------|------|------|
| 壊れたFITS | スキップ | - |
| Truncated file | スキップ | - |
| Header error | スキップ | - |
| A-B shape mismatch | スキップ | - |
| Trace推定失敗 | 1D処理をスキップ | 2D PNG のみ |
| 1D extraction失敗 | 1D処理をスキップ | 2D PNG のみ |

全てのエラーでパイプラインが継続します。

### コマンドラインオプション

```
主要オプション：
  --spectra-dir     入力FITSディレクトリ（デフォルト: ./spectra）
  --output-dir      出力ディレクトリ（デフォルト: ./quicklook）
  --pattern-a       Aフレームパターン（デフォルト: _A.fits）
  --pattern-b       Bフレームパターン（デフォルト: _B.fits）
  --aperture        Aperture幅（ピクセル、デフォルト: 10）
  --z               赤方偏移（デフォルト: 0）
  -v, --verbose     詳細ログ
  -h, --help        ヘルプ表示
```

### 波長範囲対応の輝線

デフォルトで表示される輝線（赤方偏移対応）：
- H-alpha (6563 Å)
- H-beta (4861 Å)
- He I 5875 (5876 Å)
- Na I D (5890 Å)
- [O I] 6300 (6300 Å)
- Ca II NIR triplet (8498, 8542, 8662 Å)

追加可能な輝線（calibration.py の EMISSION_LINES に全て定義）：
- Hgamma, Hdelta
- He I 6678, He II 4686
- Na I D1/D2 (個別)
- [O I] 6364, [O III] 5007/4959
- Ca II H/K, Ca II NIR (全3本)
- [S II] 6717/6731
- [N II] 6548/6584

## パイプライン処理フロー

```
Input FITS (A/B pairs)
       ↓
[FITS Loading] → Load A frame & Load B frame
       ↓
[A-B Difference] → Subtract dark/sky
       ↓
[2D Visualization] → Plot 2D spectrum with auto-scale
       ↓
[Trace Estimation] → Estimate spatial position via Gaussian fit
       ↓
[Aperture Extraction] → Extract 1D spectrum (fixed-width aperture)
       ↓
[Sky Subtraction] → Estimate & subtract sky level
       ↓
[Wavelength Calibration] → Apply polynomial wavelength solution
       ↓
[Redshift Correction] → Apply z if specified
       ↓
[1D Visualization] → Plot with emission lines (if z given)
       ↓
[Save Results] → Text + PNG output
       ↓
Output → quicklook/{object}_2d.png, _1d.txt, _1d.png
```

## パフォーマンス特性

- **処理速度**: ~1秒/オブジェクト
- **メモリ使用量**: ~10 MB/FITS
- **並列処理**: 現在は順序処理（将来の最適化対象）
- **大量処理**: 100オブジェクト ≈ 2-3分

## 出力ファイル仕様

### 2D Spectrum PNG

```
File: {object}_2d.png
Dimensions: Auto (aspect=auto)
Colors: RdBu_r (対称スケール)
Overlay: Trace位置（黄色破線）+ Aperture（緑色矩形）
Scale: 2%-98% percentile (自動スケール)
```

### 1D Spectrum Text

```
File: {object}_1d.txt
Format: 3列（カンマ区切りではなく空白区切り）
Header: 処理パラメータ（コメント行）
Columns: Wavelength(Å), Flux(ADU), Sky(ADU)
Precision: %.6e (科学記法)
```

### 1D Spectrum PNG

```
File: {object}_1d.png
Plot: Wavelength vs Flux
Blue line: 抽出スペクトル
Red dashed line: Sky成分（存在する場合）
Vertical lines: 輝線位置（赤方偏移対応）
Grid: Enabled with alpha=0.3
Legend: 自動生成
```

## 設定ファイル

`config.yaml` で以下の設定可能（将来の拡張用）：
- 入出力ディレクトリ
- Extraction パラメータ
- 波長較正タイプ
- Trace推定手法
- 表示輝線リスト

現在は参考用、実装は未完。

## 将来の拡張点

### 高優先度
- [ ] Arc line fitting による精密波長較正
- [ ] Config.yaml の完全統合
- [ ] Summary HTML 生成
- [ ] 処理ログの保存

### 中優先度
- [ ] Bad pixel mask 対応
- [ ] Standard star flux calibration
- [ ] Slit function fitting
- [ ] 2D wavelength mapping

### 低優先度
- [ ] Optimal extraction
- [ ] CCD cosmetic correction
- [ ] バッチ処理の最適化
- [ ] リアルタイム監視モード

## テスト状況

### 正常系テスト ✓
- [x] 標準的なFITS入力
- [x] 複数オブジェクト処理
- [x] 赤方偏移オプション
- [x] カスタムAperture幅

### 異常系テスト ✓
- [x] 壊れたFITSファイル
- [x] Truncated file
- [x] 空のファイル
- [x] 非FITS形式
- [x] Shape mismatch (A-B)

### エラーハンドリング ✓
- [x] Trace推定失敗
- [x] Shape mismatch での段階的失敗
- [x] パイプライン全体の継続性確認

## ファイル構成

```
KOspec_pipeline/
├── main.py                 # メインスクリプト
├── config.yaml             # 設定ファイル（テンプレート）
├── requirements.txt        # 依存パッケージ
├── README.md              # 基本README
├── USAGE_GUIDE.md         # 詳細ガイド
├── IMPLEMENTATION.md      # 実装仕様書（このファイル）
├── create_test_data.py    # テストデータ生成
├── create_broken_data.py  # 異常データ生成
│
├── pipeline/              # メインパッケージ
│   ├── __init__.py       # パッケージ初期化
│   ├── loader.py         # FITS読み込み
│   ├── preprocessing.py  # A-B差分処理
│   ├── extraction.py     # Trace推定・Aperture抽出
│   ├── calibration.py    # 波長較正
│   ├── visualization.py  # 可視化
│   └── utils.py          # ユーティリティ
│
├── spectra/              # 入力FITSディレクトリ
│   ├── {object}_A.fits
│   ├── {object}_B.fits
│   └── ...
│
└── quicklook/            # 出力ディレクトリ
    ├── {object}_2d.png
    ├── {object}_1d.txt
    ├── {object}_1d.png
    └── ...
```

## 依存パッケージ

| パッケージ | バージョン | 用途 |
|-----------|-----------|------|
| numpy | ≥1.19 | 数値計算 |
| scipy | ≥1.5 | 信号処理、最適化 |
| astropy | ≥4.0 | FITS I/O |
| matplotlib | ≥3.0 | 可視化 |

## 使用例

### 基本実行
```bash
python3 main.py
```

### 赤方偏移z=0.05で実行
```bash
python3 main.py --z 0.05
```

### Aperture幅15ピクセルで実行
```bash
python3 main.py --aperture 15 -v
```

### カスタムディレクトリ指定
```bash
python3 main.py --spectra-dir /data/obs2024 --output-dir /data/results
```

## ライセンス

MIT License

## 開発者向けメモ

### コード構造の特徴

1. **モジュール化**: 各機能が独立したモジュールで実装
2. **エラーハンドリング**: try-catch で例外を処理
3. **ロギング**: logging モジュール使用
4. **拡張性**: 新機能追加が容易な設計

### 今後の拡張時の注意点

- `calibration.py` の EMISSION_LINES に新しい線を追加可能
- `extraction.py` に新しいtrace推定手法を追加可能
- `visualization.py` に新しいプロット形式を追加可能
- Config体系の統合（YAMLパース実装）

### デバッグ方法

1. 詳細ログを有効化：`python3 main.py -v`
2. 特定オブジェクトの処理確認
3. 出力PNG/TXTの内容確認
4. エラーメッセージからボトルネック特定
