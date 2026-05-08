# 実験用 強度校正

このディレクトリは、安定版の quicklook pipeline とは意図的に分けています。
`main.py` や `pipeline/` パッケージを壊さずに、強度校正まわりを試すための場所です。

## 現在の対象

入力には quicklook pipeline が出力した 1D テキストスペクトルを使います。

```text
quicklook/{object}_1d.txt
```

必要なものは以下です。

- quicklook pipeline で抽出した標準星スペクトル
- 標準星の reference flux table、または built-in の近似テンプレート
- 必要なら、同じ quicklook pipeline で抽出した science target のスペクトル

標準星の実測カウントと reference/template flux を比較して、感度曲線を作ります。

```text
sensitivity(lambda) = reference_flux(lambda) / standard_counts_rate(lambda)
target_flux(lambda) = target_counts_rate(lambda) * sensitivity(lambda)
```

target を指定した場合は、デフォルトで FITS header の `AIRMASS` も読み、
標準星と target の airmass 差による大気減光の差を補正します。

```text
target_flux_corrected(lambda)
  = target_flux(lambda) * 10^(0.4 * k(lambda) * (X_target - X_standard))
```

ここで `k(lambda)` は mag/airmass の大気減光係数、`X` は airmass です。
現在の `k(lambda)` は内蔵の簡易 optical extinction curve です。観測所の実測値を使う場合は
`--extinction-curve` で 2 列の table を指定できます。

デフォルトでは、名目上の有効波長域を 4500-8500 Å としています。この範囲外は感度曲線の fitting から外し、図では薄い斜線つきの領域として表示します。

これはまだ production reduction ではありません。slit loss、時変の大気透過、
差分大気分散、
order overlap、telluric correction などはモデル化していません。

## 標準星 Reference の形式

reference file は 2 列以上の plain text を想定しています。

```text
# wavelength_A  flux_erg_s_cm2_A
4000.0  1.23e-13
4010.0  1.22e-13
...
```

波長は Angstrom 単位です。flux は `erg s^-1 cm^-2 A^-1` のような物理 flux density を想定しています。

## 使い方

### HR4468 の built-in 黒体テンプレートで感度曲線を作る

```bash
python3 experimental_flux_calibration/flux_calibrate.py \
  --standard quicklook/HR4468_1d.txt \
  --output-dir flux_calibrated \
  --save-template
```

この実行では、HR4468 の実測 quicklook spectrum と黒体テンプレートを比較し、感度曲線を保存します。
露出時間は、デフォルトでは `spectra/` 内の FITS header から `EXPTIME` を読んで使います。
target まで処理する場合は、同じ FITS header から `AIRMASS` も読み、差分 airmass 補正を行います。

### 同じ実行で target に感度曲線を適用する

```bash
python3 experimental_flux_calibration/flux_calibrate.py \
  --standard quicklook/HR4468_1d.txt \
  --target quicklook/SN2026acd_1d.txt \
  --output-dir flux_calibrated
```

line ID も表示する場合:

```bash
python3 experimental_flux_calibration/flux_calibrate.py \
  --standard quicklook/HR4468_1d.txt \
  --target quicklook/SN2026acd_1d.txt \
  --output-dir flux_calibrated \
  --line-z 0.01 \
  --line-list H_Balmer He_I He_II C_III C_IV N_III N_IV N_V
```

target の flux calibrated plot には、大気吸収帯も陰影で表示されます。

有効波長域を変更する場合:

```bash
python experimental_flux_calibration/flux_calibrate.py \
  --standard quicklook/HR4468_1d.txt \
  --output-dir flux_calibrated \
  --valid-wave-min 4500 \
  --valid-wave-max 8500
```

FITS file を明示する場合:

```bash
python3 experimental_flux_calibration/flux_calibrate.py \
  --standard quicklook/HR4468_1d.txt \
  --standard-fits spectra/spec260212-0289.fits spectra/spec260212-0292.fits \
  --target quicklook/SN2026acd_1d.txt \
  --target-fits spectra/spec260212-0287_copy.fits spectra/spec260212-0288_copy.fits \
  --output-dir flux_calibrated
```

### 外部 reference table を使う

```bash
python3 experimental_flux_calibration/flux_calibrate.py \
  --standard quicklook/HR4468_1d.txt \
  --reference experimental_flux_calibration/reference_standards/HR4468_reference.txt \
  --target quicklook/SN2026dix_1d.txt \
  --output-dir flux_calibrated
```

## 出力

```text
flux_calibrated/{standard}_sensitivity.txt
flux_calibrated/{standard}_sensitivity.png
flux_calibrated/{standard}_standard_vs_template.png
flux_calibrated/{target}_flux_calibrated.txt       # --target 指定時のみ
flux_calibrated/{target}_flux_calibrated.png       # --target 指定時のみ
flux_calibrated/{target}_airmass_correction.txt    # --target 指定時のみ
flux_calibrated/{target}_airmass_correction.png    # --target 指定時のみ
```

## 注意

- 露出時間は FITS header の `EXPTIME` をデフォルトで使います。quicklook は A-B 差分の正像側を抽出しているため、A/B の 2 枚が見つかった場合は正像側の A frame の `EXPTIME` を使います。
- 露出時間を手で上書きしたい場合は `--standard-exptime` / `--target-exptime` を指定してください。
- airmass も同じく正像側の A frame の `AIRMASS` を使います。手で上書きする場合は `--standard-airmass` / `--target-airmass` を指定してください。
- 差分 airmass 補正を切りたい場合は `--no-airmass-correction` を指定してください。
- 内蔵の大気減光係数は簡易値です。観測所ごとの extinction curve がある場合は `--extinction-curve wavelength_A k_mag_per_airmass` 形式の table を指定してください。
- 感度曲線の fitting/smoothing では、telluric band をデフォルトで mask します。
- 標準星の Balmer absorption 付近もデフォルトで mask します。
- 標準星カウントが低すぎる波長端は、感度曲線が暴れやすいのでデフォルトで mask します。
- デフォルトの有効波長域は 4500-8500 Å です。範囲外は図で薄い斜線つきになります。
- smoothing は単純な median smoothing です。`--smooth-window` を変えながら `*_sensitivity.png` を確認してください。
- built-in の HR4468 template は、V 等級で規格化した滑らかな黒体近似です。実験には便利ですが、実測の spectrophotometric standard table の代わりにはなりません。
