# Reference Standard Tables

ここには、標準星の reference flux table を置きます。

想定フォーマットは 2 列以上の plain text です。

```text
# wavelength_A  flux_erg_s_cm2_A
4000.0  1.23e-13
4010.0  1.22e-13
...
```

波長は Angstrom 単位です。flux は `erg s^-1 cm^-2 A^-1` のような物理 flux density を想定しています。

placeholder や手作りのファイルを、校正済み reference として扱わないでください。science 用の強度校正には、実測された spectrophotometric standard table を使ってください。
