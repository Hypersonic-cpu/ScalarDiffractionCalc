# ScalarDiffractionCalc
Scalar diffraction calc using Angular Spectrum Method.

---------

你需要在此目录下添加 `MklLoader.fsx` 文件, 内容为 
```fsharp
let MKLProviderPath = "/folder/of/MKL-Provider-lib"
```
如果不在 Intel 处理器上运行, 似乎不需要 MKL, 可以新建一个 branch 删掉这些.
