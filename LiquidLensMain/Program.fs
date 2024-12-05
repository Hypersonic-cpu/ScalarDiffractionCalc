open MathNet.Numerics
open MathNet.Numerics.LinearAlgebra
open MathNet.Numerics.IntegralTransforms
open Plotly.NET
open Plotly.NET.LayoutObjects

open LiquidLens.Optics

// Control.NativeProviderPath <- @"C:\Users\74674\.nuget\packages\mathnet.numerics.mkl.win-x64\3.0.0\runtimes\win-x64\native\"
// Control.UseNativeMKL()
let pi = System.Math.PI
let J_ = complex 0. 1.
let waveLength = 595e-9             // SI =595 nm.
let waveNumber = 2.*pi/waveLength
// Delta f * N = 1. / Delat X = const.
let freqRange = 2.0 / waveLength
let sptlRange = 3.0e-3              // SI
printfn $"{Optical.WaveLeastSampleNumber waveLength sptlRange}"
let N = stdin.ReadLine() |> int |> fun x -> (2*x+1)/2

//stdin.ReadLine() |> ignore
 
(* Plot *)
let qResize (mat:Matrix<complex>) =
    let N = mat.ColumnCount
    Array2D.init (N/8-1) (N/8-1) (fun i j -> mat.[N/2 + (i-N/16) * 8, N/2 + (j-N/16) * 8])
    |> Matrix.Build.DenseOfArray

let plotPhase (name:string) (mat:Matrix<complex>) =
    Chart.Heatmap(
        mat 
            |> Matrix.map(fun x -> x.Phase) 
            |> (fun x -> Matrix.toRowArrays x)
    )
    |> Chart.withTitle (name+" phase")
    |> Chart.show
let plotMag (name:string) (mat:Matrix<complex>) =
    Chart.Heatmap(
        mat 
            |> Matrix.map(fun x -> x.Magnitude) 
            |> (fun x -> Matrix.toRowArrays x)
    )
    |> Chart.withTitle (name+" mag")
    |> Chart.show

(* Generate original wave *)
let waveRadius = 4e-3
let f x y = 
    let rho = sqrt (waveRadius**2. + x**2. + y**2.)
    if x**2. + y**2. <= 1.5e-3**2. then 
        exp(J_*waveNumber*rho)
            * (-waveRadius)/rho * (1. + 1./J_/waveNumber/rho) / J_ / waveLength / rho
            * exp(J_ * 0.5e5 * x)
    else complex 0. 0.

//let f x y =
//    // if x**2. + y**2. < 1e-9 then complex 1e6 0.
//    if abs y < 0.1e-3 && abs x < 0.1e-3 + 0.6*y then complex 1e6 0.
//    else complex 0. 0.
let spatialX = [| for i in 0..N-1 -> (float i / (float N-1.)-0.5) * sptlRange|]
let waveFunc =
    Array2D.init N N 
        // x = i / (N-1) * W - W/2
        (fun i j ->
            (spatialX.[i], spatialX.[j]) ||> f
        )
    |> Matrix.Build.DenseOfArray

waveFunc |> qResize |> plotPhase "Original"
//stdin.ReadLine() |> ignore
//Optical.TransferWaveResizeInPlace sptlRange (0.20 * sptlRange) waveLength (-0.92 * waveRadius) waveFunc
//waveFunc |> qResize |> plotMag "92%"
//Optical.TransferWaveResizeInPlace (0.20*sptlRange) (0.04*sptlRange) waveLength (-0.08 * waveRadius) waveFunc
//waveFunc |> qResize |> plotMag "Final"
